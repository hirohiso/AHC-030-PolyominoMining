import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Predicate;

@SuppressWarnings("unchecked")
public class Main {
    public static void main(String[] args) {
        solve(System.in, System.out);
    }

    private static void solve(PrintWriter pw, FastScanner fs) {
        //==================
        //初期構築フェーズ

        //初期値
        var N = fs.ni();
        var M = fs.ni();
        var e = Double.parseDouble(fs.n());
        var rand = new Random(0);

        var connector = new InteractiveConnector(pw, fs);

        var oilList = new LinkedList<OilFieldSet>();
        var max = 0;
        var minsize = Integer.MAX_VALUE;
        for (int i = 0; i < M; i++) {
            var d = fs.ni();
            if (d < minsize) {
                minsize = d;
            }
            var list = new LinkedList<Pair>();
            for (int j = 0; j < d; j++) {
                var x = fs.ni();
                var y = fs.ni();
                list.add(new Pair(x, y));
                max++;
            }
            oilList.add(new OilFieldSet(list));
        }
        //確定マス
        var result = new MiningResult(N);

        //推測フェーズ
        var macroGuesser = new MacroGuesser(N, oilList, e);
        macroGuesser.setup(connector);
        var cnt = 0;
        var flg = false;
        var greedyMode = false;
        LOOP:
        while (connector.counter < 2 * N * N - 1) {
            if (max - cnt <= 2 * minsize) {
                //全探索?
                if (macroGuesser.trySearch(connector, result)) {
                    System.out.println("おしまい");
                    return;
                }
            }
            macroGuesser.recal(result);
            //確度の高い答えっぽいものを見つけたら博打回答(打ちすぎないようにする)
            var maybeans = macroGuesser.pollMaybeAnswer();
            if (maybeans != null && connector.counter < N * N) {
                if (connector.answer(maybeans)) {
                    return;
                }
            }
            System.out.println(connector.counter);

            var list = macroGuesser.highPoint(max, result);

            for (Pair pair : list) {
                var v = connector.mine(pair);
                result.set(pair.a, pair.b, v);
                cnt += v;
                //油田を全部見つけたら終了
                if (cnt >= max) {
                    break LOOP;
                }
                if (v == 0) {
                    macroGuesser.setNoOil(pair.a, pair.b);
                } else {
                    macroGuesser.setOil(pair.a, pair.b);
                }
                //確率的に抜ける（エラー率が低いほど早く抜けて低い回数で生成回数を行う）
                if (rand.nextDouble() < 1 - e) {
                    break;
                }
            }

        }
        //回答フェーズ
        var list = result.resultList();
        connector.answer(list);
        /**/
    }

    //深さ優先探索をして油田の個数を見つけたか確認する
    private static LinkedList<LinkedList<Pair>> searchOil(MiningResult result, MacroGuesser macroGuesser, int islandSize) {
        var dif = new int[][]{{-1, 0}, {0, -1}, {1, 0}, {0, 1}};
        //発見している油田の数
        var islandCount = 0;
        var list = new LinkedList<LinkedList<Pair>>();
        var visited = new boolean[islandSize][islandSize];
        for (int i = 0; i < islandSize; i++) {
            for (int j = 0; j < islandSize; j++) {
                if (!result.isMined(i, j)) {
                    visited[i][j] = true;
                    continue;
                }
                if (result.isNoOilField(i, j)) {
                    visited[i][j] = true;
                    continue;
                }
                if (visited[i][j]) {
                    continue;
                }
                var stack = new LinkedList<Pair>();
                stack.addLast(new Pair(i, j));
                //周りの未探索領域をスタックに入れる
                Pair target;
                var size = 0;
                var mitansakuset = new LinkedList<Pair>();
                while ((target = stack.pollFirst()) != null) {
                    var v = result.get(target.a, target.b);
                    size = Math.max(size, v);
                    for (int[] d : dif) {
                        var x = target.a + d[0];
                        var y = target.b + d[1];
                        if (result.isRange(x, y) && result.isMined(x, y) && !visited[x][y]) {
                            if (result.isNoOilField(x, y)) {
                                continue;
                            }
                            visited[x][y] = true;
                            stack.addLast(new Pair(x, y));
                        }
                        if (result.isRange(x, y) && result.isNotMined(x, y)) {
                            mitansakuset.add(new Pair(x, y));
                        }
                    }
                }
                list.add(mitansakuset);
                islandCount += size;
            }
        }
        return list;
    }


    private static void dumpOffsetcount(MacroGuesser macroGuesser) {
        //todo:実験
        for (int i = 0; i < macroGuesser.offsetCount.length; i++) {
            System.out.println("油田:" + i);
            var pq = new PriorityQueue<Triple>(Comparator.<Triple>comparingLong(t -> t.c).reversed());
            for (int j = 0; j < macroGuesser.offsetCount[0].length; j++) {
                for (int k = 0; k < macroGuesser.offsetCount[0].length; k++) {
                    pq.add(new Triple(j, k, macroGuesser.offsetCount[i][j][k]));
                }
            }
            System.out.println(pq);
        }
    }

    public static class MacroGuesser {
        private final double[][] table;
        private final LinkedList<OilFieldSet> oilList;

        private final Random random = new Random(0);
        private Generator[] generators;

        private int dx = 0;
        private int dy = 0;
        private int[][] baseTable;

        //島のサイズの何割を良い解とするか
        private static double V_THRETHOLD_RATIO = 0.85;
        //許容違反数
        private static int ALLOWED_VIOLATION = 2;

        private static int GENERATE_NUMBER = 1200;

        //todo:実験
        private long[][][] offsetCount;

        private final double errorRatio;

        private List<Pair> maybeAnswer = null;
        private final HashSet<String> oldanswers = new HashSet<>();

        private int[][] probCount;
        private int probMax;

        private void genNumCal(int n, int m, double error) {
            //mが10以上なら、生成数を落とす
            GENERATE_NUMBER = (m < 8 ? 15000 : (m < 16 ? 13000 : 10000));
            //エラー率に応じて閾値を調整する
            ALLOWED_VIOLATION = 0;
            V_THRETHOLD_RATIO = (m < 5 ? (30 * error) : (m < 10 ? (30 * error) : (30 * error)));
        }

        public MacroGuesser(int n, LinkedList<OilFieldSet> oilList, double errorRatio) {
            this.table = new double[n][n];
            this.offsetCount = new long[oilList.size()][n][n];
            this.oilList = oilList;
            this.errorRatio = errorRatio;
            int size = oilList.stream().mapToInt(set -> set.size()).sum();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    table[i][j] = (double) size / (n * n);
                }
            }

            this.generators = new Generator[oilList.size()];
            for (int i = 0; i < oilList.size(); i++) {
                this.generators[i] = new Generator(n, oilList.get(i));
            }
            genNumCal(n, oilList.size(), errorRatio);

        }

        private void printout() {
            for (int i = 0; i < table.length; i++) {
                System.out.println(Arrays.toString(table[i]));
            }
        }

        //行列から確率を更新する.まだ未採掘状態のときに利用することを想定
        public void setup(InteractiveConnector connector) {
            var x = table.length;
            var y = table[0].length;

            OilFieldSet mini = null;
            //一番小さくて、縦横比が１に近いのを選択する
            for (OilFieldSet set : oilList) {
                if (mini == null) {
                    mini = set;
                    continue;
                }
                if (set.oil.length * set.oil[0].length < mini.oil.length * mini.oil[0].length) {
                    mini = set;
                } else if (set.oil.length * set.oil[0].length == mini.oil.length * mini.oil[0].length) {
                    var ds = (double) set.oil.length / set.oil[0].length;
                    var ms = (double) mini.oil.length / mini.oil[0].length;
                    if (Math.abs(1 - ds) < Math.abs(1 - ms)) {
                        mini = set;
                    }
                }
            }
            var dx = mini.oil.length;
            var dy = mini.oil[0].length;
            var baseTable = new int[x][y];

            this.dx = dx;
            this.dy = dy;
            this.baseTable = baseTable;

            for (int i = 0; i * dx < x; i++) {
                for (int j = 0; j * dy < y; j++) {
                    var list = new LinkedList<Pair>();
                    //はみ出てるならはみ出た分戻す
                    var hamix = Math.max(0, (i + 1) * dx - x);
                    var hamiy = Math.max(0, (j + 1) * dy - y);
                    for (int k = i * dx; k < (i + 1) * dx; k++) {
                        for (int l = j * dy; l < (j + 1) * dy; l++) {
                            list.add(new Pair(k - hamix, l - hamiy));
                        }
                    }
                    var v = connector.foresee(list);
                    baseTable[(i + 1) * dx - hamix - 1][(j + 1) * dy - 1 - hamiy] = v;
                }
            }


            //ランダムに候補を生成して初期の油田の位置を決める
            var genes = new PriorityQueue<OilV>(Comparator.<OilV>comparingLong(oilV -> oilV.value));
            for (int i = 0; i < 100000; i++) {
                var temp = generate();
                var v = calgosa(baseTable, dx, dy, temp);
                genes.add(new OilV(temp, v));
            }
            //上位100件を使って、それっぽい期待値を作成する
            var select = new LinkedList<OilV>();
            for (int i = 0; i < 30; i++) {
                select.add(genes.poll());
            }
            var simu = new int[x][y];
            for (OilV v : select) {
                var oil = v.sets.actual;
                for (int i = 0; i < oil.length; i++) {
                    for (int j = 0; j < oil[0].length; j++) {
                        simu[i][j] += oil[i][j];
                    }
                }
            }
            //良い解候補は積極的に確率を上げる
            for (var v : select) {
                //良い解かどうか
                if (v.value >= x * y * V_THRETHOLD_RATIO) {
                    continue;
                }
                //良い解は次でも選ばれやすくしたい
                var sets = v.sets.oilOffsets;
                var setidx = 0;
                for (OilOffset offset : sets) {
                    offsetCount[setidx][offset.offsetX][offset.offsetY]++;
                    setidx++;
                }
            }
            System.out.println(select.size());
            System.out.println("========初期解========");
            for (int i = 0; i < x; i++) {
                System.out.println(Arrays.toString(simu[i]));
            }

        }

        public void recal(MiningResult result) {
            var x = table.length;
            var y = table[0].length;

            var select = new TreeMap<OilV, Long>(Comparator.<OilV>comparingLong(oilV -> oilV.value).thenComparing(oilV -> oilV.sets.hash));
            var queue = new PriorityQueue<OilV>(Comparator.<OilV>comparingLong(oilV -> oilV.value));
            updateGenerators();
            for (int i = 0; i < 5000; i++) {
                var temp = generate();
                //resultと矛盾してないか確認
                var tgrid = temp.actual;
                var rgrid = result.grid;

                var v = calgosa(baseTable, dx, dy, temp);
                if (validSet(tgrid, rgrid) > ALLOWED_VIOLATION) {
                    //候補に入れない
                    queue.add(new OilV(temp, v));
                    continue;
                }
                select.merge(new OilV(temp, v), 1l, Long::sum);
            }



            var simu = new int[x][y];
            var sum = 0;

            //リセット
            if (select.size() != 0) {
                for (int i = 0; i < offsetCount.length; i++) {
                    for (int j = 0; j < offsetCount[i].length; j++) {
                        for (int k = 0; k < offsetCount[i][j].length; k++) {
                            offsetCount[i][j][k] = 0;
                        }
                    }
                }
            } else {
                for (int i = 0; i < offsetCount.length; i++) {
                    System.out.println("オフセット:" + i);
                    for (int j = 0; j < offsetCount[i].length; j++) {
                        System.out.println(Arrays.toString(offsetCount[i][j]));
                    }
                }
            }
            {
                //行けそうならやってみる誤差が許容範囲なら答えてみる
                if (validSet(select.firstKey().sets.actual, result.grid) == 0
                        && select.firstKey().value < x * y * V_THRETHOLD_RATIO) {
                    var grid = select.firstKey().sets.actual;
                    var list = new ArrayList<Pair>();
                    for (int i = 0; i < grid.length; i++) {
                        for (int j = 0; j < grid[0].length; j++) {
                            if (grid[i][j] != 0) {
                                list.add(new Pair(i, j));
                            }
                        }
                    }
                    this.maybeAnswer = list;
                    this.oldanswers.add(select.firstKey().sets.hash);
                }
            }

            for (var en : select.entrySet()) {
                var v = en.getKey();
                var k = en.getValue();
                var oil = v.sets.actual;
                for (int i = 0; i < oil.length; i++) {
                    for (int j = 0; j < oil[0].length; j++) {
                        simu[i][j] += (oil[i][j] > 0 ? 1 : 0) * k;
                        sum += (oil[i][j] > 0 ? 1 : 0) * k;
                    }
                }

                var sets = v.sets.oilOffsets;
                var setidx = 0;
                for (OilOffset offset : sets) {
                    offsetCount[setidx][offset.offsetX][offset.offsetY]++;
                    setidx++;
                }
            }

            System.out.println(select.size());
            System.out.println("================");
            for (int i = 0; i < x; i++) {
                System.out.println(Arrays.toString(simu[i]));
                for (int j = 0; j < y; j++) {
                    table[i][j] = (double) simu[i][j] / sum;
                }
            }
            probCount = simu;
        }

        //違反数チェック
        private int validSet(int[][] target, int[][] result) {
            var ret = 0;
            for (int i = 0; i < target.length; i++) {
                for (int j = 0; j < target[0].length; j++) {
                    if (result[i][j] == Integer.MAX_VALUE) {
                        //未探索はスキップ
                        continue;
                    }
                    if (result[i][j] != target[i][j]) {
                        //不一致
                        ret += Math.abs(result[i][j] - target[i][j]);
                    }
                }
            }
            return ret;
        }

        private List<Pair> pollMaybeAnswer() {
            var result = this.maybeAnswer;
            this.maybeAnswer = null;
            return result;
        }

        private long calgosa(int[][] basetable, int dx, int dy, SetOfOilOffset set) {
            var x = basetable.length;
            var y = basetable[0].length;
            var result = 0l;
            for (int i = 0; i * dx < x; i++) {
                for (int j = 0; j * dy < y; j++) {
                    var list = new LinkedList<Pair>();
                    //はみ出てるならはみ出た分戻す
                    var hamix = Math.max(0, (i + 1) * dx - x);
                    var hamiy = Math.max(0, (j + 1) * dy - y);
                    var b = basetable[(i + 1) * dx - hamix - 1][(j + 1) * dy - 1 - hamiy];
                    var t = set.count(i * dx - hamix, j * dy - hamiy, (i + 1) * dx - hamix - 1, (j + 1) * dy - 1 - hamiy);

                    var p = random.nextGaussian();
                    var k = dx * dy;
                    var mean = (k - t) * errorRatio + t * (1 - errorRatio);
                    var vari = k * (1 - errorRatio) * errorRatio;
                    var vs = Math.max(0, Math.round(p * vari + mean));
                    result += (b - vs) * (b - vs);
                }
            }
            return result;
        }

        //油田をランダムに島に配置していく
        private SetOfOilOffset generate() {
            var result = new OilOffset[this.generators.length];
            for (int i = 0; i < this.generators.length; i++) {
                result[i] = this.generators[i].generate(random);
            }
            return new SetOfOilOffset(result, this.table.length);
        }

        private void updateGenerators() {
            var result = new OilOffset[this.generators.length];
            var tatal = 1l;
            for (int i = 0; i < this.generators.length; i++) {
                this.generators[i].update(offsetCount[i]);
                tatal *= this.generators[i].size();
            }
            var count = 1l;
            for (int i = 0; i < this.generators.length; i++) {
                System.out.print("[" + i + "]" + this.generators[i].countProb() + ",");
                count *= this.generators[i].countProb();
            }
            System.out.println();
            System.out.println(count);
        }

        public void setNoOil(int i, int j) {
            for (Generator generator : generators) {
                generator.setImpossibeOffset(i, j);
            }
        }


        //確度の高い採掘地点を返す
        public Pair pollMiningPoint() {
            var i = 0;
            var j = 0;
            var prov = 0.0;
            for (int k = 0; k < table.length; k++) {
                for (int l = 0; l < table[0].length; l++) {
                    if (prov < table[k][l]) {
                        prov = table[k][l];
                        i = k;
                        j = l;
                    }
                }
            }
            return new Pair(i, j);
        }


        //確度の高いベストNを返却する。
        public List<Pair> highPoint(int n, MiningResult result) {
            //消しこむとパターンを絞りこめる候補を返す
            var max = 0;
            for (int i = 0; i < this.probCount.length; i++) {
                for (int j = 0; j < this.probCount.length; j++) {
                    if (max < this.probCount[i][j]) {
                        max = this.probCount[i][j];
                    }
                }
            }
            //半分くらいにパターンが割れているポイントを探す
            var half = max / 2;
            Pair ret = null;
            var now = Integer.MAX_VALUE;
            for (int i = 0; i < this.probCount.length; i++) {
                for (int j = 0; j < this.probCount.length; j++) {
                    if (Math.abs(half - now) >= Math.abs(half - this.probCount[i][j])) {
                        if (result.isNotMined(i, j)) {
                            if (Math.abs(half - now) == 0 && ret != null) {
                                if (random.nextDouble() <= 0.5) {
                                    continue;
                                }
                            }
                            now = this.probCount[i][j];
                            ret = new Pair(i, j);
                        }
                    }
                }
            }
            var list = new LinkedList<Pair>();
            list.add(ret);
            return list;
        }

        public void setOil(int i, int j) {
            for (Generator generator : generators) {
                generator.setWeighting(i, j);
            }
        }

        public boolean trySearch(InteractiveConnector connector, MiningResult result) {
            var n = table.length;
            var table = new int[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    table[i][j] = result.get(i, j);
                }
                System.out.println(Arrays.toString(table[i]));
            }
            var index = 0;
            var list = new LinkedList<OilOffset>();
            var genelist = new LinkedList<Generator>();
            for (Generator generator : this.generators) {
                genelist.add(generator);
            }
            Collections.sort(genelist, Comparator.comparingInt(generator -> generator.countProb()));

            for (Generator generator : genelist) {
                System.out.print(generator.countProb() + ",");
            }
            System.out.println();
            return dfs(genelist, table, list, connector);
        }

        private boolean dfs(LinkedList<Generator> generators, int[][] target, LinkedList<OilOffset> list, InteractiveConnector connector) {
            if (generators.size() == 0) {
                for (int i = 0; i < target.length; i++) {
                    for (int j = 0; j < target[0].length; j++) {
                        if (target[i][j] != Integer.MAX_VALUE && target[i][j] > 0) {
                            return false;
                        }
                    }
                }
                System.out.println("解を見つけた");
                System.out.println(list);
                var array = new OilOffset[list.size()];
                for (int i = 0; i < list.size(); i++) {
                    array[i] = list.get(i);
                }
                var test = new SetOfOilOffset(array, target.length);
                var table = test.actual;
                var ans = new LinkedList<Pair>();
                for (int i = 0; i < table.length; i++) {
                    for (int j = 0; j < table.length; j++) {
                        if (table[i][j] > 0) {
                            ans.add(new Pair(i, j));
                        }
                    }
                }
                var v = connector.answer(ans);
                return v;
            }


            var generator = generators.pollFirst();
            var offset = generator.oilFieldSet;
            var flg = generator.flg;
            var b = false;
            for (int i = 0; i < flg.length; i++) {
                for (int j = 0; j < flg[0].length; j++) {
                    if (flg[i][j]) {
                        var ok = true;
                        for (Pair p : offset.list) {
                            var x = p.a + i;
                            var y = p.b + j;
                            if (target[x][y] != Integer.MAX_VALUE) {
                                target[x][y]--;
                                if (target[x][y] < 0) {
                                    ok = false;
                                }
                            }
                        }
                        //okフラグ立っているならdfsする
                        if (ok) {
                            list.addLast(new OilOffset(generator.oilFieldSet, i, j));
                            b = dfs(generators, target, list, connector);
                            list.pollLast();
                        }
                        //戻ってきたら値を戻す
                        for (Pair p : offset.list) {
                            var x = p.a + i;
                            var y = p.b + j;
                            if (target[x][y] != Integer.MAX_VALUE) {
                                target[x][y]++;
                            }
                        }
                        if (b) {
                            generators.addFirst(generator);
                            return true;
                        }
                    }

                }
            }
            generators.addFirst(generator);
            return false;
        }
    }

    //ランダム油田
    record OilV(SetOfOilOffset sets, long value) {
    }

    public static class SetOfOilOffset {
        private final OilOffset[] oilOffsets;
        private final int[][] actual;
        private final int[][] acc;

        private String hash = "";

        public SetOfOilOffset(OilOffset[] oilOffsets, int n) {
            this.oilOffsets = oilOffsets;
            this.actual = new int[n][n];
            this.acc = new int[n + 1][n + 1];
            for (OilOffset offset : oilOffsets) {
                hash += "[" + offset.offsetX + ":" + offset.offsetY + "]";
            }
            for (OilOffset offset : oilOffsets) {
                //todo: ここの最適化
                var xoff = offset.offsetX;
                var yoff = offset.offsetY;
                var oilgird = offset.oilFieldSet.oil;
                for (int i = 0; i < oilgird.length; i++) {
                    for (int j = 0; j < oilgird[0].length; j++) {
                        actual[i + xoff][j + yoff] += oilgird[i][j] ? 1 : 0;
                    }
                }
            }
            for (int i = 1; i < acc.length; i++) {
                for (int j = 1; j < acc[0].length; j++) {
                    acc[i][j] = actual[i - 1][j - 1] + acc[i][j - 1];
                }
            }
            for (int i = 1; i < acc[0].length; i++) {
                for (int j = 1; j < acc.length; j++) {
                    acc[j][i] += acc[j - 1][i];
                }
            }
        }

        public int count(int sx, int sy, int ex, int ey) {
            return acc[ex + 1][ey + 1] + acc[sx][sy] - acc[ex + 1][sy] - acc[sx][ey + 1];
        }
    }

    record OilOffset(OilFieldSet oilFieldSet, int offsetX, int offsetY) {
        public boolean isOil(int i, int j) {
            return oilFieldSet.isOil(i - offsetX, j - offsetY);
        }
    }

    public static class Generator {
        private final int n;
        private final OilFieldSet oilFieldSet;

        private final double[][] table;
        //ありえるオフセットがどこかを保持しておく
        private final boolean[][] flg;
        //そのためには島の油田が0であった場所の情報が必要
        private final int[][] weighting;

        public int[][] getWeighting() {
            return weighting;
        }

        public Generator(int n, OilFieldSet oilFieldSet) {
            this.n = n;
            this.oilFieldSet = oilFieldSet;
            var xo = n - oilFieldSet.oil.length + 1;
            var yo = n - oilFieldSet.oil[0].length + 1;
            table = new double[xo][yo];
            flg = new boolean[xo][yo];
            weighting = new int[xo][yo];
            for (int i = 0; i < xo; i++) {
                for (int j = 0; j < yo; j++) {
                    table[i][j] = 1.0 / (xo * yo);
                    flg[i][j] = true;
                    weighting[i][j] = 1;
                }
            }
        }

        public int countProb() {
            var cnt = 0;
            for (int i = 0; i < flg.length; i++) {
                for (int j = 0; j < flg[0].length; j++) {
                    cnt += flg[i][j] ? 1 : 0;
                }
            }
            return cnt;
        }

        public OilOffset generate(Random random) {
            var d = random.nextDouble();
            var sum = 0.0d;
            var offsetX = 0;
            var offsetY = 0;
            //todo:二分探索で高速化する
            LOOP:
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < table[0].length; j++) {
                    sum += table[i][j];
                    if (d <= sum) {
                        offsetX = i;
                        offsetY = j;
                        break LOOP;
                    }
                }
            }
            return new OilOffset(oilFieldSet, offsetX, offsetY);
        }

        /**
         * 油田の存在しない島上の地点i,jが与えられたときに
         * オフセットとして不可なポイントを設定する
         */
        public void setImpossibeOffset(int i, int j) {
            for (Pair p : oilFieldSet.list) {
                var x = i - p.a;
                var y = j - p.b;
                if (0 <= x && x < flg.length && 0 <= y && y < flg[0].length) {
                    flg[x][y] = false;
                }
            }
        }

        /**
         * 油田が存在している島上の地点i,jが与えられたときに
         * そこを含むオフセットが高確率で選ばれるように重みつけを行う
         */
        public void setWeighting(int i, int j) {
            for (Pair p : oilFieldSet.list) {
                var x = i - p.a;
                var y = j - p.b;
                if (0 <= x && x < flg.length && 0 <= y && y < flg[0].length) {
                    weighting[x][y] += 1;
                }
            }
        }

        private void update(long[][] count) {
            var cnt = 0l;

            for (int i = 0; i < flg.length; i++) {
                for (int j = 0; j < flg[0].length; j++) {
                    var w1 = 1;//weighting[i][j];
                    var w2 = (long) count[i][j] != 0 ? count[i][j] : 1;
                    cnt += flg[i][j] ? ((double) w1 * w2) * 1.0 : 0;
                }
            }
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < table[0].length; j++) {
                    var w1 = 1;//weighting[i][j];
                    var w2 = (long) count[i][j] != 0 ? count[i][j] : 1;
                    if (flg[i][j]) {
                        table[i][j] = ((double) w1 * w2) * 1.0 / cnt;
                    } else {
                        table[i][j] = 0.0;
                    }
                }
            }
        }

        public long size() {
            var size = 0;
            for (int i = 0; i < flg.length; i++) {
                for (int j = 0; j < flg[0].length; j++) {
                    size += flg[i][j] && weighting[i][j] >= this.oilFieldSet.size ? 1 : 0;
                }
            }
            return size;
        }
    }

    public static class MicroGuesser {

    }

    public static class OilFieldSet {
        private final boolean[][] oil;
        private int size;
        private final List<Pair> list;

        public OilFieldSet(List<Pair> ps) {
            list = ps;
            var maxX = ps.stream().mapToInt(p -> p.a).max().getAsInt();
            var maxY = ps.stream().mapToInt(p -> p.b).max().getAsInt();
            this.oil = new boolean[maxX + 1][maxY + 1];
            for (Pair p : ps) {
                oil[p.a][p.b] = true;
            }
            size = ps.size();
        }

        public int size() {
            return size;
        }

        public boolean isOil(int i, int j) {
            if (0 <= i && i < oil.length && 0 <= j && j < oil[0].length) {
                return oil[i][j];
            }
            return false;
        }
    }

    public static class MiningResult {
        private final int[][] grid;

        public MiningResult(int n) {
            this.grid = new int[n][n];
            for (int i = 0; i < n; i++) {
                Arrays.fill(grid[i], Integer.MAX_VALUE);
            }
        }

        public void set(int i, int j, int v) {
            valid(i, j);
            if (0 <= i && i < grid.length && 0 <= j && j < grid[0].length) {
                grid[i][j] = v;
            }
        }

        public int get(int i, int j) {
            valid(i, j);
            return grid[i][j];
        }

        public boolean isRange(int i, int j) {
            return 0 <= i && i < grid.length && 0 <= j && j < grid[0].length;
        }

        public boolean isMined(int i, int j) {
            valid(i, j);
            return grid[i][j] != Integer.MAX_VALUE;
        }

        public boolean isNotMined(int i, int j) {
            valid(i, j);
            return !isMined(i, j);
        }

        public boolean isOilField(int i, int j) {
            valid(i, j);
            return isMined(i, j) && grid[i][j] > 0;
        }

        public boolean isNoOilField(int i, int j) {
            valid(i, j);
            return isMined(i, j) && grid[i][j] == 0;
        }

        private void valid(int i, int j) {
            assert 0 <= i && i < grid.length && 0 <= j && j < grid[0].length;
        }

        public List<Pair> resultList() {
            var list = new LinkedList<Pair>();
            for (int i = 0; i < grid.length; i++) {
                for (int j = 0; j < grid.length; j++) {
                    if (isOilField(i, j)) {
                        list.add(new Pair(i, j));
                    }
                }
            }
            return list;
        }
    }

    public static class InteractiveConnector {
        private final PrintWriter pw;
        private final FastScanner fs;

        private int counter = 0;

        public InteractiveConnector(PrintWriter pw, FastScanner fs) {
            this.pw = pw;
            this.fs = fs;
        }

        //占う
        public int foresee(List<Pair> list) {
            pw.print("q " + list.size() + " ");
            for (Pair p : list) {
                pw.print(p.a + " " + p.b + " ");
            }
            pw.println();
            pw.flush();
            counter++;
            return fs.ni();
        }

        //採掘する
        public int mine(Pair p) {
            pw.println("q 1 " + p.a + " " + p.b);
            pw.flush();
            counter++;
            return fs.ni();
        }

        //回答する

        public boolean answer(List<Pair> list) {
            pw.print("a " + list.size() + " ");
            for (Pair p : list) {
                pw.print(p.a + " " + p.b + " ");
            }
            pw.println();
            pw.flush();
            counter++;
            return fs.ni() == 1;
        }
    }

    //--------------
    record TPair<S, T>(S a, T b) {
    }

    record TTri<S, T, U>(S a, T b, U c) {
    }

    record Pair(int a, int b) {
    }

    record Triple(int a, int b, long c) {
    }

    public static int[] concat(int[] a, int[] b) {
        var ret = new int[a.length + b.length];
        for (int i = 0; i < a.length; i++) {
            ret[i] = a[i];
        }
        for (int i = 0; i < b.length; i++) {
            ret[i + a.length] = b[i];
        }
        return ret;
    }

    public static long[] concat(long[] a, long[] b) {
        var ret = new long[a.length + b.length];
        for (int i = 0; i < a.length; i++) {
            ret[i] = a[i];
        }
        for (int i = 0; i < b.length; i++) {
            ret[i + a.length] = b[i];
        }
        return ret;
    }


    public static void solve(InputStream in, PrintStream out) {
        PrintWriter pw = new PrintWriter(out);
        FastScanner fs = new FastScanner(in);
        try {
            solve(pw, fs);
        } finally {
            pw.flush();
        }
    }


//-------------------------------------------------------------------

    /**
     * 各インデックスが配列の長さ以内に収まっているか境界チェックを行う
     * <p>
     * 多次元配列のチェックをいちいち書くのがしんどい時に。
     * arrayBound(new int[]{1,2,3} , new int[]{3,4,3})
     *
     * @param target 配列に設定したいインデックス
     * @param len    　配列の長さ
     * @return 配列の長さ内に収まってる時 true
     */
    private static boolean inarr(int[] target, int[] len) {
        var b = true;
        if (target.length != len.length) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < target.length; i++) {
            b &= (0 <= target[i] && target[i] < len[i]);
        }
        return b;
    }

    private static int[] arr(int... a) {
        return Arrays.copyOf(a, a.length);
    }

    private static long[] arr(long... a) {
        return Arrays.copyOf(a, a.length);
    }

    //半時計90度回転
    private static int[][] rot(int[][] grid) {
        var h = grid.length;
        var w = grid[0].length;

        var result = new int[w][h];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                result[w - 1 - j][i] = grid[i][j];
            }
        }
        return result;
    }


//http://fantom1x.blog130.fc2.com/blog-entry-194.html

    /**
     * <h1>指定した値以上の先頭のインデクスを返す</h1>
     * <p>配列要素が０のときは、０が返る。</p>
     *
     * @param arr   ： 探索対象配列(単調増加であること)
     * @param value ： 探索する値
     * @return<b>int</b> ： 探索した値以上で、先頭になるインデクス
     */
    public static final int lowerBound(final long[] arr, final long value) {
        int low = 0;
        int high = arr.length;
        int mid;
        while (low < high) {
            mid = ((high - low) >>> 1) + low;    //(low + high) / 2 (オーバーフロー対策)
            if (arr[mid] < value) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        return low;
    }

    /**
     * <h1>指定した値より大きい先頭のインデクスを返す</h1>
     * <p>配列要素が０のときは、０が返る。</p>
     *
     * @param arr   ： 探索対象配列(単調増加であること)
     * @param value ： 探索する値
     * @return<b>int</b> ： 探索した値より上で、先頭になるインデクス
     */
    public static final int upperBound(final long[] arr, final long value) {
        int low = 0;
        int high = arr.length;
        int mid;
        while (low < high) {
            mid = ((high - low) >>> 1) + low;    //(low + high) / 2 (オーバーフロー対策)
            if (arr[mid] <= value) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        return low;
    }

    //----------------
    public static class FastScanner {
        InputStream in;
        byte[] buffer = new byte[1 << 10];
        int length = 0;
        int ptr = 0;
        private final Predicate<Byte> isPrintable;


        public FastScanner(InputStream in) {
            this.in = in;
            this.isPrintable = b -> (33 <= b && b <= 126);
        }

        public FastScanner(InputStream in, Predicate<Byte> predicate) {
            this.in = in;
            this.isPrintable = predicate;
        }

        private boolean hasNextByte() {
            if (ptr < length) {
                return true;
            }
            try {
                length = in.read(buffer);
            } catch (IOException e) {
                e.printStackTrace();
            }
            ptr = 0;
            return length != 0;
        }


        private byte read() {
            if (hasNextByte()) {
                return buffer[ptr++];
            }
            return 0;
        }

        private void skip() {
            while (hasNextByte() && !isPrintable(buffer[ptr])) {
                ptr++;
            }
        }

        private boolean hasNext() {
            skip();
            return hasNextByte();
        }

        private boolean isPrintable(byte b) {
            return 33 <= b && b <= 126;
        }


        private String innerNext(Predicate<Byte> isReadable) {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            StringBuilder sb = new StringBuilder();
            byte b = read();
            while (isReadable.test(b)) {
                sb.appendCodePoint(b);
                b = read();
            }
            return sb.toString();
        }

        public String n() {
            return innerNext(b -> (33 <= b && b <= 126));
        }

        public int ni() {
            return (int) nl();
        }

        public char[][] ncaa(int n, int m) {
            var grid = new char[n][m];
            for (int i = 0; i < n; i++) {
                grid[i] = n().toCharArray();
            }
            return grid;
        }

        public int[] nia(int n) {
            int[] result = new int[n];
            for (int i = 0; i < n; i++) {
                result[i] = ni();
            }
            return result;
        }

        public int[][] niaa(int h, int w) {
            int[][] result = new int[h][w];
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    result[i][j] = ni();
                }
            }
            return result;
        }

        public long[][] nlaa(int h, int w) {
            long[][] result = new long[h][w];
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    result[i][j] = nl();
                }
            }
            return result;
        }

        public String[] na(int n) {
            String[] result = new String[n];
            for (int i = 0; i < n; i++) {
                result[i] = n();
            }
            return result;
        }

        public long[] nla(int n) {
            long[] result = new long[n];
            for (int i = 0; i < n; i++) {
                result[i] = nl();
            }
            return result;
        }

        public long nl() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            long result = 0;
            boolean minus = false;
            byte b;

            b = read();
            if (b == '-') {
                minus = true;
                b = read();
            }

            while (isPrintable(b)) {
                if (b < '0' || b > '9') {
                    throw new NumberFormatException();
                }
                result *= 10;
                result += (b - '0');
                b = read();
            }

            return minus ? -result : result;
        }
    }

//-------------------------------------------------------------------
}
