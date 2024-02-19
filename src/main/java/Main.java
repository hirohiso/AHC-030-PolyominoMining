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
        while (connector.counter < 2 * N * N - 1) {
            break;
        }
        //回答フェーズ
        var list = result.resultList();
        //connector.answer(list);
        pw.println("#break");
        /**/
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
        private int[][] countTable;

        private final double errorRatio;

        private List<Pair> maybeAnswer = null;
        private final HashSet<String> oldanswers = new HashSet<>();

        private void genNumCal(int n, int m, double error) {
            //mが10以上なら、生成数を落とす
            GENERATE_NUMBER = 1000;
            //エラー率に応じて閾値を調整する
            //V_THRETHOLD_RATIO = 1 - error;
            ALLOWED_VIOLATION = 0;
            V_THRETHOLD_RATIO = (m < 5 ? 0.8 : (m < 10 ? 0.9 : 1));
            //System.out.println(GENERATE_NUMBER);
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
            var totalsize = 0;
            //一番小さくて、縦横比が１に近いのを選択する
            for (OilFieldSet set : oilList) {
                totalsize += set.size;
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
            var dx = Math.max(mini.oil.length, 3);
            var dy = Math.max(mini.oil[0].length, 3);
            var baseTable = new int[x][y];

            this.dx = dx;
            this.dy = dy;
            this.baseTable = baseTable;

            var rowcount = new int[x];
            var colcount = new int[x];
            for (int i = 0; i < x; i++) {
                var list = new LinkedList<Pair>();
                for (int j = 0; j < y; j++) {
                    list.add(new Pair(i, j));
                }
                rowcount[i] = connector.foresee(list);
            }
            for (int i = 0; i < x; i++) {
                var list = new LinkedList<Pair>();
                for (int j = 0; j < y; j++) {
                    list.add(new Pair(j, i));
                }
                colcount[i] = connector.foresee(list);
            }


            //行単独で推定
            var test = new int[x * y][x * y * 2];
            for (int i = 0; i < x * y; i++) {
                for (int j = 0; j < 1000; j++) {
                    var v = (int) test(i, x);
                    test[i][v]++;
                }
                //System.out.println(i + ":" + Arrays.toString(test[i]));
            }
            //minの推移
            var mintable = new int[x];
            var maxtable = new int[x];
            Arrays.fill(maxtable, Integer.MAX_VALUE);
            {
                var min = 0;
                for (int i = 0; i < rowcount.length; i++) {
                    for (int j = 0; j < test[i].length; j++) {
                        if (test[rowcount[i]][j] != 0) {
                            min += j;
                            mintable[i] = Math.max(mintable[i], min);
                            break;
                        }
                    }
                }
            }
            {
                var max = 0;
                for (int i = 0; i < rowcount.length; i++) {
                    for (int j = test[i].length - 1; j >= 0; j--) {
                        if (test[rowcount[i]][j] != 0) {
                            max += j;
                            maxtable[i] = Math.min(maxtable[i], max);
                            break;
                        }
                    }
                }
            }
            //逆方向からの繊維
            {
                var max = totalsize;
                maxtable[x - 1] = Math.min(maxtable[x - 1], max);
                for (int i = rowcount.length - 1; i > 0; i--) {
                    for (int j = 0; j < test[i].length; j++) {
                        if (test[rowcount[i]][j] != 0) {
                            max -= j;
                            maxtable[i - 1] = Math.min(maxtable[i - 1], max);
                            break;
                        }
                    }
                }
            }
            {
                var max = totalsize;
                mintable[x - 1] = Math.max(mintable[x - 1], max);
                for (int i = rowcount.length - 1; i > 0; i--) {
                    for (int j = test[i].length - 1; j >= 0; j--) {
                        if (test[rowcount[i]][j] != 0) {
                            max -= j;
                            mintable[i - 1] = Math.max(mintable[i - 1], max);
                            break;
                        }
                    }
                }
            }

            var rowprob = new double[x][y * x];
            //1行目
            for (int i = 0; i < x; i++) {
                rowprob[0][i] = (double) test[rowcount[0]][i] / 1000;
            }

            for (int i = 0; i < rowprob.length - 1; i++) {
                for (int j = 0; j < rowprob[i].length; j++) {
                    for (int k = 0; k < test.length; k++) {
                        var prob = (double) test[rowcount[i + 1]][k] / 1000;
                        rowprob[i + 1][Math.min(j + k, rowprob[i].length - 1)] += (prob * rowprob[i][j]);
                    }
                }

                var sum = 0.0;
                for (int j = mintable[i + 1]; j <= maxtable[i + 1]; j++) {
                    sum += rowprob[i + 1][j];
                }
                for (int j = 0; j < rowprob[i + 1].length; j++) {
                    if (j >= mintable[i + 1] && j <= maxtable[i + 1]) {
                        rowprob[i + 1][j] /= sum;
                    } else {
                        rowprob[i + 1][j] = 0;
                    }

                }
            }


            //10000件ランダムに候補を生成してそれっぽいものを使って初期の確度を決める
            var genes = new PriorityQueue<OilV>(Comparator.<OilV>comparingLong(oilV -> oilV.value));
            //todo: 初期に何件生成するか

            var marCnt = new int[x * y + 1][x * y + 1];
            var marProv = new double[x * y + 1][x * y + 1];
            for (int i = 0; i < 100000; i++) {
                var temp = generate();
                var rowSum = temp.rowSum;
                //初期
                marCnt[0][rowSum[0]]++;
                for (int j = 1; j < rowcount.length; j++) {
                    marCnt[rowSum[j - 1]][rowSum[j]]++;
                }
            }
            for (int i = 0; i < marCnt.length; i++) {
                var cnt = 0;
                for (int j = 0; j < marCnt[i].length; j++) {
                    cnt += marCnt[i][j];
                }
                if (cnt == 0) {
                    continue;
                }
                for (int j = 0; j < marCnt[i].length; j++) {
                    marProv[i][j] = (double) marCnt[i][j] / cnt;
                }
            }

            //初期化
            var probfromMar = new double[x][x * y + 1];
            for (int i = 0; i < probfromMar[0].length; i++) {
                if (test[rowcount[0]][i] != 0) {
                    probfromMar[0][i] = marProv[0][i];
                }
            }

            for (int i = 1; i < rowcount.length; i++) {
                for (int j = 0; j < test[rowcount[i - 1]].length; j++) {
                    for (int k = 0; k < test[rowcount[i]].length; k++) {
                        if (test[rowcount[i - 1]][j] * test[rowcount[i]][k] != 0) {
                            probfromMar[i][k] += (marProv[j][k] * probfromMar[i - 1][j]);
                        }
                    }
                }
                var temp = 0.0;
                for (int j = 0; j < probfromMar[i].length; j++) {
                    temp += probfromMar[i][j];
                }
                for (int j = 0; j < probfromMar[i].length; j++) {
                    probfromMar[i][j] /= temp;
                }
            }

            var expectProb = new double[x][x * y];

            for (int i = 0; i < expectProb.length; i++) {
                if (i == 0) {
                    for (int j = 0; j < expectProb[i].length; j++) {
                        expectProb[i][j] = probfromMar[i][j] * rowprob[i][j];
                    }
                    var temp = 0.0;
                    for (int j = 0; j < expectProb[i].length; j++) {
                        temp += expectProb[i][j];
                    }
                    for (int j = 0; j < expectProb[i].length; j++) {
                        expectProb[i][j] /= temp;
                    }
                } else {
                    for (int j = 0; j < expectProb[i - 1].length; j++) {
                        for (int k = 0; k < probfromMar[i].length; k++) {
                            if (j + k < expectProb[i].length) {
                                expectProb[i][j + k] += (probfromMar[i][k] * expectProb[i - 1][j]);
                            }
                        }
                    }
                    for (int j = 0; j < expectProb[i].length; j++) {
                        expectProb[i][j] = rowprob[i][j] * expectProb[i][j];
                    }
                    var temp = 0.0;
                    for (int j = 0; j < expectProb[i].length; j++) {
                        temp += expectProb[i][j];
                    }
                    for (int j = 0; j < expectProb[i].length; j++) {
                        expectProb[i][j] /= temp;
                    }
                }
            }

            System.out.println("=========");
            for (int i = 0; i < expectProb.length; i++) {
                for (int j = 0; j < expectProb[i].length; j++) {
                    System.out.print(String.format("%.2f,", expectProb[i][j] * 100));
                }
                System.out.println();
            }
        }


        private long test(int sv, int k) {
            var p = random.nextGaussian();
            var mean = (k - sv) * errorRatio + sv * (1.0 - errorRatio);
            var vari = k * (1.0 - errorRatio) * errorRatio;
            var vs = Math.max(0, Math.round(p * vari + mean));
            return vs;
        }

        //油田をランダムに島に配置していく
        private SetOfOilOffset generate() {
            var result = new OilOffset[this.generators.length];
            for (int i = 0; i < this.generators.length; i++) {
                result[i] = this.generators[i].generate(random);
            }
            return new SetOfOilOffset(result, this.table.length);
        }


        public void setNoOil(int i, int j) {
            for (Generator generator : generators) {
                generator.setImpossibeOffset(i, j);
            }
        }


        public double ratio(int i, int j) {
            return table[i][j];
        }

        public void setOil(int i, int j) {
            for (Generator generator : generators) {
                generator.setWeighting(i, j);
            }
        }
    }

    //ランダム油田
    record OilV(SetOfOilOffset sets, long value) {
    }

    public static class SetOfOilOffset {
        private final OilOffset[] oilOffsets;
        private final int[][] actual;

        private final int[] rowSum;
        private final int[] colSum;

        private String hash = "";

        public SetOfOilOffset(OilOffset[] oilOffsets, int n) {
            this.oilOffsets = oilOffsets;
            this.actual = new int[n][n];
            this.rowSum = new int[n];
            this.colSum = new int[n];
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

            //rowsum/colsum
            for (int i = 0; i < actual.length; i++) {
                for (int j = 0; j < actual[i].length; j++) {
                    rowSum[i] += actual[i][j];
                    colSum[j] += actual[i][j];
                }
            }
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

        public void printTable() {
            System.out.println("=================");
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < table[0].length; j++) {
                    System.out.print(String.format("%.2f,", table[i][j] * 100));
                }
                System.out.println();
            }
            System.out.println("=================");
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
