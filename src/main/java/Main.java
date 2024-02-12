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
        var e = fs.n();
        var rand = new Random(0);

        var connector = new InteractiveConnector(pw, fs);

        var oilList = new LinkedList<OilFieldSet>();
        var max = 0;
        for (int i = 0; i < M; i++) {
            var d = fs.ni();
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
        var macroGuesser = new MacroGuesser(N, oilList);
        macroGuesser.setup(connector);
        var cnt = 0;
        var stack = new LinkedList<Pair>();
        var counter = 2 * N * N;
        LOOP:
        while (counter > 0) {
            var list = macroGuesser.highPoint(max, result);
            for (Pair pair : list) {
                var v = connector.mine(pair);
                counter--;
                result.set(pair.a, pair.b, v);
                cnt += v;
                if( v == 0){
                    macroGuesser.setNoOil(pair.a, pair.b);
                    //最初は積極的に更新しつつ、最後の方はあまり更新しないようにする
                    macroGuesser.recal();
                    break;
                }else {
                    macroGuesser.setOil(pair.a, pair.b);
                }
                //油田を全部見つけたら終了
                if( cnt >= max){
                    break LOOP;
                }

            }
        }
        //回答フェーズ
        var list = result.resultList();
        connector.answer(list);
    }

    public static class MacroGuesser {
        private final double[][] table;
        private final LinkedList<OilFieldSet> oilList;

        private final Random random = new Random(0);
        private Generator[] generators;

        private int dx = 0;
        private int dy = 0;
        private int[][] baseTable;

        public MacroGuesser(int n, LinkedList<OilFieldSet> oilList) {
            this.table = new double[n][n];
            this.oilList = oilList;
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


            //1000件ランダムに候補を生成してそれっぽいものを使って初期の確度を決める
            var genes = new PriorityQueue<OilV>(Comparator.<OilV>comparingLong(oilV -> oilV.value));
            //todo: 初期に何件生成するか
            for (int i = 0; i < 3000; i++) {
                var temp = generate();
                var v = calgosa(baseTable, dx, dy, temp);
                genes.add(new OilV(temp, v));
            }
            //上位３０件を使って、それっぽい期待値を作成する
            var select = new LinkedList<OilV>();
            //todo: 上位何件利用するか
            for (int i = 0; i < 30; i++) {
                select.add(genes.poll());
            }
            var simu = new int[x][y];
            var sum = 0;
            for (OilV v : select) {
                var oil = v.sets.actual;
                for (int i = 0; i < oil.length; i++) {
                    for (int j = 0; j < oil[0].length; j++) {
                        simu[i][j] += oil[i][j];
                        sum += oil[i][j];
                    }
                }
            }
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    table[i][j] = (double) simu[i][j] / sum;
                }
            }
        }

        public void recal(){
            var x = table.length;
            var y = table[0].length;

            //todo: 何件生成するか
            var genes = new PriorityQueue<OilV>(Comparator.<OilV>comparingLong(oilV -> oilV.value));
            for (int i = 0; i < 3000; i++) {
                var temp = generate();
                var v = calgosa(baseTable, dx, dy, temp);
                genes.add(new OilV(temp, v));
            }
            //上位10件を使って、それっぽい期待値を作成する
            //todo: 上位何件利用するか
            var select = new LinkedList<OilV>();
            for (int i = 0; i < 30; i++) {
                select.add(genes.poll());
            }
            var max = select.getLast().value;
            var simu = new int[x][y];
            var sum = 0;
            for (OilV v : select) {
                var oil = v.sets.actual;
                for (int i = 0; i < oil.length; i++) {
                    for (int j = 0; j < oil[0].length; j++) {
                        var w = max / (v.value + 1);
                        simu[i][j] += oil[i][j] * w;
                        sum += oil[i][j] * w;
                    }
                }
            }
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    table[i][j] = (double) simu[i][j] / sum;
                }
            }
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
                    var t = set.count(i * dx, j * dy, (i + 1) * dx - hamix - 1, (j + 1) * dy - 1 - hamiy);
                    result += (b - t) * (b - t);
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

        public void setNoOil(int i, int j){
            for (Generator generator : generators){
                generator.setImpossibeOffset(i,j);
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
            var x = table.length;
            var y = table[0].length;
            var array = new Pair[x * y];
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    array[i * y + j] = new Pair(i, j);
                }
            }
            Arrays.sort(array, Comparator.<Pair>comparingDouble(p -> table[p.a][p.b]).reversed());
            var list = new LinkedList<Pair>();
            var cnt = 0;
            for (int i = 0; i < array.length; i++) {
                if (result.isNotMined(array[i].a, array[i].b)) {
                    list.add(array[i]);
                    cnt++;
                    if (cnt == n) {
                        break;
                    }
                }
            }
            return list;
        }

        public void setOil(int i, int j) {
            for (Generator generator : generators){
                generator.setWeighting(i,j);
            }
        }
    }

    //ランダム油田
    record OilV(SetOfOilOffset sets, long value) {
    }

    public static class SetOfOilOffset {
        private final OilOffset[] oilOffsets;
        private final int[][] actual;
        private final int[][] acc;

        public SetOfOilOffset(OilOffset[] oilOffsets, int n) {
            this.oilOffsets = oilOffsets;
            this.actual = new int[n][n];
            this.acc = new int[n + 1][n + 1];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    for (OilOffset offset : oilOffsets) {
                        actual[i][j] += offset.isOil(i, j) ? 1 : 0;
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

        public OilOffset generate(Random random) {
            update();
            var d = random.nextDouble();
            var sum = 0.0d;
            var offsetX = 0;
            var offsetY = 0;
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
                    weighting[x][y] += 10;
                }
            }
        }

        private void update() {
            var cnt = 0;
            for (int i = 0; i < flg.length; i++) {
                for (int j = 0; j < flg[0].length; j++) {
                    cnt += flg[i][j] ? weighting[i][j] : 0;
                }
            }
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < table[0].length; j++) {
                    if(flg[i][j]){
                        table[i][j] = weighting[i][j] * 1.0 / (double)cnt;
                    }else {
                        table[i][j] = 0.0;
                    }
                }
            }
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
            return fs.ni();
        }

        //採掘する
        public int mine(Pair p) {
            pw.println("q 1 " + p.a + " " + p.b);
            pw.flush();
            return fs.ni();
        }

        //回答する

        public void answer(List<Pair> list) {
            pw.print("a " + list.size() + " ");
            for (Pair p : list) {
                pw.print(p.a + " " + p.b + " ");
            }
            pw.println();
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
