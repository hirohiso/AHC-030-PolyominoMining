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
        while (true) {
            var i = -1;
            var j = -1;
            if (stack.size() == 0) {
                //マクロ探索
                var p = macroGuesser.highPoint(1, result).get(0);
                i = p.a;
                j = p.b;
            } else {
                //マイクロ探索
                var p = stack.pollFirst();
                i = p.a;
                j = p.b;
            }
            if (result.isMined(i, j)) {
                continue;
            }

            var v = connector.mine(new Pair(i, j));
            result.set(i, j, v);
            macroGuesser.update(result);

            if (v != 0) {
                //マイクロ探索：見つけた油田の上下左右を確率的に探索し続ける
                var d = new int[][]{{-1, 0}, {0, -1}, {1, 0}, {0, 1}};
                for (int k = 0; k < d.length; k++) {
                    var x = i + d[k][0];
                    var y = j + d[k][1];
                    if (0 <= x && x < N && 0 <= y && y < N) {
                        if (result.isMined(x, y)) {
                            continue;
                        }
                        if(macroGuesser.table[x][y] > 0){
                            stack.addLast(new Pair(x, y));
                        }
                    }
                }
            }

            //油田を全部見つけたら終了
            cnt += v;
            if (cnt == max) {
                break;
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

        private int[] setupLines;
        private int[] setupColumns;

        private int[][] roundzero;

        public MacroGuesser(int n, LinkedList<OilFieldSet> oilList) {
            this.table = new double[n][n];
            this.oilList = oilList;
            this.roundzero = new int[n][n];
            int size = oilList.stream().mapToInt(set -> set.size()).sum();
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    table[i][j] = (double) size / (n * n);
                }
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

            //行列ごとに占って油田あるかもマスを初期化
            var lines = new int[x];
            {
                for (int i = 0; i < x; i++) {
                    var list = new LinkedList<Pair>();
                    for (int j = 0; j < y; j++) {
                        list.add(new Pair(i, j));
                    }
                    lines[i] = connector.foresee(list);
                }
            }

            var columns = new int[y];
            {
                for (int i = 0; i < y; i++) {
                    var list = new LinkedList<Pair>();
                    for (int j = 0; j < x; j++) {
                        list.add(new Pair(j, i));
                    }
                    columns[i] = connector.foresee(list);
                }
            }

            //1000件ランダムに候補を生成してそれっぽいものを使って初期の確度を決める
            var genes = new LinkedList<OilV>();
            for (int i = 0; i < 1000; i++) {
                var temp = generate();
                //縦横の二乗誤差を計算
                var xgosa = 0;
                for (int j = 0; j < temp.length; j++) {
                    var cnt = 0;
                    for (int k = 0; k < temp[0].length; k++) {
                        cnt += temp[j][k];
                    }
                    var d = lines[j] - cnt;
                    xgosa += d * d;
                }
                var ygosa = 0;
                for (int j = 0; j < temp.length; j++) {
                    var cnt = 0;
                    for (int k = 0; k < temp[0].length; k++) {
                        cnt += temp[k][j];
                    }
                    var d = columns[j] - cnt;
                    ygosa += d * d;
                }
                genes.add(new OilV(temp, xgosa, ygosa));
            }
            Collections.sort(genes, Comparator.comparingLong(o -> o.xvalue * o.yvalue));
            //上位３０件を使って、それっぽい期待値を作成する
            var select = new LinkedList<OilV>();
            for (int i = 0; i < 30; i++) {
                select.add(genes.get(i));
            }
            var temp = new long[x][y];
            var count = 0l;
            for (OilV o : select) {
                var t = o.olis;
                for (int i = 0; i < t.length; i++) {
                    for (int j = 0; j < t[0].length; j++) {
                        //誤差が少ないほど高い重みつけをする
                        var v = t[i][j];
                        temp[i][j] += v;
                        count += v;
                    }
                }
            }
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < table[0].length; j++) {
                    table[i][j] = (double) temp[i][j] / count;
                }
            }

            //あとで使うので占い結果は保持しておく
            setupColumns = columns;
            setupLines = lines;
        }

        //油田をランダムに島に配置していく
        private int[][] generate() {
            var xi = table.length;
            var yi = table[0].length;
            var result = new int[table.length][table[0].length];
            for (OilFieldSet oilf : oilList) {
                var xo = oilf.oil.length;
                var yo = oilf.oil[0].length;
                var offset = random.nextInt((xi - xo + 1) * (yi - yo + 1));
                var xoff = offset / (yi - yo + 1);
                var yoff = offset % (yi - yo + 1);
                for (int i = 0; i < xo; i++) {
                    for (int j = 0; j < yo; j++) {
                        result[i + xoff][j + yoff] += oilf.oil[i][j] ? 1 : 0;
                    }
                }
            }
            return result;
        }

        //採掘結果を元に確率を更新する
        public void update(MiningResult result) {
            var x = table.length;
            var y = table[0].length;
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    var d = new int[][]{{-1, 0}, {0, -1}, {1, 0}, {0, 1}};
                    var cntz = 0;
                    var mi = 0;
                    for (int k = 0; k < d.length; k++) {
                        var ti = i + d[k][0];
                        var tj = j + d[k][1];
                        if (result.isRange(ti, tj)) {
                            if (result.isNoOilField(ti, tj)) {
                                cntz++;
                            } else if (result.isNotMined(ti, tj)) {
                                mi++;
                            }
                        }
                    }
                    if (roundzero[i][j] != cntz) {
                        table[i][j] *= Math.pow(table[i][j], cntz);
                        roundzero[i][j] = cntz;
                    }
                }
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
    }

    record OilV(int[][] olis, int xvalue, int yvalue) {

    }

    public static class MicroGuesser {

    }

    public static class OilFieldSet {
        private final boolean[][] oil;
        private int size;

        public OilFieldSet(List<Pair> ps) {
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
