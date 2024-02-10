import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.text.Format;
import java.util.*;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.IntPredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

        var max = 0;
        for (int i = 0; i < M; i++) {
            var d = fs.ni();
            for (int j = 0; j < d; j++) {
                var x = fs.ni();
                var y = fs.ni();
                max++;
            }
        }
        //確定マス
        var grid = new int[N][N];
        for (int i = 0; i < grid.length; i++) {
            Arrays.fill(grid[i], Integer.MAX_VALUE);
        }

        //推測フェーズ
        var guesstable = new int[N][N];
        for (int i = 0; i < N; i++) {
            var list = new LinkedList<Pair>();
            for (int j = 0; j < N; j++) {
                list.add(new Pair(i, j));
            }
            pw.print("q " + list.size() + " ");
            for (Pair p : list) {
                pw.print(p.a + " " + p.b + " ");
            }
            pw.println();
            pw.flush();
            var vs = fs.ni();
            for (int j = 0; j < N; j++) {
                guesstable[i][j] += vs;
            }
        }
        for (int i = 0; i < N; i++) {
            var list = new LinkedList<Pair>();
            for (int j = 0; j < N; j++) {
                list.add(new Pair(j, i));
            }
            pw.print("q " + list.size() + " ");
            for (Pair p : list) {
                pw.print(p.a + " " + p.b + " ");
            }
            pw.println();
            pw.flush();
            var vs = fs.ni();
            for (int j = 0; j < N; j++) {
                guesstable[j][i] += vs;
            }
        }
        var cnt = 0;
        var stack = new LinkedList<Pair>();
        while (true) {
            var i = -1;
            var j = -1;
            if (stack.size() == 0) {
                var target = -1;
                updateGuess(guesstable, grid);
                for (int k = 0; k < N; k++) {
                    for (int l = 0; l < N; l++) {
                        if (grid[k][l] != Integer.MAX_VALUE) {
                            continue;
                        }
                        if (target < guesstable[k][l]) {
                            target = guesstable[k][l];
                            i = k;
                            j = l;
                        }
                    }
                }
            } else {
                var p = stack.pollFirst();
                i = p.a;
                j = p.b;
            }

            if (grid[i][j] != Integer.MAX_VALUE) {
                continue;
            }

            pw.println("q 1 " + i + " " + j);
            pw.flush();
            var v = fs.ni();
            grid[i][j] = v;
            //guess更新
            for (int k = 0; k < N; k++) {
                if (Math.abs(i - k) <= 1) {
                    guesstable[k][j] += (v != 0 ? +v : -1);
                }
                if (Math.abs(j - k) <= 1) {
                    guesstable[i][k] += (v != 0 ? +v : -1);
                }
            }
            if (v != 0) {
                var d = new int[][]{{-1, 0}, {0, -1}, {1, 0}, {0, 1}};
                for (int k = 0; k < d.length; k++) {
                    var x = i + d[k][0];
                    var y = j + d[k][1];
                    if (0 <= x && x < N && 0 <= y && y < N) {
                        if (grid[x][y] != Integer.MAX_VALUE) {
                            continue;
                        }
                        stack.addLast(new Pair(x, y));
                    }
                }
            }
            cnt += v;
            if (cnt == max) {
                break;
            }
        }

        //回答フェーズ
        var list = new LinkedList<Pair>();
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (grid[i][j] > 0 && grid[i][j] != Integer.MAX_VALUE) {
                    list.add(new Pair(i, j));
                }
            }
        }
        pw.print("a " + list.size() + " ");
        for (Pair p : list) {
            pw.print(p.a + " " + p.b + " ");
        }
        pw.println();
    }

    private static void updateGuess(int[][] table, int[][] wakwak) {
        var d = new int[][]{{-1, 0}, {0, -1}, {1, 0}, {0, 1}};
        for (int i = 0; i < table.length; i++) {
            for (int j = 0; j < table[0].length; j++) {
                var cnt = 0;
                for (int k = 0; k < d.length; k++) {
                    var x = i + d[k][0];
                    var y = j + d[k][1];
                    if (0 <= x && x < wakwak.length && 0 <= y && y < wakwak[0].length) {
                        if (wakwak[x][y] == 0) {
                            cnt++;
                        }
                    }
                }
                table[i][j] = Math.max(0, table[i][j] - cnt);
            }
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
