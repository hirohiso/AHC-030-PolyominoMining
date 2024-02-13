import spock.lang.Specification

class MainTest extends Specification {

    def "test"() {
        given:
        def input = new File(inFile).readLines()
        def result = new File(inFile.replace('/in/', '/out/'))
        result.parentFile.mkdirs()
        result.text = ""
        def (String ns, String ms, String es) = input[0].split(/\s+/)
        def n = ns.toInteger()
        def m = ms.toInteger()
        def e = es.toDouble()

        //読み込み
        int[][] grid = input[2 * m + 1..2 * m + n].collect {
            it -> it.split(' ').collect { it as int }
        } as int[][]
        double[] normalize = input[2 * m + n + 1..<input.size()].collect { it -> it.toDouble() } as double[]

        def outs = new PipedOutputStream()
        def ins = new PipedInputStream()
        Thread thread = new Thread({
            Main.solve(new PipedInputStream(outs), new PrintStream(new PipedOutputStream(ins)))
        } as Runnable)

        when:
        thread.start()
        sleep(50)
        def reader = ins.newReader()
        def sb = new StringBuilder()
        for (i in 0..m) {
            sb.append(input[i]).append("\n")
        }
        //初期入力
        outs.write(sb.toString().getBytes("UTF-8"))
        outs.flush()
        //対話開始
        String str;
        def countNo = 0;

        def cost = 0.0;
        def totalmine = 0;
        def zeromine = 0;

        def actual = true;
        while ((str = reader.readLine()) != null) {
            result << (str + "\n")
            def token = str.split(" ")
            if (token[0] == "q") {
                if (token[1].toInteger() == 1) {
                    cost += 1
                    //採掘
                    totalmine++
                    def x = token[2].toInteger()
                    def y = token[3].toInteger()
                    def temp = new StringBuilder()
                    temp.append(grid[x][y]).append("\n")
                    outs.write(temp.toString().getBytes("UTF-8"))
                    outs.flush()
                    if (grid[x][y] == 0) {
                        zeromine++;
                    }
                } else {
                    //占い
                    def p = normalize[countNo]
                    def size = token[1].toInteger()
                    cost += (1 / Math.sqrt(size))
                    def sumv = 0;
                    //v(i,j)の合計値を求める
                    for (i in 0..<size) {
                        def x = token[2 * i + 2].toInteger()
                        def y = token[2 * i + 1 + 2].toInteger()
                        sumv += grid[x][y]
                    }
                    //
                    def temp = new StringBuilder()
                    def mean = (size - sumv) * e + sumv * (1 - e)
                    def vari = size * (1 - e) * e
                    def vs = Math.max(0, Math.round(p * vari + mean));
                    temp.append(vs).append("\n")
                    outs.write(temp.toString().getBytes("UTF-8"))
                    outs.flush()
                    countNo++;
                }

            } else {
                //回答
                actual = valid(str, grid)
                if (actual) {
                    outs.write("1\n".getBytes("UTF-8"))
                    outs.flush()
                    //todo: 答え合わせ
                    //全体的な失敗率やコスト改善率の平均、分散などの統計値が知りたい
                    //コスト改善率　= 実際のコスト /全部めくった時のコスト(N^2)
                    print result.name + "," + (n * n) + "," + m + "," + grid.collect { it -> it.sum() }.sum() + "," + cost + "," + zeromine + "," + totalmine
                    println()
                    break;
                } else {
                    outs.write("0\n".getBytes("UTF-8"))
                    outs.flush()
                    cost += 1
                }
            }
        }

        then:
        //actual
        thread.join()

        where:
        inFile << getTxtFilesFromFolder("./src/test/resources/in/")
    }

    //answerがgridと一致しているか判定する
    static boolean valid(String answer, int[][] grid) {
        def tokens = answer.split(" ")
        def size = tokens[1].toInteger()
        def result = true

        def gridsize = 0
        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[i].length; j++) {
                if (grid[i][j] > 0) {
                    gridsize++
                }
            }
        }
        result &= (gridsize == size)
        for (i in 0..<size) {
            def x = tokens[2 * i + 2].toInteger()
            def y = tokens[2 * i + 1 + 2].toInteger()
            result &= (grid[x][y] > 0)
        }
        return result
    }

    // 指定されたフォルダから.txtファイル一覧を取得するヘルパーメソッド
    static List<String> getTxtFilesFromFolder(String folderPath) {
        return new File(folderPath).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".txt")
            }
        })*.path.sort()[0..20] //todo:一旦絞る
    }
}
