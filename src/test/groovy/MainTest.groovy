import spock.lang.Specification

class MainTest extends Specification {

    def "test"() {
        given:
        def input = new File("./src/test/resources/sample/sample1.txt").readLines()
        def (String ns, String ms, String es) = input[0].split(/\s+/)
        def n = ns.toInteger()
        def m = ms.toInteger()
        def e = es.toDouble()
        println(n)
        println(m)
        println(e)

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
        while ((str = reader.readLine()) != null) {
            println(str)
            def token = str.split(" ")
            if (token[0] == "q") {
                if (token[1].toInteger() == 1) {
                    //採掘
                    def x = token[2].toInteger()
                    def y = token[3].toInteger()
                    def temp = new StringBuilder()
                    temp.append(grid[x][y]).append("\n")

                    outs.write(temp.toString().getBytes("UTF-8"))
                    outs.flush()
                } else {
                    //占い
                    def p = normalize[countNo]
                    def size = token[1].toInteger()
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
                    def vs = Math.max(0 ,Math.round(p * vari + mean));
                    temp.append(vs).append("\n")
                    outs.write(temp.toString().getBytes("UTF-8"))
                    outs.flush()
                    countNo++;
                }

            } else {
                //回答
                //todo: 答え合わせ
                break;
            }
        }

        then:
        0 == 0
    }
}
