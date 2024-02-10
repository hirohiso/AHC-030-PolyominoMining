import spock.lang.Specification

class MainTest extends Specification {

    def "test"(){
        given:
        def input = new File("./src/test/resources/sample/sample1.txt").readLines()
        def (String ns,String ms,String es) = input[0].split(/\s+/)
        def n = ns.toInteger()
        def m = ms.toInteger()
        def e = es.toDouble()
        println(n)
        println(m)
        println(e)


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
        outs.write(sb.toString().getBytes("UTF-8"))
        outs.flush()
        println(reader.readLine())

        then:
        0 == 0
    }
}
