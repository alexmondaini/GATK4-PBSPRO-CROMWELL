version 1.0

workflow Test {
    call Clipper {

    }

}


task Clipper {

    command <<<
    clipper -h
    >>>
    runtime {
        cpu: 4
        memory: "6 GB"
        docker: "brianyee/clipper@sha256:094ede2a0ee7a6f2c2e07f436a8b63486dc4a072dbccad136b7a450363ab1876"
    }
    output {
        String hello = stdout()
    }

}

