version 1.0

workflow HelloWorld_Scatter {

    input {
        Array[File] collection
    }

    scatter (item in collection) {
    call WriteGreeting {
        input:
        item = item
    }
    }
}

task WriteGreeting {
    input {
        File item
    }
    command {
        cat ~{item}
    }
    output {
        File output_greeting = stdout()
    }
}
