version 1.0

workflow HelloWorld {

    input {
        File change
        File constant
    }

    call WriteGreeting {
        input:
        change = change,
        constant = constant
    }
}

task WriteGreeting {
    input {
        File change
        File constant
    }
    command {
        cat ~{change}
    }
    output {
        File output_greeting = stdout()
    }
}
