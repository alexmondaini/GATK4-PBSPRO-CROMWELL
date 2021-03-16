version 1.0

workflow HelloWorld {

    input {
        File here
    }

    call WriteGreeting {
        input:
        here = here
    }
}

task WriteGreeting {
    input {
        File here
    }
    command {
        cat ~{here}
    }
    output {
        File output_greeting = stdout()
    }
}
