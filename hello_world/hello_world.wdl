version 1.0

workflow HelloWorld {
    call WriteGreeting
}

task WriteGreeting {
    command {
        echo "Hello World"
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:33574f446ac991f77bac125fbf6a2340e6db972a3f334e6c61bff94740165938"
        memory: "3 GB"
        cpu: 4
    }
    output {
        File output_greeting = stdout()
    }
}
