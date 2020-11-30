version 1.0

workflow HelloWorld {
    call WriteGreeting
}

task WriteGreeting {
    command {
        echo "Hello World"
    }
  #  runtime {
  #      docker: "broadinstitute/gatk:latest"
  #      memory: "3 GB"
  #      cpu: 4
  #  }
    output {
        File output_greeting = stdout()
    }
}
