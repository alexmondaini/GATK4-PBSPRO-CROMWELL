version 1.0

workflow ValidateVariants {
    input {
        File reference
        File reference_fai
        Array[File] samples
    }
    
    scatter (sample in samples) {
    call Validate {
        input:
        reference = reference,
        sample = sample,
        reference_fai = reference_fai
    }
    }
}

task Validate {
    input {
        File reference
        File reference_fai
        File sample
    }
    command {
        gatk --java-options "-Xms6G -Xmx6G" ValidateVariants \
        -R ~{reference} \
        -V ~{sample} \
        --validation-type-to-exclude ALL
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09"
        cpu: 6
        memory: "8 GB"
    }
    output {
        File result = stdout()
    }

}
