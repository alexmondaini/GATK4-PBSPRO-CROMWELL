version 1.0

workflow SelectVariants {
    input {
        Array[File] vcfs
    }
    scatter (vcf in vcfs) {
        call Select_INDELS {
            input:
            vcf = vcf
        }
    }
}

task Select_INDELS {
    input {
        File vcf
    }
    String output_vcf = basename(vcf)

    command {
        module load gatk
        gatk SelectVariants -V ~{vcf} \
        -O ~{output_vcf} \
        --select-type-to-include INDEL
    }
    runtime {
        cpu: 4
        memory: "8 GB"
    }
    output {
        File out_vcf = "${output_vcf}"
    }
}