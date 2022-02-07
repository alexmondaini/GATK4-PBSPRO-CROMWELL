version 1.0

workflow SelectVariants {
    input {
        Array[Pair[File,File]] vcf_list
    }
    scatter (pair in vcf_list) {
        call Select_INDELS {
            input:
            vcf = pair.left,
            vcf_index = pair.right
        }
    }
}

task Select_INDELS {
    input {
        File vcf
        File vcf_index
    }
    String output_vcf = basename(vcf)


    command {
        module load java/11.0.9
        module load gatk
        gatk SelectVariants -V ~{vcf} -O ~{output_vcf} --select-type-to-include INDEL
    }
    runtime {
        cpu: 4
        memory: "8 GB"
    }
    output {
        File out_vcf = "${output_vcf}"
    }
}