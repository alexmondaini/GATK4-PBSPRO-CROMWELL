version 1.0

workflow Create_Pon {
    input {
        File bam
        File bam_index
        File reference
        File ref_fai
        File dict_file
    }
    call ThePon {
        input:
        bam = bam,
        bam_index = bam_index,
        reference = reference,
        ref_fai = ref_fai,
        dict_file = dict_file
    }
}

task ThePon {
    input {
        File bam
        File bam_index
        File reference
        File ref_fai
        File dict_file
    }

    
    String strip_bam_extension = basename(bam,".bam")

    command {
        /gatk/gatk Mutect2 \
        -R ${reference} \
        -I ${bam} \
        -O ${strip_bam_extension}.vcf.gz
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09"
        cpu: 4
        memory: "16 GB"
    }
    output {
        File pon_vcf_file = "${strip_bam_extension}.vcf.gz"
    }
}