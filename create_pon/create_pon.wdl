version 1.0

struct BamSample {
    File bam
    File bam_index
}

workflow Create_Pon {
    input {
        Array[BamSample] samples
        File reference
        File ref_fai
        File dict_file
    }
    scatter (sample in samples) {
    call ThePon {
        input:
        bam = sample.bam,
        bam_index = sample.bam_index,
        reference = reference,
        ref_fai = ref_fai,
        dict_file = dict_file
    }
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
        gatk --java-options "-Xmx17G" Mutect2 \
        -R ${reference} \
        -I ${bam} \
        -O ${strip_bam_extension}.vcf.gz
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09"
        cpu: 4
        memory: "20 GB"
    }
    output {
        File pon_vcf_file = "${strip_bam_extension}.vcf.gz"
    }
}