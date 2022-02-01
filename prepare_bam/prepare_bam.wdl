version 1.0

workflow PrepareBam {
    input {
    Array[File] bams
    }

    scatter (bam in bams) {
        call Sort {
            input:
            bam = bam
        }
    }

    scatter (bam in sorted_bams) {
        call AddReadGroup {
            input:
            bam = bam
        }
    }
    output {
    Array[File] sorted_bams = Sort.out
    }
}

task Sort {
    input {
        File bam
    }
    String output_bam = basename(bam,'.round2.sorted_STAR_hg19Aligned.out.bam') + "_sorted.bam"

    command {
    module load samtools
    samtools sort ~{bam} \
    -o ~{output_bam}
    }
    runtime {
        cpu: 8
        memory: "16 GB"
    }

    output {
        File out = "${output_bam}"
    }
}

task AddReadGroup {
    input {
        File bam
    }
    String output_read = basename(bam,"_sorted")

    command <<<
    module load java
    module load Picard
    java -jar /software/Picard/2.25.2/picard.jar AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=LENIENT \ 
    I=~{bam} \
    O=~{output_read} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=sample
    >>>
    runtime {
        cpu: 8
        memory: "16 GB"
    }
    output {
        File out_final = "${output_read}"
    }
}