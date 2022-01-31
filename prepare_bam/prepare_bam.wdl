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
    Array[File] sorted_bams = Sort.out
    Array[File] final_bams_with_read_group = AddReadGroup.out
}

task Sort {
    input {
        File bam
    }
    String output_bam = basename(bam,'.sorted_STAR_hg19Aligned.out.bam') + "_sorted.bam"

    command {
    module load samtools
    samtools sort ${bam} \
    -o ${output_bam}
    }

    output {
        File out = "${output_bam}"
    }
}

task AddReadGroup {
    input {
        File bam
    }
    String output_bam = basename(bam,"_sorted.bam") + ".bam"

    command <<<
    modue load java
    module load Picard
    java -jar /software/Picard/2.25.2/picard.jar AddOrReplaceReadGroups \
    -I=${bam} \
    -O=${output_bam} \
    RGPL=ILLUMINA
    >>>

    output {
        File out = "${output_bam}"
    }
}