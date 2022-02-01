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

    scatter (bam_sorted in Sort.out) {
        call AddReadGroup {
            input:
            bam_sorted = bam_sorted
        }
    }
    output {
    Array[File] sorted_bams = Sort.out
    Array[File] final_bam_read_group = AddReadGroup.out_final
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
        File bam_sorted
    }
    
    String output_read = sub(bam_sorted,"_sorted","")

    command {
    module load java/11.0.9
    module load Picard
    java -jar /software/Picard/2.25.2/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=~{bam_sorted} \
    O=~{output_read} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=sample
    }

    runtime {
        cpu: 8
        memory: "16 GB"
    }
    output {
        File out_final = "${output_read}"
    }
}