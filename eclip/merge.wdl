version 1.0

struct ToMerge {
    File bam_rep_1
    File bai_rep_1
    File bam_rep_2
    File bai_rep_2
}

workflow SamtoolsMerge {
    input {
        Array[ToMerge] samples
    }
    scatter (sample in samples) {
        call Merge {
            input:
            bam_rep_1 = sample.bam_rep_1,
            bai_rep_1 = sample.bai_rep_1,
            bam_rep_2 = sample.bam_rep_2,
            bai_rep_2 = sample.bai_rep_2
        }
    }
}

task Merge {
    input {
        File bam_rep_1
        File bai_rep_1
        File bam_rep_2
        File bai_rep_2
    }
    String merged_bam = "merged.bam"
    command <<<
    samtools merge ~{merged_bam} ~{bam_rep_1} ~{bam_rep_2}
    >>>
    runtime {
        cpu: 6
        memory: "10 GB"
    }
    output {
        File result = "${merged_bam}"
    }
}