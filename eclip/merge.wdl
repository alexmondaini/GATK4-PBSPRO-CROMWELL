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
    scatter (merged_bam in Merge.result) {

        call Index {
            input:
            merged_bam = merged_bam
        }
        call View {
            input:
            final_bam = Index.result_bam,
            final_bai = Index.result_bai
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
    String merged_bam = basename(bam_rep_1,'Aligned.sortedByCoord.out.bam') + "_" + basename(bam_rep_2,'Aligned.sortedByCoord.out.bam') + ".bam"
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools merge ~{merged_bam} ~{bam_rep_1} ~{bam_rep_2}
    >>>
    runtime {
        cpu: 4
        memory: "15 GB"
    }
    output {
        File result = "${merged_bam}"
    }
}

task Index {
    input {
        File merged_bam
    }
    
    String merged_bai = basename(merged_bam,".bam") + '.bai'
    
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools index ~{merged_bam} > ~{merged_bai}
    >>>
    runtime {
        cpu: 3
        memory: "8 GB"
    }
    output {
        File result_bai = "${merged_bai}"
        File result_bam = "${merged_bam}"
    }
}

task View {
    input {
        File final_bam
        File final_bai
    }
    String ready_to_peak_call = basename(final_bam,'.bam') + '_ready.bam'
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools view -f 128 -b -o ~{ready_to_peak_call} ~{final_bam}
    >>>
    runtime {
        cpu: 3
        memory: "8 GB"
    }
    output {
        File ready_bam = "${ready_to_peak_call}"
    }
}