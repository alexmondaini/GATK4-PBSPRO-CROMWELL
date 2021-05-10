version 1.0

workflow Tosort {
    input {
        File sort_star_bam
    }
    call Sort_Bam {
        input:
        sort_star_bam = sort_star_bam
    }
}



task Sort_Bam {
    input {
        File sort_star_bam
    }
    String sort_star_bam_from_hg19 = sub(basename(sort_star_bam),'sorted_STAR_hg19_Aligned\\.out','')

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools sort -o ~{sort_star_bam_from_hg19}  -T " foo"  ~{sort_star_bam} 
    >>>
    runtime {
        cpu: 3
        memory: "7 GB"
    }
    output {
        File result_name_sort = "${sort_star_bam_from_hg19}"
    }
 }