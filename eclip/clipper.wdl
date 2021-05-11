version 1.0

workflow Call_Peaks {
    input {
        Array[File] samples
        File chr_size
    }
    scatter (sample in samples) {
        call Sort_Bam {
            input:
            sort_star_bam = sample
        }
        call Index {
            input:
            index_after_sort = Sort_Bam.result_name_sort
        }
        call Clipper {
            input:
            call_peak_bam = Index.result_bam,
            call_peak_bai = Index.result_bai
        }
        call Wigs {
            input:
            wigs_bam = Index.result_bam,
            wigs_bai = Index.result_bai,
            chr_size = chr_size
        }
    }
}

task Sort_Bam {
    input {
        File sort_star_bam
    }
    String sort_star_bam_from_hg19 = basename(sort_star_bam,'.bam') + '_sorted_again.bam'

    command <<<
    mkdir tmp
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools sort -o ~{sort_star_bam} > 'sample.bam'
    >>>
    runtime {
        cpu: 3
        memory: "7 GB"
    }
    output {
        File result_name_sort = stdout()
    }
 }

task Index {
     input {
        File index_after_sort 
     }
     String result_bai_index = basename(index_after_sort)
     command <<<
     ln ~{index_after_sort}
     source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
     conda activate stepbystep
     samtools index ~{index_after_sort} > ~{result_bai_index}
     >>>
     runtime {
         cpu: 3
         memory: "6 GB"
     }
     output {
         File result_bai = "${result_bai_index}"
         File result_bam = "${index_after_sort}"
     }
}

task Clipper {
    input {
        File call_peak_bam
        File call_peak_bai
    }
    String bed_peak_intervals = basename(call_peak_bam,'.bam') + '.bed'
    
    command <<<
    clipper \
    --species hg19 \
    --bam ~{call_peak_bam} \
    --outfile ~{bed_peak_intervals}
    >>>
    
    runtime {
        cpu: 8
        memory: "50 GB"
        docker: "brianyee/clipper@sha256:094ede2a0ee7a6f2c2e07f436a8b63486dc4a072dbccad136b7a450363ab1876"
    }

}

task Wigs {
    input {
        File wigs_bam
        File wigs_bai
        File chr_size
        String? direction
    }
    String bw_pos = basename(wigs_bam,'.bam') + '_norm_pos.bw'
    String bw_neg = basename(wigs_bam,'.bam') + '_norm_neg.bw'

    command <<<
    makebigwigfiles \
    --bw_pos  ~{bw_pos} \
    --bw_neg ~{bw_neg}  \
    --bam  ~{wigs_bam} \
    --genome ~{chr_size} \
    --direction ~{default="r" direction}
    >>>

    runtime {
        cpu: 4
        memory: "30 GB"
        docker: "brianyee/makebigwigfiles@sha256:8d67afc36e388aa12f1b6d2bed8ea3b6ddaa9ec4296a93d5fa9f31a5b1ff16d4"
    }

    output {
        File result_bw_pos = "${bw_pos}"
        File result_bw_neg = "${bw_neg}"
    }
}
