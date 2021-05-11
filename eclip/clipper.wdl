version 1.0

workflow Call_Peaks {
    input {
        Array[File] samples
        File chr_size
    }
    scatter (sample in samples) {
        
        String output_name = basename(sample,'.sorted_STAR_hg19Aligned.out.bam')

        call Sort_and_Index_Bam {
            input:
            sort_star_bam = sample,
            result_bam = output_name + '.bam',
            result_view = output_name + '_final.bam'
        }
        call Clipper {
            input:
            call_peak_bam = Sort_and_Index_Bam.result_sorted_indexed_bam,
            call_peak_bai = Sort_and_Index_Bam.result_sorted_indexed_bai
        }
        call Wigs {
            input:
            wigs_bam = Sort_and_Index_Bam.result_sorted_indexed_bam,
            wigs_bai = Sort_and_Index_Bam.result_sorted_indexed_bai,
            chr_size = chr_size
        }
    }
}

task Sort_and_Index_Bam {
    input {
        File sort_star_bam
        String result_bam
        String result_view
    }

    command <<<
    set -e
    ln ~{sort_star_bam} ~{basename(sort_star_bam)}
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools sort -o ~{result_bam} ~{basename(sort_star_bam)} 
    samtools index ~{result_bam}
    samtools view -f 128 -b ~{result_bam}
    >>>

    runtime {
        cpu: 3
        memory: "7 GB"
    }
    output {
        File result_sorted_indexed_bam = result_bam
        File result_sorted_indexed_bai = "~{result_bam}.bai"
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
        cpu: 20
        memory: "30 GB"
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
        cpu: 20
        memory: "30 GB"
        docker: "brianyee/makebigwigfiles@sha256:8d67afc36e388aa12f1b6d2bed8ea3b6ddaa9ec4296a93d5fa9f31a5b1ff16d4"
    }

    output {
        File result_bw_pos = "${bw_pos}"
        File result_bw_neg = "${bw_neg}"
    }
}
