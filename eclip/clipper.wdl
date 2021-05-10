version 1.0

workflow Call_Peaks {
    input {
        Array[File] samples
        File chr_size
    }
    scatter (sample in samples) {
        call Clipper {
            input:
            call_peak_bam = sample
        }
        call Wigs {
            input:
            bam = sample,
            chr_size = chr_size
        }
    }
}

task Clipper {
    input {
        File call_peak_bam
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
        File bam
        File chr_size
        String? direction
    }
    String bw_pos = basename(bam,'.bam') + '_norm_pos.bw'
    String bw_neg = basename(bam,'.bam') + '_norm_neg.bw'

    command <<<
    makebigwigfiles \
    --bw_pos  ~{bw_pos} \
    --bw_neg ~{bw_neg}  \
    --bam  ~{bam} \
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