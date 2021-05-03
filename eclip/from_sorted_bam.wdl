version 1.0

workflow Eclip {
    
    input {
        Array[File] bams
    }
    
    scatter (bam in bams) {
        call Index {
            input:
            bam = bam
        }
    }
    #call Merge {
    #    input:
    #    pairs = bam.result
    #}
}

task Index {
    input {
        File bam
    }

    String bai = basename(bam,'.bam') + ".bai"

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    ln ~{bam}
    samtools index ~{bam} > ~{bai}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        File result = stdout()
    }

}

task Merge {
    input {
        Array[Pair[File,File]] pairs
    }
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    samtools merge ~{pairs.left} ~{pairs.right}
    >>>

    runtime {
        cpu: 3
        memory: "8 GB"
    }

}
