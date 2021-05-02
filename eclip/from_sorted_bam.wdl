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
}

task Index {
    input {
        File bam
    }

    String bai = basename(bam,'bam') + ".bai"

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
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