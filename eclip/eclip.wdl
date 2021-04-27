version 1.0

workflow Eclip {
    
    input {
        Array[Pair[File,File]]  fastqs
    }
    
    scatter (fast_pair in fastqs) {
    
    call CutAdapt {
        input:
        fastq_r1 = fast_pair.left,
        fastq_r2 = fast_pair.right
    }
    }
}

task CutAdapt {
    
    input {
        File fastq_r1
        File fastq_r2
    }
    
    String left_fasta = basename(fastq_r1,'.gz')
    String right_fasta = basename(fastq_r2,'.gz')

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
    -a ATCACG \
    -g ATCACG \
    -A ATCACG \
    -G ATCACG \
    -o ~{left_fasta} \
    -p ~{right_fasta} \
    ~{fastq_r1} \
    ~{fastq_r2}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        File read_left_r1 = "${left_fasta}"
        File read_right_2 = "${right_fasta}"
    }

}