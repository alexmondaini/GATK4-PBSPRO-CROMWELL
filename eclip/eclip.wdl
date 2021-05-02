version 1.0

struct FastaSamples {
    Pair[File,File]  pair
    String barcode
}



workflow Eclip {
    
    input {
        Array[FastaSamples] samples
    }
    
    scatter (sample in samples) {
    call CutAdapt {
        input:
        fastq_r1 = sample.pair.left,
        fastq_r2 = sample.pair.right,
        barcode = sample.barcode
    }
    }
}

task CutAdapt {
    
    input {
        File fastq_r1
        File fastq_r2
        String barcode
    }
    
    String left_r1 = basename(fastq_r1,'.gz')
    String right_r2 = basename(fastq_r2,'.gz')

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
    -a ~{barcode} \
    -g ~{barcode} \
    -A ~{barcode} \
    -G ~{barcode} \
    -o ~{left_r1} \
    -p ~{right_r2} \
    ~{fastq_r1} \
    ~{fastq_r2}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        Array[Pair[File,File]] output_fq = "round1.fastq"
    }

}

#task CutAdapt_round2 {
#    input {
#        File left_r1
#        File right_r2
#    }
#    String round2_left_r1 = 'round2' + left_r1
#    String round2_right_r2 = 'round2' + right_r2
#
#    command <<<
#    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
#    conda activate stepbystep
#    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
#    -a ~{barcode} \
#    -g ~{barcode} \
#    -A ~{barcode} \
#    -G ~{barcode} \
#    -o ~{round2_left_r1} \
#    -p ~{round2_right_r2} \
#    ~{left_r1} \
#    ~{right_r2}
#    >>>
#
#    runtime {
#        cpu: 3
#        memory: "6 GB"
#    }
#
#    output {
#        File round2_r1 = "${round2_left_r1}"
#        File round2_r2 = "${round2_right_r2}"
#    }
#}