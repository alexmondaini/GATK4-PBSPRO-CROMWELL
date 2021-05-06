version 1.0

struct FastaSamples {
    File fastq_r1
    File fastq_r2
    String barcode
}

workflow Eclip {
    
    input {
        Array[FastaSamples] samples
    }
    
    scatter (sample in samples) {
    call CutAdapt {
        input:
        fastq_r1 = sample.fastq_r1,
        fastq_r2 = sample.fastq_r2,
        barcode = sample.barcode
    }
    call FasQC_round1 {
        input:
        fastqc_r1 = CutAdapt.result_round1_cutadapt_left,
        fastqc_r2 = CutAdapt.result_round1_cutadapt_right
    }
    call CutAdapt_round2 {
        input:
        round1_left_r1 = CutAdapt.result_round1_cutadapt_left,
        round1_right_r2 = CutAdapt.result_round1_cutadapt_right,
        barcode = sample.barcode
    }
    call FasQC_round2 {
        input:
        fastqc_round2_r1 = CutAdapt_round2.result_round2_cutadapt_left,
        fastqc_round2_r2 = CutAdapt_round2.result_round2_cutadapt_right
    }
    call FastQ_sort {
        input:
        fastq_sort_r1 = CutAdapt_round2.result_round2_cutadapt_left,
        fastq_sort_r2 = CutAdapt_round2.result_round2_cutadapt_right
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
        File result_round1_cutadapt_left = "${left_r1}"
        File result_round1_cutadapt_right = "${right_r2}"
    }

}

task FasQC_round1 {
    input {
        File fastqc_r1
        File fastqc_r2
    }
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastqc -t 2 --extract -k 7 ~{fastqc_r1} -o .
    fastqc -t 2 --extract -k 7 ~{fastqc_r2} -o .
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
}

task CutAdapt_round2 {
    input {
        File round1_left_r1
        File round1_right_r2
        String barcode
    }

    String round2_left_r1 = basename(round1_left_r1,'fq') + 'round2.fq'
    String round2_right_r2 = basename(round1_right_r2,'fq') + 'round2.fq'

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
    -a ~{barcode} \
    -g ~{barcode} \
    -A ~{barcode} \
    -G ~{barcode} \
    -o ~{round2_left_r1} \
    -p ~{round2_right_r2} \
    ~{round1_left_r1} \
    ~{round1_right_r2}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        File result_round2_cutadapt_left = "${round2_left_r1}" 
        File result_round2_cutadapt_right = "${round2_right_r2}"
    }
}

task FasQC_round2 {
    input {
        File fastqc_round2_r1
        File fastqc_round2_r2
    }
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastqc -t 2 --extract -k 7 ~{fastqc_round2_r1} -o .
    fastqc -t 2 --extract -k 7 ~{fastqc_round2_r2} -o .
    >>>
    
    runtime {
        cpu: 3
        memory: "5 GB"
    }

}

task FastQ_sort {
    input {
        File fastq_sort_r1
        File fastq_sort_r2
    }
    String sorted_r1 = basename(fastq_sort_r1,'.fq') + 'sorted.fq'
    String sorted_r2 = basename(fastq_sort_r2,'.fq') + 'sorted.fq' 

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastq-sort --id ~{fastq_sort_r1} > ~{sorted_r1}
    fastq-sort --id ~{fastq_sort_r2} > ~{sorted_r2}
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
    output {
        File result_fastq_sort_left = "${sorted_r1}"
        File result_fastq_sort_right = "${sorted_r2}"
     }
}