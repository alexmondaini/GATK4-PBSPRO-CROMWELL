version 1.0 

workflow Test {
    input {
        File chr_len
        File chr_name_len
        File chr_name
        File chr_start
        File genome
        File genome_parameters
        File sa
        File sa_index
        File smallrep_fasta
        File fastq_starrep_r1
        File fastq_starrep_r2
    }
    call Star {
        input:
        chr_len = chr_len,
        chr_name_len = chr_name_len,
        chr_name = chr_name,
        chr_start = chr_start,
        genome = genome,
        genome_parameters = genome_parameters,
        sa = sa,
        sa_index = sa_index,
        smallrep_fasta = smallrep_fasta,
        fastq_starrep_r1 = fastq_starrep_r1,
        fastq_starrep_r2 =fastq_starrep_r2
        }
    }

task Star {
    input {
        File chr_len
        File chr_name_len
        File chr_name
        File chr_start
        File genome
        File genome_parameters
        File sa
        File sa_index
        File smallrep_fasta
        File fastq_starrep_r1
        File fastq_starrep_r2
    }

    String prefix = basename(fastq_starrep_r1,'_r1.fq') + "_STAR"

    command <<<
    mkdir RepElements
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    STAR \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir RepElements \
    --genomeLoad NoSharedMemory \
    --alignEndsType EndToEnd \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 30 \
    --outFilterMultimapScoreRange 1 \
    --outFileNamePrefix ~{prefix} \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --outBAMcompression 10 \
    --outReadsUnmapped Fastx \
    --outFilterScoreMin 10 \
    --outSAMattrRGline ID:foo \
    --outSAMattributes All \
    --outSAMmode Full \
    --outStd Log \
    --readFilesIn ~{fastq_starrep_r1} ~{fastq_starrep_r2}
    >>>
    runtime {
        cpu: 1
        memory: "2 GB"
    }
    output {
        File result_fq = glob('*fq')
        File result_bam = glob('*.bam')
    }
}