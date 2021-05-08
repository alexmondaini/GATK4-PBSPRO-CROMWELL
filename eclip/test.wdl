version 1.0 

workflow Test {
    input {
        File zipped_star_files
        File fastq_starrep_r1
        File fastq_starrep_r2
    }
    call Star {
        input:
        zipped_star_files = zipped_star_files,
        fastq_starrep_r1 = fastq_starrep_r1,
        fastq_starrep_r2 =fastq_starrep_r2
        }
    }

task Star {
    input {
        File zipped_star_files
        File fastq_starrep_r1
        File fastq_starrep_r2
    }

    String prefix = sub(basename(fastq_starrep_r1,'.fq.gz'),'_r1','') + "_STAR"

    command <<<
    mkdir RepElements
    tar -xzf ~{zipped_star_files} -C RepElements
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
    --readFilesCommand zcat \
    --outStd Log \
    --readFilesIn ~{fastq_starrep_r1} ~{fastq_starrep_r2}
    >>>
    runtime {
        cpu: 1
        memory: "2 GB"
    }
    output {
        File result_fq_r1 = "${prefix}Unmapped.out.mate1"
        File result_fq_r2 = "${prefix}Unmapped.out.mate2"
        File result_bam = "${prefix}Aligned.out.bam"
    }
}