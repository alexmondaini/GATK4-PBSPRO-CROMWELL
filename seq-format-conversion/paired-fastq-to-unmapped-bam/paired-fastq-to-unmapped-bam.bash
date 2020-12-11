#!/bin/bash
#PBS -N job_paired-fastq-to-unmapped-bam
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-53.1.jar run \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/paired-fastq-to-unmapped-bam/example.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/paired-fastq-to-unmapped-bam/example.json