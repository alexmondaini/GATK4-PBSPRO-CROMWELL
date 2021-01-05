#!/bin/bash
#PBS -N interleaved-fastq-to-paired-fastq
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/interleaved-fastq-to-paired-fastq/interleaved-fastq-to-paired-fastq.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/interleaved-fastq-to-paired-fastq/interleaved-fastq-to-paired-fastq.inputs.json