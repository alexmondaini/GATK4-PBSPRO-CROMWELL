#!/bin/bash
#PBS -N job_parallel_paired-fastq-to-unmapped-bam
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=6gb
#PBS -q cgsd
#PBS -v USER
#PBS -k oed 

module load java/13.0.2
CROMWELL_VERSION='83'
cd $PBS_O_WORKDIR

java -Dconfig.file=/groups/cgsd/$USER/GATK_workflows/application.conf -Xmx5g \
-jar /groups/cgsd/$USER/cromwell-${CROMWELL_VERSION}.jar run \
parallel-paired-fastq-ubam.wdl --inputs liver.json