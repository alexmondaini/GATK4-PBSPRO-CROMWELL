#!/bin/bash
#PBS -N job_haplotypecaller
#PBS -l walltime=40:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/gatk4-germline-snps-indels/haplotypecaller-gvcf-gatk4.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/gatk4-germline-snps-indels/haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.json