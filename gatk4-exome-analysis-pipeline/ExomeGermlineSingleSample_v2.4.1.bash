#!/bin/bash
#PBS -N job.targetsequencing
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/gatk4-exome-analysis-pipeline/ExomeGermlineSingleSample_v2.4.1.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/gatk4-exome-analysis-pipeline/ExomeGermlineSingleSample_v2.4.1.AML.inputs.json