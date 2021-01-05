#!/bin/bash
#PBS -N validatebam
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/seq-format-validation/validate-bam.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/seq-format-validation/validate-bam.json
