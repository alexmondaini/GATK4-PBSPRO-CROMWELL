#!/bin/bash
#PBS -N to_unmapped_bam
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb:host=compute-0-7
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-53.1.jar run \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/bam-to-unmapped-bams/bam-to-unmapped-bams.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/seq-format-conversion/bam-to-unmapped-bams/bam-to-unmapped-bams.inputs.json