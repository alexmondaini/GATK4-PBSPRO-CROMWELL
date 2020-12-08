#!/bin/bash
#PBS -N job_mutect2
#PBS -l walltime=40:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-53.1.jar run \
/groups/cgsd/$USER/gatk-workflows/mutect2/mutect2.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/mutect2/mutect2.inputs.json