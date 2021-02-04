#!/bin/bash
#PBS -N job.targetsequencing
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=2:mem=6gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Xmx25g -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/tmp/ExomeGermlineSingleSample_v2.4.1.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/tmp/ExomeGermlineSingleSample_v2.4.1.inputs.json