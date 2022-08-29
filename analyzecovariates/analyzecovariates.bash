#!/bin/bash
#PBS -N job_analyzecovariates
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/analyzecovariates/analyzecovariates.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/analyzecovariates/analyzecovariates.json