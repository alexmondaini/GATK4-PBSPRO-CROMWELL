#!/bin/bash
#PBS -N job_create_pon
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/db.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/create_pon/create_pon.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/create_pon/test.json