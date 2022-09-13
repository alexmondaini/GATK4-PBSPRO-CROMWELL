#!/bin/bash
#PBS -N job_helloworld
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd

module load java/11.0.9


java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/application.conf \
-jar /groups/cgsd/$USER/cromwell-83.jar \
run /groups/cgsd/$USER/gatk-workflows/hello_world/hello.wdl