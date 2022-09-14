#!/bin/bash
#PBS -N job_helloworld
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -k oed

module load java/11.0.9

java -Dconfig.file=/groups/cgsd/alexandre/GATK_workflows/test.conf \
-jar /groups/cgsd/alexandre/cromwell-83.jar \
run /groups/cgsd/alexandre/GATK_workflows/hello_world/hello_world.wdl \
--inputs /groups/cgsd/alexandre/GATK_workflows/hello_world/hello_world.json