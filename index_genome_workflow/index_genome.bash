#!/bin/bash
#PBS -N job_index_dict_interval_list
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-54.jar run \
/groups/cgsd/$USER/gatk-workflows/index_genome_workflow/index_genome.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/index_genome_workflow/index_genome.json