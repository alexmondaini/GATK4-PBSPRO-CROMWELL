#!/bin/bash
#PBS -N job_BWAindex_GATKDict
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q cgsd
#PBS -v USER
module load java
java -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf \
-jar /groups/cgsd/$USER/cromwell-53.1.jar run \
/groups/cgsd/$USER/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/nick_data/index_genome_workflow/index_genome.wdl \
--inputs \
/groups/cgsd/$USER/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/nick_data/index_genome_workflow/index_genome.json
