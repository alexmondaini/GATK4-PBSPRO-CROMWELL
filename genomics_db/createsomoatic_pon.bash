#!/bin/bash
#PBS -N job_create_pon
#PBS -l walltime=240:00:00
#PBS -l select=1:ncpus=12:mem=60gb
#PBS -q cgsd
#PBS -v USER
cd $PBS_O_WORKDIR
export SINGULARITY_CACHEDIR=/groups/cgsd/$USER/.singularity
module load singularity
module load java
singularity exec -C --bind /groups/cgsd/alexandre/gatk-workflows/ \
docker://broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09 \
gatk --java-options "-Xms50G -Xmx50G" CreateSomaticPanelOfNormals \
-R /groups/cgsd/alexandre/gatk-workflows/src/ref_Homo38/Homo_sapiens_assembly38_chrHPV.fasta \
-V gendb:///groups/cgsd/alexandre/gatk-workflows/genomics_db/genome_db \
--tmp-dir /groups/cgsd/alexandre/gatk-workflows/genomics_db/pon_tmp_dir \
-O /groups/cgsd/alexandre/gatk-workflows/genomics_db/SY_pon.vcf.gz