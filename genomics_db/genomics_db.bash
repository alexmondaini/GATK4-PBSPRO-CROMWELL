#!/bin/bash
#PBS -N job_genomics_db
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
gatk --java-options "-Xms50G -Xmx50G" GenomicsDBImport \
--genomicsdb-workspace-path /groups/cgsd/alexandre/gatk-workflows/genomics_db/genome_db/ \
--batch-size 50 \
-R /groups/cgsd/alexandre/gatk-workflows/src/ref_Homo38/Homo_sapiens_assembly38_chrHPV.fasta \
-L /groups/cgsd/alexandre/gatk-workflows/src/interval_list/sy/Basic_Core_xGen_MSI_TERT_HPV_EBV_targets_hg38_sorted_merge.interval_list \
--sample-name-map /groups/cgsd/alexandre/gatk-workflows/src/pon/my_panel_of_normals/vcf_from_mutect/cohort.sample_map \
--tmp-dir /groups/cgsd/alexandre/gatk-workflows/genomics_db/tmp_dir 

