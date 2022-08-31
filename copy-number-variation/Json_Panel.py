from pathlib import Path
import json

bams = sorted(Path('/groups/cgsd/alexandre/liver/bams/').glob('*-1.hg38.bam'))
bais = sorted(Path('/groups/cgsd/alexandre/liver/bams/').glob('*-1.hg38.bai'))
print(bams)

data = {
  "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.2.0.0",

  "CNVSomaticPanelWorkflow.normal_bams": [str(x) for x in bams],
  "CNVSomaticPanelWorkflow.normal_bais": [str(x) for x in bais],
  
  "CNVSomaticPanelWorkflow.ref_fasta": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta",
  "CNVSomaticPanelWorkflow.ref_fasta_fai": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta.fai",
  "CNVSomaticPanelWorkflow.ref_fasta_dict": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.dict",
  
  "CNVSomaticPanelWorkflow.pon_entity_id":"Liver",

  "CNVSomaticPanelWorkflow.padding": "100",
  "CNVSomaticPanelWorkflow.bin_length": "0",
  "CNVSomaticPanelWorkflow.do_explicit_gc_correction": "true",
  "CNVSomaticPanelWorkflow.mem_gb_for_annotate_intervals": "6",
  
  
  "CNVSomaticPanelWorkflow.segmental_duplication_track_bed": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/segmental_dup.sorted.merged.bed.gz",
  "CNVSomaticPanelWorkflow.segmental_duplication_track_bed_idx": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/segmental_dup.sorted.merged.bed.gz.tbi",
  
  "CNVSomaticPanelWorkflow.mappability_track_bed": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/k36.umap.bed.gz",
  "CNVSomaticPanelWorkflow.mappability_track_bed_idx": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/k36.umap.bed.gz.tbi",
  
  "CNVSomaticPanelWorkflow.intervals": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/S07604715_Regions.interval_list",
  "CNVSomaticPanelWorkflow.blacklist_intervals": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/CNV_and_centromere_blacklist.hg38liftover.list"
}

if __name__=='__main__':
    output = Path('LiverPanel.json')
    with output.open('w') as f:
        json.dump(data,f,indent=4)
