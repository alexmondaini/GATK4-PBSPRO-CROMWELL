from pathlib import Path
import json

bams = sorted(Path('/groups/cgsd/alexandre/cromwell-executions/PreProcessingForVariantDiscovery_GATK4/').glob('*/call-GatherBamFiles/execution/*ba?'))
print(bams)

def pairwise(iterator):
    a = iter(iterator)
    return zip(a,a)

normals = list(pairwise([x for x in bams if x.name.split('.')[0].split('-')[1] == '1']))
tumors  = list(pairwise([x for x in bams if x.name.split('.')[0].split('-')[1] != '1']))

# tumor_samples = [t[1].stem.split('.')[0] for t in tumors]
# # check done samples
# p = sorted(Path('/groups/cgsd/alexandre/liver/mutect2_filter_artifacts_vcf/').glob('*vcf'))
# done_samples = [x.stem.split('.')[0] for x in p]


data = []
for n in normals:
    for t in tumors:
        if n[1].stem.split('.')[0].split('-')[0] == t[1].stem.split('.')[0].split('-')[0]:
            # if t[1].stem.split('.')[0] not in done_samples:
            data.append((n[1],n[0],t[1],t[0]))


def create_data():
    result = []
    for file in data:

        result.append(
            {
    "Mutect2.gatk_docker": "broadinstitute/gatk:4.2.0.0",
    "Mutect2.scatter_count": 30,
    "Mutect2.filter_funcotations": "true",
    "Mutect2.run_funcotator": "true",
    "Mutect2.make_bamout": "true",
    "Mutect2.m2_extra_args": "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",
    
    "Mutect2.intervals": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/S07604715_Regions.interval_list",
    
    "Mutect2.funco_reference_version": "hg38",
    "Mutect2.funco_data_sources_tar_gz": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/funcotator_dataSources.v1.6.20190124s.tar.gz",
    "Mutect2.funco_transcript_selection_list": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt",
    
  
    "Mutect2.ref_fasta": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta",
    "Mutect2.ref_fai": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta.fai",
    "Mutect2.ref_dict": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.dict",
    "Mutect2.normal_reads": f"{file[0]}",
    "Mutect2.normal_reads_index": f"{file[1]}",
    "Mutect2.tumor_reads": f"{file[2]}",
    "Mutect2.tumor_reads_index": f"{file[3]}",
  
    "Mutect2.pon": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/1000g_pon.hg38.vcf.gz",
    "Mutect2.pon_idx": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/1000g_pon.hg38.vcf.gz.tbi",
    "Mutect2.gnomad": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/af-only-gnomad.hg38.vcf.gz",
    "Mutect2.gnomad_idx": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/af-only-gnomad.hg38.vcf.gz.tbi",
    "Mutect2.variants_for_contamination": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/small_exac_common_3.hg38.vcf.gz",
    "Mutect2.variants_for_contamination_idx": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/small_exac_common_3.hg38.vcf.gz.tbi",
    "Mutect2.realignment_index_bundle": "/groups/cgsd/alexandre/GATK_workflows/mutect2/inputs/Homo_sapiens_assembly38.index_bundle"  
    }
        )
    
    return result

if __name__=='__main__':
    output = Path('liver.json')
    with output.open('w') as f:
        json.dump(create_data(),f,indent=4)        