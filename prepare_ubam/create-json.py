from pprint import pprint
from pathlib import Path
import json
from io import StringIO 

p = Path('/groups','cgsd','alexandre','cromwell-executions','ConvertPairedFastQsToUnmappedBamWf','d8e57ecb-a44d-4d7e-bc17-2a8acaacb7df','call-PairedFastQsToUnmappedBAM')

samples = sorted(p.glob('shard-*/execution/*.unmapped.bam'))
samples_names = [x.name.split('.')[0] for x in samples]

def create_data():
    data = []
    for name,path in zip(samples_names,samples):
        data.append(
        {
        "PreProcessingForVariantDiscovery_GATK4.sample_name": f"{name}",
        "PreProcessingForVariantDiscovery_GATK4.ref_name": "hg38",
        "PreProcessingForVariantDiscovery_GATK4.unmapped_bam": f"{path}",
        "PreProcessingForVariantDiscovery_GATK4.unmapped_bam_suffix": ".unmapped.bam",
        "PreProcessingForVariantDiscovery_GATK4.ref_dict": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.dict",
        "PreProcessingForVariantDiscovery_GATK4.ref_fasta": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta",
        "PreProcessingForVariantDiscovery_GATK4.ref_fasta_index": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.fai",
        "PreProcessingForVariantDiscovery_GATK4.ref_alt": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.alt",
        "PreProcessingForVariantDiscovery_GATK4.ref_sa": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.sa",
        "PreProcessingForVariantDiscovery_GATK4.ref_amb": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.amb",
        "PreProcessingForVariantDiscovery_GATK4.ref_bwt": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.bwt",
        "PreProcessingForVariantDiscovery_GATK4.ref_ann": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.ann",
        "PreProcessingForVariantDiscovery_GATK4.ref_pac": "/groups/cgsd/alexandre/gatk-workflows/src/ref_without_HPV/Homo_sapiens_assembly38.fasta.64.pac",
        "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf": "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Homo_sapiens_assembly38.dbsnp138.vcf",
        "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index": "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs": [
            "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Homo_sapiens_assembly38.known_indels.vcf.gz"
        ],
        "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices": [
            "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
            "/groups/cgsd/alexandre/gatk-workflows/gatk4-exome-analysis-pipeline/inputs/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
        ]
        }
        )
    return data


if __name__=='__main__':
    output = Path('liver.json')
    with output.open('w') as f:
        json.dump(create_data(),f,indent=4)
