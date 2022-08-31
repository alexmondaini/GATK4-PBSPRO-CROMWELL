from pathlib import Path
import json

bams = sorted(Path('/groups/cgsd/alexandre/liver/bams').glob('600-*.hg38.ba?'))
print(bams)

def pairwise(iterator):
    a = iter(iterator)
    return zip(a,a)

normals = list(pairwise([x for x in bams if x.name.split('.')[0].split('-')[1] == '1']))
tumors  = list(pairwise([x for x in bams if x.name.split('.')[0].split('-')[1] != '1']))


data = []
for n in normals:
    for t in tumors:
        if n[1].stem.split('.')[0].split('-')[0] == t[1].stem.split('.')[0].split('-')[0]:
            data.append((n[1],n[0],t[1],t[0]))


def create_data():
    result = []
    for file in data:
        result.append( {
            "CNVSomaticPairWorkflow.gatk_docker": "broadinstitute/gatk:4.2.0.0",
            "CNVSomaticPairWorkflow.common_sites": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/S07604715_Regions.interval_list",

            "CNVSomaticPairWorkflow.normal_bam": f"{file[0]}",
            "CNVSomaticPairWorkflow.normal_bam_idx": f"{file[1]}",
            "CNVSomaticPairWorkflow.tumor_bam": f"{file[2]}",
            "CNVSomaticPairWorkflow.tumor_bam_idx": f"{file[3]}",

            "CNVSomaticPairWorkflow.blacklist_intervals": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/CNV_and_centromere_blacklist.hg38liftover.list",
            "CNVSomaticPairWorkflow.padding": 100,
            "CNVSomaticPairWorkflow.bin_length": 0,

            "CNVSomaticPairWorkflow.ref_fasta": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta",
            "CNVSomaticPairWorkflow.ref_fasta_fai": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.fasta.fai",
            "CNVSomaticPairWorkflow.ref_fasta_dict": "/groups/cgsd/alexandre/GATK_workflows/src/hg38/Homo_sapiens_assembly38.dict",

            "CNVSomaticPairWorkflow.intervals": "/groups/cgsd/alexandre/GATK_workflows/copy-number-variation/inputs/S07604715_Regions.interval_list",
            "CNVSomaticPairWorkflow.read_count_pon": "/groups/cgsd/alexandre/cromwell-executions/CNVSomaticPanelWorkflow/LIVER/call-CreateReadCountPanelOfNormals/execution/Liver.pon.hdf5"
        }
        )
    return result

if __name__=='__main__':
    output = Path('LiverPair.json')
    with output.open('w') as f:
        json.dump(create_data(),f,indent=4) 