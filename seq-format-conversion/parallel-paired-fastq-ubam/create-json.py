from pprint import pprint
from pathlib import Path
import json

p = Path('/groups','cgsd','javed','liver')


def pairwise(iterable):
    obj = iter(iterable)
    return zip(obj,obj)

def altelement(a):
    return a[::2]


samples = sorted(p.glob('*/*fastq.gz'))
sample_names = altelement([x.name.split('_')[0] for x in samples])

print(len(sample_names))

data = []
for name,sample in zip(sample_names,pairwise(samples)):
    data.append({"name":f"{name}","readgroup":f"{name}","fastq_1":f"{sample[0]}","fastq_2":f"{sample[1]}","library":f"{name}"})
    

config = {
    "ConvertPairedFastQsToUnmappedBamWf.platform_name": "illumina",
    "ConvertPairedFastQsToUnmappedBamWf.samples": data,
    "ConvertPairedFastQsToUnmappedBamWf.sequencing_center": "CGS"
}



if __name__=='__main__':
    pass

    # output = Path('liver.json')
    # with output.open('w') as file:
    #     json.dump(config,file,indent=4)


