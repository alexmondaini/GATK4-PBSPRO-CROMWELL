version 1.0

# this workflow will create a <ref.fasta>.fai with task CreateFaidx and all the other indices files (fasta.alt,fasta.amb,fasta.ann,fasta.bwt,fasta.pac,fasta.sa) with task IndexBWA,
# after that it will create a sequence dictionary (.dict) from the fasta file with task CreateSequenceDictionary,
# and once the dictionary is created it will create an interval_list with dictionary output along with a supplied .bed file.

#Inputs:
# fasta file is compulsory
# bed file is optional in case you don't need an interval list

#Outputs:
# (fasta.alt,fasta.amb,fasta.ann,fasta.bwt,fasta.pac,fasta.sa) and (fasta.fai) from fasta
# (.dict and .interval_list) from bed

# Tools used:
# samtools v1.9-4
# bwa v0.7.17
# picard v2.23.8


workflow IndexGenome {

    input {
        File fasta
        File? bed_file
        File? probe_file
        Boolean create_fasta_index = false
        Boolean index_the_reference = false
        Boolean create_interval_list = false
    }

    if (create_fasta_index) {
    call CreateFaidx {
        input:
        fasta = fasta
    }
    }

    if (index_the_reference) {
        call IndexBWA {
            input:
            fasta = fasta
        }
    }

    call CreateSequenceDictionary {
        input:
        fasta = fasta
    }

    if (create_interval_list) {
        call BedToIntervalList {
            input:
            bed_file = bed_file,
            dict_file = CreateSequenceDictionary.output_dict
        }
    }
    if (create_interval_list) {
        call ProbeToIntervalList {
            input:
            probe_file = probe_file,
            dict_file = CreateSequenceDictionary.output_dict
        }
    }
}

task CreateFaidx {
    input {
        File fasta
    }
    command {
        samtools faidx ~{fasta}
    }
    runtime {
        docker: "biocontainers/samtools@sha256:da61624fda230e94867c9429ca1112e1e77c24e500b52dfc84eaf2f5820b4a2a"
        cpu: 2
        memory: '4 GB'
    }
    output {
        File output_fai = stdout()
    }
}

task IndexBWA {
    input {
        File fasta
    }
    command {
        bwa index -a bwtsw ~{fasta}
    }
    runtime {
        docker: "biocontainers/bwa@sha256:9479b73e108ded3c12cb88bb4e918a5bf720d7861d6d8cdbb46d78a972b6ff1b"
        cpu: 12
        memory: "50 GB"
    }
    output {
        File output_indexes = stdout()
    }
}

task CreateSequenceDictionary {
    input {
        File fasta
    }

    String strip_fasta_extension = sub(basename(fasta),"\\.f.*$","")

    command {
        java -Xmx4g -jar \
        /usr/picard/picard.jar CreateSequenceDictionary \
        R=~{fasta} \
        O=~{strip_fasta_extension}.dict
    }
    runtime {
        docker : "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:bb9207a31bdbdd96e34623ecde44d8a45e11bcbabacb34bbf860aa99d841cfea"
        cpu: 4
        memory: "4 GB"
    }
    output {
        File output_dict = "~{strip_fasta_extension}.dict"
    }
}

task BedToIntervalList {
    input {
        File? bed_file
        File dict_file
    }

    String output_file = basename(dict_file,'.dict') + ".interval_list"

    command {
        java -Xmx4g -jar \
        /usr/picard/picard.jar BedToIntervalList \
        I=~{bed_file} \
        O=~{output_file} \
        SD=~{dict_file}
    }
    runtime {
        docker : "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:bb9207a31bdbdd96e34623ecde44d8a45e11bcbabacb34bbf860aa99d841cfea"
        cpu: 4
        memory: "4 GB"
    }
    output {
        File output_interval_list = output_file
    }
}

task ProbeToIntervalList {
    input {
        File? probe_file
        File dict_file
    }

    String output_file = basename(dict_file,'.dict') + ".interval_list"

    command {
        java -Xmx4g -jar \
        /usr/picard/picard.jar BedToIntervalList \
        I=~{probe_file} \
        O=~{output_file} \
        SD=~{dict_file}
    }
    runtime {
        docker : "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:bb9207a31bdbdd96e34623ecde44d8a45e11bcbabacb34bbf860aa99d841cfea"
        cpu: 4
        memory: "4 GB"
    }
    output {
        File output_interval_list = output_file
    }
}