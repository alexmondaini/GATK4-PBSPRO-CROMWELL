# this workflow wil create a <ref.fasta>.fai
# Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence.
# link for tools search http://www.htslib.org/doc/samtools-faidx.html
version 1.0

workflow IndexGenome {

    input {
        File fasta
        File? bed_file
        Boolean index_the_reference = false
        Boolean create_interval_list = false
    }

    call CreateFaidx {
        input:
        fasta = fasta
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

    String bed_basename = basename(bed_file,'.bed')

    command {
        java -Xmx4g -jar \
        /usr/picard/picard.jar BedToIntervalList \
        I=~{bed_file} \
        O=~{bed_basename}.interval_list \
        SD=~{dict_file}
    }
    runtime {
        docker : "us.gcr.io/broad-gotc-prod/picard-cloud@sha256:bb9207a31bdbdd96e34623ecde44d8a45e11bcbabacb34bbf860aa99d841cfea"
        cpu: 4
        memory: "4 GB"
    }
    output {
        File output_interval_list = "~{bed_basename}.interval_list"
    }
}