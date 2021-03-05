version 1.0

workflow GenomicsDB {
    input {
        File cohort
        File reference
        File ref_fai
        File intervals
    }

    call DBImport {
        input:
        cohort = cohort,
        reference = reference,
        ref_fai = ref_fai,
        intervals = intervals
    }
}

task DBImport {
    input {
        File cohort
        File reference
        File ref_fai
        File intervals
    }

    String genomics_db_dir = genomics_db_dir

    command {
        gatk --java-options "-Xmx 17G -Xmx17G" GenomicsDBImport \
        --genomicsdb-workspace-path ${genomics_db_dir} \
        -R ${reference} \
        -L ${intervals} \
        --sample-name-map ${cohort}
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09"
        cpu: 4
        memory: "20 GB"
    }
    output {
        File GenomicsDB_workspace = "${genomics_db_dir}"
    }
}