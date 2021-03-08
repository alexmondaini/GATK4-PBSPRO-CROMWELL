version 1.0

workflow GenomicsDB {
    input {
        File reference
        File ref_fai
        File intervals
        Array[File] intervals_idx
        File cohort
    }

    call DBImport {
        input:
        reference = reference,
        ref_fai = ref_fai,
        intervals = intervals,
        cohort = cohort,
        intervals_idx = intervals_idx
    }
}

task DBImport {
    input {
        File reference
        File ref_fai
        File intervals
        Array[File] intervals_idx
        File cohort
    }

    command <<<
        mkdir genome_db

        gatk --java-options "-Xms17G -Xmx17G" GenomicsDBImport \
        --genomicsdb-workspace-path genome_db \
        -R ~{reference} \
        -L ~{intervals} \
        --sample-name-map ~{cohort}
    >>>
    runtime {
        docker: "broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09"
        cpu: 6
        memory: "20 GB"
    }
    output {
        File GenomicsDB_workspace = stdout()
    }
}