version 1.0

workflow AnalyzeCovariates {
    input {
        File table
        String gatk_path = "/gatk/gatk"
    }
    call Producepdf {
        input:
        table = table,
        gatk_path = gatk_path
    }
}

task Producepdf {
    input {
        File table
        String gatk_path
    }

    String out_pdf = basename(table,'.csv') + '.pdf'

    command {
        ~{gatk_path} --java-options "-Xmx3g" \
        AnalyzeCovariates \
        -bqsr ~{table} \
        -plots ~{out_pdf}
    }
    runtime {
        docker: "broadinstitute/gatk@sha256:33574f446ac991f77bac125fbf6a2340e6db972a3f334e6c61bff94740165938"
        cpu: 4
        memory: "4 GB"
    }
    output {
        File output_pdf_comparison = out_pdf
    }
}