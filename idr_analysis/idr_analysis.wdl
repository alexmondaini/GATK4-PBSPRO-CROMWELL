version 1.0

workflow IDR {
    input {
        Array[File] ep_293T_siNeg
        Array[File] ep_293T_siFus
    }
    Array[Pair[File,File]] crossed = cross(ep_293T_siNeg,ep_293T_siFus)

    scatter (sample in crossed) {
        call Idr_Analysis {
            input:
            sample_left = sample.left,
            sample_right = sample.right
        }
    }
}

task Idr_Analysis {
    input {
        File sample_left
        File sample_right
    }
    String out_file = basename(sample_left,'.round2') + '_vs_' + basename(sample_right,'.round2')
    
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate idr
    idr --samples ~{sample_left} ~{sample_right} \
    --input-file-type bed \
    --rank p.value \
    --output-file ~{out_file} \
    --output-file-type bed \
    --verbose
    >>>

    runtime {
        cpu: 6
        memory: "10 GB"
    }

    output {
        File output_bed_file = "~{out_file}"
    }
}