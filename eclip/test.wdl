version 1.0 

workflow Test {
    call CreateDir {

    }
}

task CreateDir {
    command <<<
    mkdir outputdir
    touch outputdir/hello.fq
    touch outputdir/hello.bam
    >>>
    runtime {
        cpu: 1
        memory: "2 GB"
    }
    output {
        File result_fq = glob('outputdir/*fq')
        File result_bam = glob('outputdir/*.bam')
    }
}