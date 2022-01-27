workflow RNAseq {
    # Bams from STAR
	File inputBam
	String sampleName = basename(inputBam,".bam")

    # References
	File refFasta
	File refFastaIndex
	File refDict

    # Docker config
	String? gatk4_docker_override
	String gatk4_docker = select_first([gatk4_docker_override, "broadinstitute/gatk:latest"])
	String? gatk_path_override
	String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

    # Known indel sites in hg19
	Array[File] knownVcfs
	Array[File] knownVcfsIndices

    # SNP database
	File dbSnpVcf
	File dbSnpVcfIndex

    # The minimum phred-scaled confidence threshold at which variants should be called
	Int? minConfidenceForVariantCalling

	# Inputs for STAR
	File fusgenes
  
	## Optional user optimizations
	Int? haplotypeScatterCount
	Int scatterCount = select_first([haplotypeScatterCount, 6])

	call bedToCallingIntervals {
	    input:
	        bed = fusgenes,
	        ref_dict = refDict,
	        gatk_path = gatk_path,
	        docker = gatk4_docker
	}

	call MarkDuplicates {
		input:
			input_bam = inputBam,
			base_name = sampleName + ".dedupped",
			
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

    call SplitNCigarReads {
        input:
            input_bam = MarkDuplicates.output_bam,
            input_bam_index = MarkDuplicates.output_bam_index,
            base_name = sampleName + ".split",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


	call BaseRecalibrator {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

	call ApplyBQSR {
		input:
			input_bam =  SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}


    call ScatterIntervalList {
        input:
            interval_list = bedToCallingIntervals.interval_list,
            scatter_count = scatterCount,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


	scatter (interval in ScatterIntervalList.out) {
        call HaplotypeCaller {
            input:
                input_bam = ApplyBQSR.output_bam,
                input_bam_index = ApplyBQSR.output_bam_index,
                base_name = sampleName + ".hc",
                interval_list = interval,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                ref_dict = refDict,
                dbSNP_vcf = dbSnpVcf,
                dbSNP_vcf_index = dbSnpVcfIndex,
                stand_call_conf = minConfidenceForVariantCalling,
                docker = gatk4_docker,
                gatk_path = gatk_path
        }

		File HaplotypeCallerOutputVcf = HaplotypeCaller.output_vcf
		File HaplotypeCallerOutputVcfIndex = HaplotypeCaller.output_vcf_index
	}

    call MergeVCFs {
        input:
            input_vcfs = HaplotypeCallerOutputVcf,
            input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
            output_vcf_name = sampleName + ".g.vcf.gz",
    
            docker = gatk4_docker,
            gatk_path = gatk_path
    }
	
	call VariantFiltration {
		input:
			input_vcf = MergeVCFs.output_vcf,
			input_vcf_index = MergeVCFs.output_vcf_index,
			base_name = sampleName + ".variant_filtered.vcf.gz",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_vcf = MergeVCFs.output_vcf
		File merged_vcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_vcf = VariantFiltration.output_vcf
		File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
	}
}

task bedToCallingIntervals {
    File bed
    File ref_dict

    String output_name = basename(bed, ".txt") + ".exons.interval_list"

    String docker
    String gatk_path

    command <<<

        set -e

        Rscript --no-save -<<'RCODE'
            bed = read.table("${bed}", sep="\t", skip=1)
            write.table(data.frame(chrom=bed[,'V2'], start=bed[,'V3'], end=bed[,'V4']), "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
        RCODE

        awk '{print $1 "\t" ($2 - 1) "\t" $3}' exome.bed > exome.fixed.bed

        ${gatk_path} \
            BedToIntervalList \
            -I exome.fixed.bed \
            -O ${output_name} \
            -SD ${ref_dict}
    >>>

    output {
        File interval_list = "${output_name}"
    }

    runtime {
        docker: docker
    }
}

task MarkDuplicates {
    File input_bam
 	String base_name
 
    String gatk_path
    String docker
 
 	command <<<
     ${gatk_path} --java-options "-Xmx20g" \
     MarkDuplicates \
     --INPUT ${input_bam} \
     --OUTPUT ${base_name}.bam  \
     --CREATE_INDEX true \
     --VALIDATION_STRINGENCY SILENT \
     --METRICS_FILE ${base_name}.metrics
     >>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 		File metrics_file = "${base_name}.metrics"
 	}

	runtime {
		docker: docker
		memory: "24 GB"
	}
}

task SplitNCigarReads {
    File input_bam
    File input_bam_index
    String base_name

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String gatk_path
    String docker

    command <<<
    ${gatk_path} --java-options "-Xmx4g" \
    SplitNCigarReads \
    -R ${ref_fasta} \
    -I ${input_bam} \
    -O ${base_name}.bam 
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        docker: docker
        memory: "5 GB"
    }
}

task BaseRecalibrator {
    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path
    String docker

    command <<<
    ${gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
    -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
    -Xloggc:gc_log.log -Xms3g -Xms4g" \
    BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${input_bam} \
    --use-original-qualities \
    -O ${recal_output_file} \
    -known-sites ${dbSNP_vcf} \
    -known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
    >>>

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
        memory: "6 GB"
        docker: docker
    }
}

task ApplyBQSR {
    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path
    String docker

    command <<<
    ${gatk_path} \
    --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
    -XX:+PrintGCDetails -Xloggc:gc_log.log \
    -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3g -Xmx4g" \
    ApplyBQSR \
    --add-output-sam-program-record \
    -R ${ref_fasta} \
    -I ${input_bam} \
    --use-original-qualities \
    -O ${base_name}.bam \
    --bqsr-recal-file ${recalibration_report}
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        memory: "5 GB"
        docker: docker
    }
}

task ScatterIntervalList {
	File interval_list
	Int scatter_count
	String gatk_path
	String docker

    command <<<
        set -e
        mkdir out
        ${gatk_path} --java-options "-Xmx1g" \
            IntervalListTools \
            --SCATTER_COUNT ${scatter_count} \
            --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            --UNIQUE true \
            --SORT true \
            --INPUT ${interval_list} \
            --OUTPUT out
	
        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        if len(intervals) == 0:
          raise ValueError("Interval list produced 0 scattered interval lists. Is the gtf or input interval list empty?")
        f = open("interval_count.txt", "w+")
        f.write(str(len(intervals)))
        f.close()
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

    runtime {
        memory: "2 GB"
        docker: docker
    }
}

task HaplotypeCaller {
	File input_bam
	File input_bam_index
	String base_name

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	String gatk_path
	String docker

	Int? stand_call_conf

	command <<<
		${gatk_path} --java-options "-Xms5g -Xmx6g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O ${base_name}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
		--dbsnp ${dbSNP_vcf}
	>>>

	output {
		File output_vcf = "${base_name}.vcf.gz"
		File output_vcf_index = "${base_name}.vcf.gz.tbi"
	}

	runtime {
		docker: docker
		memory: "6.5 GB"
	}
}

task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    String gatk_path

    String docker

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        ${gatk_path} --java-options "-Xms2g -Xmx2g"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT ' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    runtime {
        memory: "3 GB"
        docker: docker
    }
}

task VariantFiltration {
	File input_vcf
	File input_vcf_index
	String base_name

 	File ref_dict
 	File ref_fasta
 	File ref_fasta_index

	String gatk_path
	String docker

	command <<<
    ${gatk_path} --java-options "-Xmx4g" \
    VariantFiltration \
    --R ${ref_fasta} \
    --V ${input_vcf} \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0" \
    -O ${base_name}
	>>>

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
		docker: docker
		memory: "6 GB"
	}
}