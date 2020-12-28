version 1.0

workflow ValidateBamsWf {
  input {
    Array[File] bam_array 
    String gatk_docker = "broadinstitute/gatk@broadinstitute/gatk@sha256:33574f446ac991f77bac125fbf6a2340e6db972a3f334e6c61bff94740165938"
    String gatk_path = "/gatk/gatk"
    String gatk_path = "/gatk/gatk"
  }

  # Process the input files in parallel

  # Process the input files in parallel
  scatter (input_bam in bam_array) {

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(input_bam, ".bam")

    # Run the validation 
    call ValidateBAM {
      input:
        input_bam = input_bam,
        output_basename = bam_basename + ".validation",
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }
  # Outputs that will be retained when execution is complete
  output {
    Array[File] validation_reports = ValidateBAM.validation_report
  }

}

# TASK DEFINITIONS

# Validate a SAM or BAM using Picard ValidateSamFile
task ValidateBAM {
  input {
    # Command parameters
    File input_bam
    String output_basename
    String? validation_mode
    String gatk_path
  
    # Runtime parameters
    String docker
    Int machine_mem_gb = 16
    Int cores = 6
  }
    
  String output_name = "${output_basename}_${validation_mode}.txt"
 
  command {
    ${gatk_path} \
      ValidateSamFile \
      --INPUT ${input_bam} \
      --OUTPUT ${output_name} \
      --MODE ${default="SUMMARY" validation_mode}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: cores
  }
  output {
    File validation_report = "${output_name}"
    Boolean done = true
  }
}