## Copyright Broad Institute, 2017
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - One or more read groups, one per pair of FASTQ files 
##
## Outputs :
## - Set of unmapped BAMs, one per read group
##
## Cromwell version support 
## - Successfully tested on v24
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {
  Array[String] readgroup_list
  Map[String, Array[File]] fastq_pairs
  Map[String, Array[String]] metadata

  Int preemptible_tries

  # Convert multiple pairs of input fastqs in parallel
  scatter (readgroup in readgroup_list) {

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        fastq_1 = fastq_pairs[readgroup][0],
        fastq_2 = fastq_pairs[readgroup][1],
        readgroup_name = readgroup,
        sample_name = metadata[readgroup][0],
        library_name = metadata[readgroup][1],
        platform_unit = metadata[readgroup][2],
        run_date = metadata[readgroup][3],
        platform_name = metadata[readgroup][4],
        sequencing_center = metadata[readgroup][5],
        platform_model = metadata[readgroup][6],

        preemptible_tries = preemptible_tries
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
  }
}

# TASK DEFINITIONS

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  File fastq_1
  File fastq_2
  String readgroup_name
  String sample_name
  String library_name
  String platform_unit
  String run_date
  String platform_name
  String sequencing_center
  String platform_model
  Int disk_size
  String mem_size
  String docker
  String gatk_path
  Int preemptible_tries

  command {
    ${gatk_path} --java-options "-Xmx3000m" \
      FastqToSam \
      --FASTQ ${fastq_1} \
      --FASTQ2 ${fastq_2} \
      --OUTPUT ${readgroup_name}.unmapped.bam \
      --READ_GROUP_NAME ${readgroup_name} \
      --SAMPLE_NAME ${sample_name} \
      --LIBRARY_NAME ${library_name} \
      --PLATFORM_UNIT ${platform_unit} \
      --RUN_DATE ${run_date} \
      --PLATFORM ${platform_name} \
      --SEQUENCING_CENTER ${sequencing_center} \
      --PLATFORM_MODEL ${platform_model} 
  }
  runtime {
    docker: docker
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}

