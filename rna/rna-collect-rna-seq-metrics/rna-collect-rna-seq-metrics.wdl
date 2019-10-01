## rna-collect-rna-seq-metrics.wdl
##
## This workflow runs Picard's CollectRnaSeqMetrics on a BAM file
## (https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectRnaSeqMetrics.php).
##
## Inputs:
## - Per-sample:
##   - Sample name - Name of the sample associated with the BAM - used in names of outputs.
##   - BAM file - Path to a BAM file.
##   - Ref flat file - Gene annotations in refFlat form.
##   - Ribosomal intervals file - Location of rRNA sequences in genome, in interval_list format.
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - VM memory
##   - Num CPU cores
##   - Runtime zones, ex: "us-central1-a us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample:
##   - <sample-id>.RNA_Metrics

workflow RNACollectRnaSeqMetrics {

  String sample_name
  File bam_file
  File ref_flat_file
  File ribosomal_intervals_file

  String picard_docker

  Int vm_disk_size_gb
  String vm_memory
  Int num_cpu_cores
  String runtime_zones
  Int preemptible_tries
  
  call collect_rna_seq_metrics {
    input:
      sample_name=sample_name,
      bam_file=bam_file,
      ref_flat_file=ref_flat_file,
  	  ribosomal_intervals_file=ribosomal_intervals_file,
      
      picard_docker=picard_docker,

      vm_disk_size_gb=vm_disk_size_gb,
      vm_memory=vm_memory,
      num_cpu_cores=num_cpu_cores,
      runtime_zones=runtime_zones,
      preemptible_tries=preemptible_tries
   }

  output {
    File metrics_output_path = collect_rna_seq_metrics.metrics_output_path
  }
}

task collect_rna_seq_metrics {

  String sample_name
  File bam_file
  File ref_flat_file
  File ribosomal_intervals_file
  
  String picard_docker

  Int vm_disk_size_gb
  String vm_memory
  Int num_cpu_cores
  String runtime_zones
  Int preemptible_tries
  
  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail

    java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      INPUT="${bam_file}" \
      OUTPUT="${sample_name}.RNA_Metrics" \
      REF_FLAT="${ref_flat_file}" \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      RIBOSOMAL_INTERVALS="${ribosomal_intervals_file}"
  >>>

  runtime {
    docker: picard_docker

    disks: "local-disk " + vm_disk_size_gb + " HDD"

    memory: vm_memory
    cpu: num_cpu_cores
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File metrics_output_path = "${sample_name}.RNA_Metrics"
  }
}