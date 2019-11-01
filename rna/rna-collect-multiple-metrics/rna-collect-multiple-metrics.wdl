## rna-collect-multiple-metrics.wdl
##
##
## This WDL runs Picard's CollectMultipleMetrics, executing the programs:
##   - CollectInsertSizeMetrics
##   - CollectAlignmentSummaryMetrics
##   - QualityScoreDistribution
##   - MeanQualityByCycle
##
## Inputs:
## - Per-sample:
##   - Sample name - Name of the sample associated with the BAM - used in names of outputs.
##   - BAM file - Path to a BAM file.
##   - Ref FASTA file - Reference Sequence FASTA file
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - Num CPU cores
##   - VM memory
##   - Runtime zones, ex: "us-central1-a us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample:
##   - <sample-id>.Multiple_Metrics.alignment_summary_metrics
##   - <sample-id>.Multiple_Metrics.base_distribution_by_cycle.pdf
##   - <sample-id>.Multiple_Metrics.base_distribution_by_cycle_metrics
##   - <sample-id>.Multiple_Metrics.insert_size_histogram.pdf
##   - <sample-id>.Multiple_Metrics.insert_size_metrics
##   - <sample-id>.Multiple_Metrics.quality_by_cycle.pdf
##   - <sample-id>.Multiple_Metrics.quality_by_cycle_metrics
##   - <sample-id>.Multiple_Metrics.quality_distribution.pdf
##   - <sample-id>.Multiple_Metrics.quality_distribution_metrics


workflow RNACollectMultipleMetrics {

  String sample_name
  File bam_file
  File reference_sequence_file

  String picard_docker

  Int vm_disk_size_gb
  Int num_cpu_cores
  String vm_memory

  String runtime_zones
  Int preemptible_tries
  
  call collect_multiple_metrics {
    input:
      sample_name=sample_name,
      bam_file=bam_file,
      reference_sequence_file=reference_sequence_file,
      
      picard_docker=picard_docker,

      vm_disk_size_gb=vm_disk_size_gb,
      num_cpu_cores=num_cpu_cores,
      vm_memory=vm_memory,

      runtime_zones=runtime_zones,
      preemptible_tries=preemptible_tries
   }

  output {
    Array[File] metrics_output_path = collect_multiple_metrics.metrics_output_path
  }
}

task collect_multiple_metrics {

  String sample_name
  File bam_file
  File reference_sequence_file
  
  String picard_docker

  Int vm_disk_size_gb
  Int num_cpu_cores
  String vm_memory

  String runtime_zones
  Int preemptible_tries
  
  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail

    java -jar /usr/picard/picard.jar CollectMultipleMetrics \
      INPUT=${bam_file} \
      OUTPUT="${sample_name}.Multiple_Metrics" \
      REFERENCE_SEQUENCE=${reference_sequence_file} \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT

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
    Array[File] metrics_output_path = glob("${sample_name}.Multiple_Metrics.*")
  }
}
