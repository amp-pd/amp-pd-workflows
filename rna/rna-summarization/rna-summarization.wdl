## rna-summarization.wdl
##
## This workflow runs the subread featureCounts command on a BAM file
## (http://subread.sourceforge.net/).
##
## Inputs:
## - Per-sample:
##   - Sample Name
##   - BAM file
##
## - Reference:
##   - Gene map (GTF) file
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - VM memory
##   - Num CPU cores
##   - Runtime zones, ex: "us-central1-c us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - <sample_name>.featureCounts.tsv
##   - <sample_name>.featureCounts.tsv.summary

workflow RNASummarization {
  File bam_file
  String sample_name

  File gene_map

  String featureCounts_docker
  Int featureCounts_vm_disk_size_gb
  String featureCounts_vm_memory
  Int num_cpu_cores
  String runtime_zones
  Int preemptible_tries

  call subread_featureCounts {
    input:
      bam_file = bam_file,
      sample_name = sample_name,

      gene_map = gene_map,

      featureCounts_docker = featureCounts_docker,
      featureCounts_vm_disk_size_gb = featureCounts_vm_disk_size_gb,
      num_cpu_cores = num_cpu_cores,
      runtime_zones = runtime_zones,
      preemptible_tries = preemptible_tries,
      featureCounts_vm_memory = featureCounts_vm_memory
   }

  output {
    File output_tsv_file = subread_featureCounts.output_tsv_file
    File summary_file = subread_featureCounts.summary_file
  }
}

task subread_featureCounts {
  File bam_file
  String sample_name

  File gene_map

  String featureCounts_docker
  Int featureCounts_vm_disk_size_gb
  Int num_cpu_cores
  String runtime_zones
  Int preemptible_tries
  String featureCounts_vm_memory

  command {
    set -o errexit
    set -o nounset
    set -o pipefail

    # Set the output file name
    OUTPUT="$(pwd)/${sample_name}.featureCounts.tsv"

    # Unzip the gene map
    gunzip "${gene_map}"
    UNZIPPED_GENE_MAP="$(echo ${gene_map} | sed -e 's/.gz$//')"

    # Run featureCounts in the directory with the BAM; 
    # otherwise the output file records the full path in the header
    cd $(dirname "${bam_file}")

    featureCounts \
      -T 2 \
      -p \
      -t exon \
      -g gene_id \
      -a "$UNZIPPED_GENE_MAP" \
      -s 2 \
      -o "$OUTPUT" \
      "$(basename "${bam_file}")"
  }
  runtime {
    docker: featureCounts_docker

    disks: "local-disk " + featureCounts_vm_disk_size_gb + " HDD"
    zones: runtime_zones
    preemptible: preemptible_tries
    memory: featureCounts_vm_memory
    cpu: num_cpu_cores
  }
  output {
    # Two outputs produced:
    # - <sample_name}.featureCounts.tsv
    # - <sample_name>.featureCounts.tsv.summary
    File output_tsv_file = "${sample_name}.featureCounts.tsv"
    File summary_file = "${sample_name}.featureCounts.tsv.summary"
  }
}
