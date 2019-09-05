## rna-quantification.wdl
##
## This WDL runs the salmon command on a pair of FASTQ files
## (https://combine-lab.github.io/salmon/).
##
## Inputs:
## - Per-sample:
##   - Sample name
##   - R1 FASTQ file list
##   - R2 FASTQ file list
##
## - Reference:
##   - Gene map (GTF) file.
##   - Salmon index directory tar file (.tar.gz)
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - VM memory
##   - Runtime zones, ex: "us-central1-c us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - Transcript Quantification files
##   - Gene Quantification files
##   - Command info
##   - Parameter files
##   - Log files
##   - Auxilliary files

workflow RNAQuantification {
  String sample_name
  Array[String] fastq_1
  Array[String] fastq_2

  File gene_map
  # WDL does not support localizing a directory.
  # Workaround that by passing a TAR file.
  File salmon_index

  String salmon_docker
  Int salmon_vm_disk_size_gb
  String salmon_vm_memory
  String runtime_zones
  Int preemptible_tries

  call salmon_quant {
    input:
      sample_name = sample_name,
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,

      gene_map = gene_map,
      salmon_index = salmon_index,
      
      salmon_docker = salmon_docker,
      salmon_vm_disk_size_gb = salmon_vm_disk_size_gb,
      salmon_vm_memory = salmon_vm_memory,
      runtime_zones = runtime_zones,
      preemptible_tries = preemptible_tries
   }

  output {
    File quant_genes_sf = salmon_quant.quant_genes_sf
    File quant_sf = salmon_quant.quant_sf
    File cmd_info_json = salmon_quant.cmd_info_json
    File lib_format_counts_json = salmon_quant.lib_format_counts_json
    Array[File] logs = salmon_quant.logs
    Array[File] libParams = salmon_quant.libParams
    Array[File] aux_info = salmon_quant.aux_info
    Array[File] bootstraps = salmon_quant.bootstraps
  }
}

task salmon_quant {
  String sample_name
  Array[File] fastq_1
  Array[File] fastq_2

  File gene_map
  File salmon_index

  Int salmon_vm_disk_size_gb
  String salmon_vm_memory
  String salmon_docker
  String runtime_zones
  Int preemptible_tries

  command {
    set -o errexit
    set -o nounset
    set -o pipefail
    
    mkdir salmon_index
    tar xfz "${salmon_index}" --directory salmon_index
    gunzip "${gene_map}"
    UNZIPPED_GENE_MAP="$(echo ${gene_map} | sed -e "s/.gz$//")"

    salmon quant \
      --no-version-check \
      --index salmon_index \
      --libType A \
      --mates1 ${sep=' ' fastq_1} \
      --mates2 ${sep=' ' fastq_2} \
      --threads 16 \
      --numBootstraps 100 \
      --seqBias \
      --gcBias \
      --dumpEq \
      --geneMap "$UNZIPPED_GENE_MAP" \
      --output output
  }
  runtime {
    docker: salmon_docker

    disks: "local-disk " + salmon_vm_disk_size_gb + " HDD"

    memory: salmon_vm_memory
    cpu: 16
    zones: runtime_zones
    preemptible: preemptible_tries
  }
  output {
    File quant_genes_sf = "output/quant.genes.sf"
    File quant_sf = "output/quant.sf"
    File cmd_info_json = "output/cmd_info.json"
    File lib_format_counts_json = "output/lib_format_counts.json"

    # Contains salmon_quant.log
    Array[File] logs = glob("output/logs/*")

    # Contains flenDist.txt
    Array[File] libParams = glob("output/libParams/*")

    # Thirteen outputs produced:
    # ambig_info.tsv, eq_classes.txt, exp3_seq.gz, exp5_seq.gz,
    # expected_bias.gz, fld.gz, meta_info.json, obs3_seq.gz, obs5_seq.gz,
    # obs_gc.gz, observed_bias.gz, observed_bias_3p.gz
    Array[File] aux_info = glob("output/aux_info/*")

    # Two outputs produced:
    # bootstraps.gz and names.tsv.gz
    Array[File] bootstraps = glob("output/aux_info/bootstrap/*")
  }
}

