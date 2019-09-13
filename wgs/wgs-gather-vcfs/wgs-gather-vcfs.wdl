## wgs-gather-vcfs.wdl
## 
## This workflow runs the GATK's `GatherVcfsCloud` on a GCS directory containing sharded VCF files
## from a joint genotyping workflow such as the one here: https://github.com/gatk-workflows/gatk4-germline-snps-indels.
##
## Each VCF file corresponds to an interval from the joint genotyping intervals file. That file is a required
## input to this workflow.
##
## The output of the workflow is 1 VCF file per chromosome requested.
##
## Inputs:
## - intervals_file: GCS path of intervals file
## - chromosomes: An array of chromosome to gather, such as ["1", "2", "21", "22", "X", "Y"]
## - shards_to_ignore: An array of numbers from the intervals file to ignore, ex [576, 1032]
## - output_suffix: Output file will look like "chr22${output_suffix}.vcf"
## - input_directory: GCS path of directory containing sharded vcf files
## - input_file_prefix: Input VCF files look like "${input_file_prefix}.n.${input_file_suffix}"
## - input_file_suffix: Input VCF files look like "${input_file_prefix}.n.${input_file_suffix}"
##
## - VM configuration
##   - Docker image url
##   - VM disk size needed for chromosome 1
##   - VM memory
##   - Num CPU cores
##   - Runtime zones, ex: "us-central1-a us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per chromosome 
##   - <chromosome>.vcf.gz
##   - <chromosome>.vcf.idx
##   - <chromosome>.vcf.gz.tbi


workflow GatherVcfsCloud {
  File intervals_file
  Array[String] chromosomes
  Array[Int] shards_to_ignore
  String output_suffix
  String input_directory
  String input_file_prefix
  String input_file_suffix

  String docker_image
  Int disk_size_gb_chr_1
  String vm_memory
  Int num_cpu_cores
  String runtime_zones
  Int num_preemptible_retries

  scatter (chromosome in chromosomes) {
    call gather_vcfs_cloud {
      input:
        intervals_file = intervals_file,
        chromosome = chromosome,
        shards_to_ignore = shards_to_ignore,
        output_suffix = output_suffix,
        input_directory = input_directory,
        input_file_prefix = input_file_prefix,
        input_file_suffix = input_file_suffix,

        docker_image = docker_image,
        disk_size_gb_chr_1 = disk_size_gb_chr_1,
        vm_memory = vm_memory,
        num_cpu_cores = num_cpu_cores,
        runtime_zones = runtime_zones,
        num_preemptible_retries = num_preemptible_retries,
    }
  }

  output {
    Array[File] combined_vcf = gather_vcfs_cloud.combined_vcf
    Array[File] combined_vcf_idx = gather_vcfs_cloud.combined_vcf_idx
    Array[File] combined_vcf_tbi = gather_vcfs_cloud.combined_vcf_tbi
  }
}

task gather_vcfs_cloud {
  File intervals_file
  String chromosome
  Array[Int] shards_to_ignore
  String output_suffix
  String input_directory
  String input_file_prefix
  String input_file_suffix

  String docker_image
  Int disk_size_gb_chr_1
  String vm_memory
  Int num_cpu_cores
  String runtime_zones
  Int num_preemptible_retries

  # This map gives an estimation of how much disk space we need per chromosome
  # Source: http://www.cshlp.org/ghg5_all/section/dna.shtml
  Map[String, Float] filesize_map = {
    '1': 1.0,
    '2': 0.98,
    '3': 0.8,
    '4': 0.77,
    '5': 0.73,
    '6': 0.69,
    '7': 0.64,
    '8': 0.59,
    '9': 0.57,
    '10': 0.55,
    '11': 0.55,
    '12': 0.54,
    '13': 0.47,
    '14': 0.44,
    '15': 0.42,
    '16': 0.37,
    '17': 0.33,
    '18': 0.32,
    '19': 0.24,
    '20': 0.26,
    '21': 0.2,
    '22': 0.21,
    'X': 0.63,
    'Y': 0.24,
  }
  # Add an additional 20% for compression
  Int vm_disk_size = ceil(filesize_map[chromosome] * disk_size_gb_chr_1 * 1.2)

  command {
    set -o errexit
    set -o nounset
    set -o pipefail

    START=$(($(grep -n "chr${chromosome}:" "${intervals_file}"| cut -d':' -f1 | head -1) - 1))
    END=$(($(grep -n "chr${chromosome}:" "${intervals_file}"  | cut -d':' -f1 | tail -1) - 1))

    INPUT_PARAMETERS=""
    for ((i=START;i<=END;i++)); do
        for IGNORED_NUMBER in ${sep=' ' shards_to_ignore}; do
        if [[ $i == $IGNORED_NUMBER ]]; then
          # Continue to the outer for loop
          continue 2
        fi
      done
      INPUT_PARAMETERS+=" -I ${input_directory}/${input_file_prefix}.$i.${input_file_suffix}"
    done

    COMMAND="gatk GatherVcfsCloud$INPUT_PARAMETERS -O chr${chromosome}${output_suffix}.vcf"

    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting gatk GatherVcfsCloud"
    $COMMAND

    echo "$(date "+%Y-%m-%d %H:%M:%S") Compressing combined vcf"
    bgzip -c "chr${chromosome}${output_suffix}.vcf" > "chr${chromosome}${output_suffix}.vcf.gz"

    echo "$(date "+%Y-%m-%d %H:%M:%S") Creating tbi file"
    tabix -p vcf "chr${chromosome}${output_suffix}.vcf.gz"
  }
  runtime {
    docker: docker_image
    disks: "local-disk " + vm_disk_size + " HDD"
    zones: runtime_zones
    preemptible: num_preemptible_retries
    memory: vm_memory
    cpu: num_cpu_cores
  }
  output {
    File combined_vcf = "chr${chromosome}${output_suffix}.vcf.gz"
    File combined_vcf_idx = "chr${chromosome}${output_suffix}.vcf.idx"
    File combined_vcf_tbi = "chr${chromosome}${output_suffix}.vcf.gz.tbi"
  }
}
