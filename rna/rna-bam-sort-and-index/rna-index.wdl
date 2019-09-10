## rna-index.wdl
##
## This workflow runs samtools sort and samtools index on a bam file
## (http://samtools.sourceforge.net/).
##
## Inputs:
## - Per-sample:
##   - Sample Name
##   - BAM file
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - VM memory
##   - Runtime zones, ex: "us-central1-a us-central1-b"
##   - Num CPU cores
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - <sample_name>.samtools.bam
##   - <sample_name>.samtools.bam.bai
##   - <sample_name>.idxstats.txt
##   - <sample_name>.flagstat.txt

workflow RNAIndex {
  String sample_name
  File input_bam_file

  String samtools_docker
  Int vm_disk_size_gb
  String vm_memory
  String runtime_zones
  Int num_cpu_cores
  Int preemptible_tries

  call samtools_index {
    input:
      sample_name=sample_name,
      input_bam_file=input_bam_file,
      samtools_docker=samtools_docker,
      vm_disk_size_gb=vm_disk_size_gb,
      vm_memory=vm_memory,
      runtime_zones=runtime_zones,
      num_cpu_cores=num_cpu_cores,
      preemptible_tries=preemptible_tries
   }

  output {
    File samtools_bam_output_path = samtools_index.samtools_bam_output_path
    File samtools_bai_output_path = samtools_index.samtools_bai_output_path
    File samtools_idxstats_output_path = samtools_index.samtools_idxstats_output_path
    File samtools_flagstat_output_path = samtools_index.samtools_flagstat_output_path
  }
}

task samtools_index {

  String sample_name
  File input_bam_file

  String samtools_docker
  Int vm_disk_size_gb
  String vm_memory
  String runtime_zones
  Int num_cpu_cores
  Int preemptible_tries

  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail

    echo "$(date "+%Y-%m-%d %H:%M:%S") Creating soft link to old bam"
    ln -s "${input_bam_file}" old_bam_file.bam

    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting samtools sort"
    samtools sort -o "${sample_name}.samtools.bam" "${input_bam_file}"
    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting samtools index on new bam"
    samtools index -b "${sample_name}.samtools.bam"

    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting samtools idxstats on new bam"
    samtools idxstats "${sample_name}.samtools.bam" > idxstats.new.txt

    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting samtools flagstat on old bam"
    samtools flagstat old_bam_file.bam > flagstat.old.txt
    echo "$(date "+%Y-%m-%d %H:%M:%S") Starting samtools flagstat on new bam"
    samtools flagstat "${sample_name}.samtools.bam" > flagstat.new.txt

    # Verifying the old and new BAMs produce the same flagstats
    if ! diff flagstat.old.txt flagstat.new.txt; then
      1>&2 echo "BAM flagstats differ!"
      exit 1
    fi
    ln -s idxstats.new.txt "${sample_name}.idxstats.txt"
    ln -s flagstat.new.txt "${sample_name}.flagstat.txt"

    rm old_bam_file.bam
  >>>

  runtime {
    docker: samtools_docker

    disks: "local-disk " + vm_disk_size_gb + " HDD"

    memory: vm_memory
    cpu: num_cpu_cores
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File samtools_bam_output_path = "${sample_name}.samtools.bam"
    File samtools_bai_output_path = "${sample_name}.samtools.bam.bai"
    File samtools_idxstats_output_path = "${sample_name}.idxstats.txt"
    File samtools_flagstat_output_path = "${sample_name}.flagstat.txt"
  }
}

