## rna-alignment.wdl
##
## This workflow runs the STAR command on a pair of fastq files.
## (https://github.com/alexdobin/STAR)
##
## Inputs:
## - Per-sample:
##   - Sample name
##   - R1 FASTQ file list
##   - R2 FASTQ file list
##
## - Reference:
##   - STAR index directory tar file (.tar.gz)
##
## - VM configuration
##   - docker image url
##   - VM disk size
##   - VM memory
##   - Runtime zones, ex: "us-central1-c us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - <sample-id>.Aligned.out.bam
##   - <sample-id>.SJ.out.tab
##   - <sample-id>.Chimeric.out.junction
##   - <sample-id>.Log.final.out
##   - <sample-id>.Log.out
##   - <sample-id>.Log.progress.out
##   - <sample-id>._STARgenome
##     - sjdbInfo.txt
##     - sjdbList.out.tab
##   - <sample-id>._STARpass1
##     - Log.final.out
##     - SJ.out.tab

workflow RNAAlignment {
  String sample_name
  Array[String] fastq_1
  Array[String] fastq_2

  # WDL does not support localizing a directory.
  # Workaround that by passing a TAR file.
  File star_index

  String star_docker
  Int star_vm_disk_size_gb
  String star_vm_memory
  String runtime_zones
  Int preemptible_tries
  Float star_timeout_hours

  call star_align {
    input:
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      sample_name = sample_name,
      star_index = star_index,
      star_vm_disk_size_gb = star_vm_disk_size_gb,
      star_vm_memory = star_vm_memory,
      star_docker = star_docker,
      runtime_zones = runtime_zones,
      preemptible_tries = preemptible_tries,
      star_timeout_hours = star_timeout_hours,
   }

  output {
    Array[File] outputFiles = star_align.outputFiles
    Array[File] outputStarGenomeFiles = star_align.outputStarGenomeFiles
    Array[File] outputStarpass1Files = star_align.outputStarpass1Files
  }
}

task star_align {

  Array[File] fastq_1
  Array[File] fastq_2
  String sample_name
  File star_index

  Int star_vm_disk_size_gb
  String star_vm_memory
  String star_docker
  String runtime_zones
  Int preemptible_tries
  Float star_timeout_hours

  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail
    set -o xtrace

    # For each FASTQ1 file, extract fields to build up a list
    # of readgroup tags. The readgroup tag takes the form:
    # <FLOWCELL>.<INDEX>.<LANE>
    RGTAGLIST=""
    for FILE in ${sep=" " fastq_1}; do
      SAMPLE_NAME="$(basename $FILE | cut -d. -f1)"

      LINE="$(head -n 1 < <(zcat $FILE))"
      FLOWCELL="$(echo $LINE | cut -d: -f3)"
      LANE="$(echo $LINE | cut -d: -f4)"
      INDEX="$(echo $LINE | cut -d: -f10)"

      if [[ -z $LINE ]]; then
        2>&1 "ERROR: Unable to parse first line from $FILE"
        exit 1
      fi
      if [[ -z $FLOWCELL ]]; then
        2>&1 "ERROR: Unable to parse flowcell from $LINE"
        exit 1
      fi
      if [[ -z $LANE ]]; then
        2>&1 "ERROR: Unable to parse lane from $LINE"
        exit 1
      fi
      if [[ -z $INDEX ]]; then
        2>&1 "ERROR: Unable to parse index from $LINE"
        exit 1
      fi

      RGTAG="ID:$FLOWCELL.$INDEX.$LANE PL:ILLUMINA SM:$SAMPLE_NAME "
      RGTAGLIST="$RGTAGLIST , $RGTAG"
    done
    RGTAGLIST="$(echo $RGTAGLIST | sed 's/^, //')"

    # Untar the STAR index
    mkdir star_index
    tar xfz "${star_index}" --directory star_index --strip-components=1

    mkdir output_dir

    # Run STAR --runMode alignReads
    # We run using timeout to prevent long-running workflows
    # from wasteful preemption loops.
    timeout "${star_timeout_hours}"h STAR \
      --genomeDir star_index \
      --runMode alignReads \
      --twopassMode Basic \
      --outFileNamePrefix output_dir/${sample_name}. \
      --readFilesCommand zcat \
      --readFilesIn ${sep=',' fastq_1} ${sep=',' fastq_2} \
      --outSAMtype BAM Unsorted \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.1 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --chimOutType Junctions \
      --chimSegmentMin 15 \
      --chimJunctionOverhangMin 15 \
      --runThreadN 16 \
      --outSAMstrandField intronMotif \
      --outSAMunmapped Within \
      --outSAMattrRGline $RGTAGLIST
  >>>

  runtime {
    docker: star_docker

    disks: "local-disk " + star_vm_disk_size_gb + " HDD"

    memory: star_vm_memory
    cpu: 16
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    # 6 outputs produced:
    #   - <sample-id>.Aligned.out.bam
    #   - <sample-id>.SJ.out.tab
    #   - <sample-id>.Chimeric.out.junction
    #   - <sample-id>.Log.final.out
    #   - <sample-id>.Log.out
    #   - <sample-id>.Log.progress.out
    Array[File] outputFiles = glob("output_dir/*")

    # 2 outputs produced:
    #   - sjdbInfo.txt
    #   - sjdbList.out.tab
    Array[File] outputStarGenomeFiles = glob("output_dir/${sample_name}._STARgenome/*")

    # 2 outputs produced:
    #     - Log.final.out
    #     - SJ.out.tab
    Array[File] outputStarpass1Files = glob("output_dir/${sample_name}._STARpass1/*")

  }
}
