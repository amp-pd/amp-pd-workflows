## update-samples-entity-from-samples-list.wdl
##
## This workflow reads a "mapping file" from Cloud Storage 
## and writes Array["gs://..."] to the samples entity
##
## Inputs:
## - Mapping (tsv) file.
## - Sample ID
## - Runtime zones, ex: "us-central1-a us-central1-b"
## - Number of times to try the workflow with a preemptible VM before
##   falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - Array[String] of fastq files

workflow SampleToPairedFastq {
  File mapping_file
  String sample_id
  String runtime_zones
  Int preemptible_tries

  call sample_to_paired_fastq_task_tsv {
    input:
      mapping_file = mapping_file,
      sample_id = sample_id,
      runtime_zones = runtime_zones,
      preemptible_tries = preemptible_tries
   }

  output {
    Array[String] fastq1 = sample_to_paired_fastq_task_tsv.fastq1
    Array[String] fastq2 = sample_to_paired_fastq_task_tsv.fastq2
  }
}

task sample_to_paired_fastq_task_tsv {
  File mapping_file
  String sample_id
  String runtime_zones
  Int preemptible_tries
  command <<<
    echo "${mapping_file}"
    python <<CODE
    import csv

    d = {}

    with open("${mapping_file}", 'r') as tsvfile:
      reader = csv.DictReader(tsvfile, delimiter='\t')
      for row in reader:
        sample_id = row['sample_id']
        fastq_1 = row['fastq_1']
        fastq_2 = row['fastq_2']
        if sample_id not in d:
          d[sample_id] = {'fastq_1': [fastq_1], 'fastq_2': [fastq_2]}
        else:
          d[sample_id]['fastq_1'].append(fastq_1)
          d[sample_id]['fastq_2'].append(fastq_2)
    sample_id = "${sample_id}"
    with open('fastq1_list.txt', 'w') as outputfile:
      for fastq in d[sample_id]['fastq_1']:
        outputfile.write(fastq + '\n')
    with open('fastq2_list.txt', 'w') as outputfile:
      for fastq in d[sample_id]['fastq_2']:
        outputfile.write(fastq + '\n')
    CODE
  >>>
  runtime {
    docker: "python:2.7-slim"
    zones: runtime_zones
    preemptible: preemptible_tries
  }
  output {
    Array[String] fastq1 = read_lines('fastq1_list.txt')
    Array[String] fastq2 = read_lines('fastq2_list.txt')
  }
}

task sample_to_paired_fastq_task_json {
  File mapping_file
  String sample_id
  String runtime_zones

  command <<<
    echo "${mapping_file}"
    python <<CODE
    import json
    with open("${mapping_file}", 'r') as inputfile:
      j = json.load(inputfile)
    sample_id = "${sample_id}"
    with open('fastq1_list.txt', 'w') as outputfile:
      for fastq in j[sample_id]['fastq1']:
        outputfile.write(fastq + '\n')
    with open('fastq2_list.txt', 'w') as outputfile:
      for fastq in j[sample_id]['fastq2']:
        outputfile.write(fastq + '\n')
    CODE
  >>>
  runtime {
    docker: "python:2.7-slim"
    zones: runtime_zones
    preemptible: 1
  }
  output {
    Array[String] fastq1 = read_lines('fastq1_list.txt')
    Array[String] fastq2 = read_lines('fastq2_list.txt')
  }
}
