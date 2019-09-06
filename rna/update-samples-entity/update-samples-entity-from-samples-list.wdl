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
    echo "Starting sample_to_paired_fastq_task_tsv on ${mapping_file}."
    python <<CODE

    import csv

    fastqs_by_sample_id = {}

    # The mapping file is a tsv with three columns: sample_id, fastq_1, and fastq_2.
    # Each row represents a pair of fastqs in GCS for a sample.
    # A single sample may have multiple pairs of fastq files.
    # Example:
    #   sample_id fastq_1 fastq_2
    #   SAMPLE_ID_1 gs://bucket/SAMPLE_ID_1.L1.R1.fastq.gz gs://bucket/SAMPLE_ID_1.L1.R2.fastq.gz
    #   SAMPLE_ID_1 gs://bucket/SAMPLE_ID_1.L2.R1.fastq.gz gs://bucket/SAMPLE_ID_1.L2.R2.fastq.gz
    with open("${mapping_file}", 'r') as tsvfile:

        # Turns each row in the TSV into a dictionary keyed by column name.
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for row in reader:
            sample_id = row['sample_id']
            fastq_1 = row['fastq_1']
            fastq_2 = row['fastq_2']

            # Add newly seen samples to the fastqs_by_sample_id dict...
            if sample_id not in fastqs_by_sample_id:
                fastqs_by_sample_id[sample_id] = {
                    'fastq_1': [fastq_1],
                    'fastq_2': [fastq_2]
                }
            # ...or append additional fastq files to previously seen samples.
            else:
                fastqs_by_sample_id[sample_id]['fastq_1'].append(fastq_1)
                fastqs_by_sample_id[sample_id]['fastq_2'].append(fastq_2)

    sample_id = "${sample_id}"

    # Write out the list of fastqs, which will be output to the sample entity.

    with open('fastq1_list.txt', 'w') as outputfile:
        for fastq in fastqs_by_sample_id[sample_id]['fastq_1']:
            outputfile.write(fastq + '\n')

    with open('fastq2_list.txt', 'w') as outputfile:
        for fastq in fastqs_by_sample_id[sample_id]['fastq_2']:
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
