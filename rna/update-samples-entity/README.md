# Update Terra `samples` entity from a samples list file

This workflow reads a "mapping file" from Cloud Storage
and writes Array["gs://..."] to the samples entity

Inputs:
- Mapping (tsv) file.
- Sample ID
- Runtime zones, ex: "us-central1-a us-central1-b"
- Number of times to try the workflow with a preemptible VM before
   falling back to a full-price VM.

Outputs:
- Per-sample
  - Array[String] of fastq files

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The mapping file is a tsv with three columns: sample_id, fastq_1, and fastq_2. Each row represents a pair of fastqs in GCS for a sample. A single sample may have multiple pairs of fastq files. As an example:
```
sample_id fastq_1 fastq_2
SAMPLE_ID_1 gs://bucket/SAMPLE_ID_1.L1.R1.fastq.gz gs://bucket/SAMPLE_ID_1.L1.R2.fastq.gz
SAMPLE_ID_1 gs://bucket/SAMPLE_ID_1.L2.R1.fastq.gz gs://bucket/SAMPLE_ID_1.L2.R2.fastq.gz
SAMPLE_ID_1 gs://bucket/SAMPLE_ID_1.L3.R1.fastq.gz gs://bucket/SAMPLE_ID_1.L3.R2.fastq.gz
SAMPLE_ID_2 gs://bucket/SAMPLE_ID_2.L1.R1.fastq.gz gs://bucket/SAMPLE_ID_2.L1.R2.fastq.gz
SAMPLE_ID_2 gs://bucket/SAMPLE_ID_2.L2.R1.fastq.gz gs://bucket/SAMPLE_ID_2.L2.R2.fastq.gz
SAMPLE_ID_3 gs://bucket/SAMPLE_ID_3.R1.fastq.gz gs://bucket/SAMPLE_ID_3.R2.fastq.gz
```

For the above example, the resulting `samples` table would then look like:
```
sample_id fastq_1 fastq_2
SAMPLE_ID_1 ["gs://bucket/SAMPLE_ID_1.L1.R1.fastq.gz", "gs://bucket/SAMPLE_ID_1.L2.R1.fastq.gz", "gs://bucket/SAMPLE_ID_1.L3.R1.fastq.gz"] ["gs://bucket/SAMPLE_ID_1.L1.R2.fastq.gz", "gs://bucket/SAMPLE_ID_1.L2.R2.fastq.gz", "gs://bucket/SAMPLE_ID_1.L3.R2.fastq.gz"]
SAMPLE_ID_2 ["gs://bucket/SAMPLE_ID_2.L1.R1.fastq.gz", "gs://bucket/SAMPLE_ID_2.L2.R1.fastq.gz"] ["gs://bucket/SAMPLE_ID_2.L1.R2.fastq.gz", "gs://bucket/SAMPLE_ID_2.L2.R2.fastq.gz"]
SAMPLE_ID_3 ["gs://bucket/SAMPLE_ID_3.R1.fastq.gz"] ["gs://bucket/SAMPLE_ID_3.R2.fastq.gz"]
