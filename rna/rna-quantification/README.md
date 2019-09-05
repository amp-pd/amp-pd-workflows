# RNAQuantification with salmon quant

This workflow runs the `salmon quant` command on a list of paired FASTQ files
(https://combine-lab.github.io/salmon/).

Inputs:
- Per-sample:
  - Sample name
  - R1 FASTQ file list
  - R2 FASTQ file list

- Reference:
  - Gene map (GTF) file.
  - Salmon index directory tar file (.tar.gz)

- VM configuration
  - Docker image url
  - VM disk size
  - VM memory
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - Transcript Quantification files
  - Gene Quantification files
  - Command info
  - Parameter files
  - Log files
  - Auxilliary files

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The gene map (GTF) file used was [GENCODE](https://www.gencodegenes.org/) v29.

The `salmon` Docker image used was from [combinelab](https://combine-lab.github.io/salmon/).