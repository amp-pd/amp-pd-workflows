# RNASummarization with subread featureCounts

This workflow runs the subread featureCounts command on a BAM file
(http://subread.sourceforge.net/).

Inputs:
- Per-sample:
  - Sample Name
  - BAM file

- Reference:
  - Gene map (GTF) file

- VM configuration
  - Docker image url
  - VM disk size
  - VM memory
  - Num CPU cores
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - &lt;sample_name&gt;.featureCounts.tsv
  - &lt;sample_name&gt;.featureCounts.tsv.summary

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).
A terra.inputs.json and a terra.outputs.json file that you can directly upload to Terra is also included here.

The gene map (GTF) file used was [GENCODE](https://www.gencodegenes.org/) v29.

The `subread` Docker image used was from [BioContainers](https://biocontainers.pro).
