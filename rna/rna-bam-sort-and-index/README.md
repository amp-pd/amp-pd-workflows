# RNABAMSortAndIndex with Samtools

This workflow runs samtools sort and samtools index on a bam file.
(http://samtools.sourceforge.net/).

Inputs:
- Per-sample:
  - Sample Name
  - BAM file
- VM configuration
  - Docker image url
  - VM disk size
  - VM memory
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Num CPU cores
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - &lt;sample_name&gt;.samtools.bam
  - &lt;sample_name&gt;.samtools.bam.bai
  - &lt;sample_name&gt;.idxstats.txt
  - &lt;sample_name&gt;.flagstat.txt

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The `samtools` Docker image used was from [BioContainers](https://biocontainers.pro).