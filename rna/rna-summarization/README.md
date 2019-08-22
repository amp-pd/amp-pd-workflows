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
  - Runtime zones, ex: "us-central1-c us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - &lt;sample_name&gt;.featureCounts.tsv
  - &lt;sample_name&gt;.featureCounts.tsv.summary

