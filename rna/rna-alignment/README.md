# RNAAlignment with STAR

This workflow runs the STAR command on a pair of fastq files.
(https://github.com/alexdobin/STAR).

Inputs:
- Per-sample:
  - Sample name
  - R1 FASTQ file list
  - R2 FASTQ file list

- Reference:
  - STAR index directory tar file (.tar.gz)

- VM configuration
  - docker image url
  - VM disk size
  - VM memory
  - Runtime zones, ex: "us-central1-c us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - &lt;sample-id&gt;.Aligned.sortedByCoord.out.bam
  - &lt;sample-id&gt;.SJ.out.tab
  - &lt;sample-id&gt;.Chimeric.out.junction
  - &lt;sample-id&gt;.Log.final.out
  - &lt;sample-id&gt;.Log.out
  - &lt;sample-id&gt;.Log.progress.out
    - &lt;sample-id&gt;.STARgenome
      - sjdbInfo.txt
      - sjdbList.out.tab
    - &lt;sample-id&gt;.STARpass1
      - Log.final.out
      - SJ.out.tab

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra.](https://app.terra.bio/)

The `star` Docker image used was from [Alex Dobin](https://hub.docker.com/r/alexdobin/star/).
