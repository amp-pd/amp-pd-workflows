# RNACollectMultipleMetrics with Picard

This workflow runs `picard CollectInsertSizeMetrics`, `picardCollectAlignmentSummaryMetrics`, `picard QualityScoreDistribution`, and `picard MeanQualityByCycle` on a BAM file
(https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_analysis_CollectInsertSizeMetrics.php).

Inputs:
- Per-sample:
  - Sample name - Name of the sample associated with the BAM - used in names of outputs.
  - BAM file - Path to a BAM file.
  - Ref FASTA file - Reference Sequence FASTA file.
- VM configuration
  - Docker image url
  - VM disk size
  - Num CPU cores
  - VM memory
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample:
  - File outputFiles - glob(*)
  - <sample-id>.Multiple_Metrics.<metric>

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The `picard` Docker image used was from [Broad Institute](https://hub.docker.com/r/broadinstitute/picard).