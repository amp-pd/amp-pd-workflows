# RNACollectMultipleMetrics with Picard

This workflow runs the following `picard` commands on a BAM file:

- [CollectInsertSizeMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_analysis_CollectInsertSizeMetrics.php)
- [CollectAlignmentSummaryMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_analysis_CollectAlignmentSummaryMetrics.php)
- [QualityScoreDistribution](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_analysis_QualityScoreDistribution.php)
- [MeanQualityByCycle](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_analysis_MeanQualityByCycle.php)

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
  - <sample-id>.Multiple_Metrics.alignment_summary_metrics
  - <sample-id>.Multiple_Metrics.base_distribution_by_cycle.pdf
  - <sample-id>.Multiple_Metrics.base_distribution_by_cycle_metrics
  - <sample-id>.Multiple_Metrics.insert_size_histogram.pdf
  - <sample-id>.Multiple_Metrics.insert_size_metrics
  - <sample-id>.Multiple_Metrics.quality_by_cycle.pdf
  - <sample-id>.Multiple_Metrics.quality_by_cycle_metrics
  - <sample-id>.Multiple_Metrics.quality_distribution.pdf
  - <sample-id>.Multiple_Metrics.quality_distribution_metrics

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The `picard` Docker image used was from [Broad Institute](https://hub.docker.com/r/broadinstitute/picard).