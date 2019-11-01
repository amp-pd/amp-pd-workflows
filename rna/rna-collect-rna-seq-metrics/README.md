# RNACollectRnaSeqMetrics with Picard

This workflow runs [picard CollectRnaSeqMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectRnaSeqMetrics.php) on a BAM file.

Inputs:
- Per-sample:
  - Sample name - Name of the sample associated with the BAM - used in names of outputs.
  - BAM file - Path to a BAM file.
  - Ref flat file - Gene annotations in refFlat form.
  - Ribosomal intervals file - Location of rRNA sequences in genome, in interval_list mat.
- VM configuration
  - Docker image url
  - VM disk size
  - VM memory
  - Num CPU cores
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample:
  - <sample-id>.RNA_Metrics

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

The `picard` Docker image used was from [Broad Institute](https://hub.docker.com/r/broadinstitute/picard).