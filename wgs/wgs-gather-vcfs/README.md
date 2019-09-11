# WGSGatherVcfs with GATK

This workflow runs the GATK's `GatherVcfsCloud` on a GCS directory containing sharded VCF files from a joint genotyping workflow such as the one here: https://github.com/gatk-workflows/gatk4-germline-snps-indels.

Each VCF file corresponds to an interval from the joint genotyping intervals file. That file is a required input to this workflow.

The output of the workflow is 1 VCF file per chromosome requested.

Inputs:
- intervals_file: GCS path of joint genotyping workflow's output intervals file.
- chromosomes: An array of chromosome string values, ["1", "2", "21", "22", "X", "Y"]
- shards_to_ignore: An array of numbers from the intervals file to ignore.
                    Some intervals from joint genotyping may have zero variants called.
                    List these shards here to prevent GatherVcfs from raising an error.

- input_directory: GCS path of directory containing sharded vcf files
- input_file_prefix: Input VCF files look like "${input_file_prefix}.n.${input_file_suffix}"
- input_file_suffix: Input VCF files look like "${input_file_prefix}.n.${input_file_suffix}"

- VM configuration
  - Docker image url
  - VM disk size needed for chromosome 1
  - VM memory
  - Num CPU cores
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per chromosome 
  - <chromosome>.vcf.gz
  - <chromosome>.vcf.idx
  - <chromosome>.vcf.gz.tbi

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).

A run on ~4000 samples was observed to have larger chromosomes take ~35 hours, so by default the workflow will run on non-preemptible machines.

To enable preemptible machines for smaller sample sets, set `GatherVcfsCloud.preemptible_tries` to the number of preemptible retry attempts desired.

The `gatk` Docker image used was from [broadinstitute](https://hub.docker.com/u/broadinstitute/).
