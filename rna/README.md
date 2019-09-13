The workflows in this repository can be run by Cromwell standalone or in the Terra environment.
The following section describes their use in the Terra environment:

Metadata about samples can be stored in Terra workspace data
(see https://support.terra.bio/hc/en-us/articles/360025758392-Populating-the-Data-table),
which can then be used to drive workflows. When workflows complete successfully,
Terra can update the `sample` data table (also known as `sample` "entities") with the location of their outputs.

Thus if you set each workflow to run from the `sample` entity and set the `outputs` to update the `sample` entity, you can run a progression of workflows such as:

1- Run `update-samples-entity` workflow. The `samples` table would then contain, for each sample:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to fastq files in GCS)
- fastq_2 (Array of STRINGs: paths to fastq files in GCS)

2- Run the `rna-alignment` (STAR) workflow. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to fastq files in GCS)
- fastq_2 (Array of STRINGs: paths to fastq files in GCS)
- STAR_bam (STRING: path to STAR output BAM in GCS)
- [other STAR output files]

3- Run the `rna-bam-sort-and-index` (samtools) workflow. The `samples` entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to fastq files in GCS)
- fastq_2 (Array of STRINGs: paths to fastq files in GCS)
- STAR_bam (STRING: path to STAR [unsorted] output BAM in GCS)
- [other STAR output files]
- samtools_bam (STRING: path to samtools [sorted] output BAM in GCS)
- [other samtools output files]

4- Run the `rna-summarization` (featureCounts) workflow. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to fastq files in GCS)
- fastq_2 (Array of STRINGs: paths to fastq files in GCS)
- STAR_bam (STRING: path to STAR [unsorted] output BAM in GCS)
- [other STAR output files]
- samtools_bam (STRING: path to samtools [sorted] output BAM in GCS)
- [other samtools output files]
- featureCounts_files (Array of STRINGS: paths to featureCounts output files)

All of your workflow outputs will be stored in your Terra workspace bucket and the `sample` entity table will contain the paths to the outputs.