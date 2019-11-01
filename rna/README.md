The workflows in this repository can be run by Cromwell standalone or in the Terra environment.
The following section describes their use in the Terra environment:

Metadata about samples can be stored in Terra workspace data
(see https://support.terra.bio/hc/en-us/articles/360025758392-Populating-the-Data-table),
which can then be used to drive workflows. When workflows complete successfully,
Terra can update the `sample` data table (also known as `sample` "entities") with the location of their outputs.

Thus if you set each workflow to run from the `sample` entity and set the `outputs` to update the `sample` entity, you can run a progression of workflows. All of your workflow outputs will be stored in your Terra workspace bucket and the `sample` entity table will contain the paths to the outputs.

1- Upload a list of samples following instructions from the help article above. The file to upload could be as simple as a 1-column CSV file:
```
entity:sample_id
SAMPLE_ID_1
SAMPLE_ID_2
SAMPLE_ID_3
```
The `samples` table would then contain, for each sample:

- sample_id (STRING)

2- Run `update-samples-entity` workflow. The `samples` table would then contain, for each sample:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to FASTQ files in GCS)
- fastq_2 (Array of STRINGs: paths to FASTQ files in GCS)

3- Run the `rna-quantification` (Salmon) workflow. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to fastq files in GCS)
- fastq_2 (Array of STRINGs: paths to fastq files in GCS)
- salmon_quant_sf (STRING: path to Salmon output quant.sf file in GCS)
- salmon_quant_genes_sf (STRING: path to Salmon output quant.genes.sf file in GCS)
- [other Salmon output files]

4- Run the `rna-alignment` (STAR) workflow. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to FASTQ files in GCS)
- fastq_2 (Array of STRINGs: paths to FASTQ files in GCS)
- salmon_quant_sf (STRING: path to Salmon output quant.sf file in GCS)
- salmon_quant_genes_sf (STRING: path to Salmon output quant.genes.sf file in GCS)
- [other Salmon output files]
- STAR_bam (STRING: path to STAR output BAM in GCS)
- [other STAR output files]

5- Run the `rna-bam-sort-and-index` (samtools) workflow. The `samples` entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to FASTQ files in GCS)
- fastq_2 (Array of STRINGs: paths to FASTQ files in GCS)
- salmon_quant_sf (STRING: path to Salmon output quant.sf file in GCS)
- salmon_quant_genes_sf (STRING: path to Salmon output quant.genes.sf file in GCS)
- [other Salmon output files]
- STAR_bam (STRING: path to STAR [unsorted] output BAM in GCS)
- [other STAR output files]
- samtools_bam (STRING: path to samtools [sorted] output BAM in GCS)
- [other samtools output files]

6- Run the `rna-summarization` (featureCounts) workflow. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to FASTQ files in GCS)
- fastq_2 (Array of STRINGs: paths to FASTQ files in GCS)
- salmon_quant_sf (STRING: path to Salmon output quant.sf file in GCS)
- salmon_quant_genes_sf (STRING: path to Salmon output quant.genes.sf file in GCS)
- [other Salmon output files]
- STAR_bam (STRING: path to STAR [unsorted] output BAM in GCS)
- [other STAR output files]
- samtools_bam (STRING: path to samtools [sorted] output BAM in GCS)
- [other samtools output files]
- featureCounts_files (Array of STRINGS: paths to featureCounts output files)

7- Run the `rna-collect-rna-seq-metrics` and the `rna-collect-multiple-metrics` workflows. The samples entity would then contain:

- sample_id (STRING)
- fastq_1 (Array of STRINGs: paths to FASTQ files in GCS)
- fastq_2 (Array of STRINGs: paths to FASTQ files in GCS)
- salmon_quant_sf (STRING: path to Salmon output quant.sf file in GCS)
- salmon_quant_genes_sf (STRING: path to Salmon output quant.genes.sf file in GCS)
- [other Salmon output files]
- STAR_bam (STRING: path to STAR [unsorted] output BAM in GCS)
- [other STAR output files]
- samtools_bam (STRING: path to samtools [sorted] output BAM in GCS)
- [other samtools output files]
- featureCounts_files (Array of STRINGS: paths to featureCounts output files)
- rnaseq_metrics_path (STRING: path to RNA_Metrics file in GCS)
- multiple_metrics_path (Array of STRINGS: paths to Multiple_metrics output files in GCS)
