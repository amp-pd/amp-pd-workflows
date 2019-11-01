# RNAQuantification with salmon quant

This workflow runs the `salmon quant` command on a list of paired FASTQ files
(https://combine-lab.github.io/salmon/).

Inputs:
- Per-sample:
  - Sample name
  - R1 FASTQ file list
  - R2 FASTQ file list

- Reference:
  - Gene map (GTF) file.
  - Salmon index directory tar file (.tar.gz)

- VM configuration
  - Docker image url
  - VM disk size
  - VM memory
  - Runtime zones, ex: "us-central1-a us-central1-b"
  - Number of times to try the workflow with a preemptible VM before
    falling back to a full-price VM.

Outputs:
- Per-sample
  - Transcript Quantification files
  - Gene Quantification files
  - Command info
  - Parameter files
  - Log files
  - Auxilliary files

## Notes
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/).
A terra.inputs.json and a terra.outputs.json file that you can directly upload to Terra is also included here.

The gene map (GTF) file used was [GENCODE](https://www.gencodegenes.org/) v29.

The Salmon index is packaged as a gzipped TAR file as WDL draft-2 does not support the WDL `Directory` type.
The Salmon index was created with the following steps:

1- Download the GENCODE v29 FASTA file:
```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz
```
2- Run "salmon index" as:
```
salmon index \
  -t gencode.v29.transcripts.fa.gz \
  --gencode \
  -i Homo_sapiens.gencode.v29.all.salmon0.11.3_gencodeOption.index \
  --type quasi \
  -k 31
```
3- Build a gzipped TAR file for the Salmon workflow to pull just a single reference file:
```
tar cvfz Homo_sapiens.gencode.v29.all.salmon0.11.3_gencodeOption.index.tar.gz Homo_sapiens.gencode.v29.all.salmon0.11.3_gencodeOption.index/
```


The `salmon` Docker image used was from [combinelab](https://combine-lab.github.io/salmon/).
