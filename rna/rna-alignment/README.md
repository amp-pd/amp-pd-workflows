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
  - Runtime zones, ex: "us-central1-a us-central1-b"
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
A sample inputs.json file is included here with values derived from running workflows for AMP PD on [Terra](https://app.terra.bio/). By default, the workflow will run on preemptible machines, but won't let the task run over 24 hours. To disable preemptible machines for large samples, set `RNAAlignment.preemptible_tries` to 0 and `RNAAlignment.star_timeout_hours` to a very large value.

The STAR index is packaged as a gzipped TAR file as WDL draft-2 does not support the WDL `Directory` type.
The STAR index was created with the following steps:

1- Download the GENCODE v29 FASTA and file from:

 ```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
```
2- Unzip the reference files:
```
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v29.primary_assembly.annotation.gtf.gz
```
3- Run "STAR ... genomeGenerate" as:

```
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir STAR_genome_GencodeV29_oh125 \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.v29.primary_assembly.annotation.gtf \
  --genomeSAindexNbases 14 \
  --sjdbOverhang 125
```
4- Build a gzipped TAR file:
```
   tar cvfz STAR_genome_GencodeV29_oh125.tar.gz STAR_genome_GencodeV29_oh125/
```

The `star` Docker image used was from [Alex Dobin](https://hub.docker.com/r/alexdobin/star/).
