# Prepare Exon based Position Annotation Table

Prepare Exon based Position Annotation Table

## Usage

``` r
prepare_exon_annotation(gene_name, file_path, file_type = c("gff", "gtf"))
```

## Arguments

- gene_name:

  String or vector of gene names to extract.

- file_path:

  Path to annotation file.

- file_type:

  A character string specifying the type of file to download. Valid
  options are `"gff"`, `"gtf"`, `"fasta"` or `"genome_fasta"`. Defaults
  to `"gff"`. **Note:** `"fasta"` refers to the transcriptome FASTA
  file. `"genome_fasta"` refers to the whole genome sequence FASTA file.
