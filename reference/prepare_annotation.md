# Prepare Annotation

Prepare annotation to be imported as `rowRanges` and `rowData` for both
Genes, Transcripts and Exons based Position Annotation Table. From a GTF
or GFF3 annotation file.

## Usage

``` r
prepare_annotation(file_path, file_type = c("gtf", "gff"))
```

## Arguments

- file_path:

  Path to annotation file.

- file_type:

  Character indicating the type of file to download. One of `"gtf"` or
  `"gff"`. Defaults to `"gtf"`.
