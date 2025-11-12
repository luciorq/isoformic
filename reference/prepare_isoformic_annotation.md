# Write Feature Annotation to Parquet Files

This function reads an annotation file and parse feature annotation to
Parquet files each level of required feature (i.e. gene, transcript, and
exon).

## Usage

``` r
prepare_isoformic_annotation(
  input_path,
  output_path = NULL,
  file_type = c("gff")
)
```

## Arguments

- input_path:

  Character string specifying the path to the input GFF file.

- output_path:

  Character string specifying the path to the output directory where
  Parquet files are written. If `NULL` or an empty string, the cache
  directory will be used.

- file_type:

  Character string specifying the type of the input file. Currently,
  only "gff" is supported (default is "gff").

## Value

Invisible path to the created Parquet file.
