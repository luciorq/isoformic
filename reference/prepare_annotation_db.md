# Write Parquet File from GFF

This function reads a GFF file and writes its contents to a Parquet file
using DuckDB.

## Usage

``` r
prepare_annotation_db(input_path, output_path = NULL, file_type = c("gff"))
```

## Arguments

- input_path:

  Character string specifying the path to the input GFF file.

- output_path:

  Character string specifying the path to the output Parquet file. If
  `NULL` or an empty string, a temporary file will be created.

- file_type:

  Character string specifying the type of the input file. Currently,
  only "gff" is supported (default is "gff").

## Value

Invisible path to the created Parquet file.
