# Convert a SummarizedExperiment to an IsoformicExperiment Object

This function converts a `SummarizedExperiment` object to an
`IsoformicExperiment` object. It extracts the assays, row data, column
data, and metadata from the input object and uses them to create a new
`IsoformicExperiment` object.

## Usage

``` r
as_isoformic(se, annot_path, annot_type = c("gff", "annot_db"))
```

## Arguments

- se:

  A `SummarizedExperiment` object to be converted.

- annot_path:

  Path to the annotation file. This can be a GFF file or the path
  pre-built annotation database created with
  `[prepare_isoformic_annotation()]`.

- annot_type:

  Type of the annotation file provided. Options are "gff" for GFF files
  and "annot_db" for pre-built annotation databases.
