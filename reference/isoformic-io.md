# Save and Restore IsoformicExperiment Objects as Parquet

`write_isoformic()` writes an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
to a directory of Parquet files and `read_isoformic()` restores it. The
layout is:

    <path>/
      experiment.parquet        key/value: experiment_name, annot_path, ...
      col_data.parquet          sample metadata
      metadata.parquet          key/value: scalar entries of `metadata`
      tx_type_palette.parquet   transcript_type, color
      assays/<name>.parquet     transcript_id + one column per sample
      gene_assays/<name>.parquet  only explicitly set gene assays
      dea/det.parquet, dea/deg.parquet, dea/cutoffs.parquet
      gsea/<name>.parquet
      annotation/               copy of the annotation directory (optional)

With `copy_annotation = TRUE` the annotation directory is copied inside
the experiment directory so that the result is self-contained and
portable. Otherwise `annot_path` is stored as an absolute path and must
still exist when the object is read back (it can be overridden with the
`annot_path` argument of `read_isoformic()`).

## Usage

``` r
write_isoformic(
  self,
  path = NULL,
  overwrite = FALSE,
  copy_annotation = FALSE,
  ...
)

read_isoformic(path, annot_path = NULL)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- path:

  Directory to write to or read from. For `write_isoformic()` defaults
  to the object's `data_path`, or to
  `get_isoformic_cache("experiments", experiment_name)`.

- overwrite:

  Logical, replace an existing directory.

- copy_annotation:

  Logical, copy the annotation directory inside `path`.

- ...:

  Passed on to methods.

- annot_path:

  Optional annotation directory overriding the stored one when reading.

## Value

`write_isoformic()` returns the object with `data_path` set, invisibly.
`read_isoformic()` returns an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Examples

``` r
iso <- example_isoformic()
dir <- file.path(tempdir(), "iso_example")
iso <- write_isoformic(iso, dir, overwrite = TRUE, copy_annotation = TRUE)
list.files(dir, recursive = TRUE)
#>  [1] "annotation/annot_data_exons.parquet"      
#>  [2] "annotation/annot_data_genes.parquet"      
#>  [3] "annotation/annot_data_sequences.parquet"  
#>  [4] "annotation/annot_data_transcripts.parquet"
#>  [5] "annotation/annot_metadata.parquet"        
#>  [6] "annotation/annot_seq_regions.parquet"     
#>  [7] "assays/abundance.parquet"                 
#>  [8] "assays/counts.parquet"                    
#>  [9] "col_data.parquet"                         
#> [10] "dea/deg.parquet"                          
#> [11] "dea/det.parquet"                          
#> [12] "experiment.parquet"                       
#> [13] "metadata.parquet"                         
#> [14] "tx_type_palette.parquet"                  
iso2 <- read_isoformic(dir)
identical(assay(iso), assay(iso2))
#> [1] TRUE
```
