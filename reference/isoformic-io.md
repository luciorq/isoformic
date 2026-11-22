# Save and Restore IsoformicExperiment Objects as Parquet

`write_isoformic()` writes an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
to a directory of Parquet files and `read_isoformic()` restores it. The
layout is:

    <path>/
      experiment.parquet        key/value: experiment_name, annot_path,
                                isoformic_version, written_at
      col_data.parquet          sample metadata
      metadata.parquet          key/value: scalar entries of `metadata`
      tx_type_palette.parquet   transcript_type, color
      assays/<name>.parquet     transcript_id + one column per sample
      gene_assays/<name>.parquet  only explicitly set gene assays
      dea/det.parquet, dea/deg.parquet, dea/cutoffs.parquet
      dea/uncertainty.parquet, dea/ambiguity.parquet, dea/confidence.parquet
      gsea/<name>.parquet       list columns joined with ";"
      inf_reps/draw_0001.parquet, ...  one draw: transcript_id + samples
      inf_reps/moments/mean.parquet, var.parquet  replicate moments
      annotation/               copy of the annotation directory (optional)

Each file is written only when the object holds that part. Entries of
`metadata` that are not single values are dropped with a warning.

The inferential replicates (`self@metadata$inf_reps`, see
[`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
and
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
are written to `inf_reps/`: one Parquet file per draw
(`draw_0001.parquet`, ...) with `transcript_id` and one column per
sample, and, when the replicates have a `moments` attribute, its `mean`
and `var` matrices in `inf_reps/moments/`. Their number, their type (the
`samp_type` attribute) and the number of draws behind the moments are
recorded in `metadata.parquet` under the keys `inf_reps_n_draws`,
`inf_reps_samp_type` and `inf_reps_moments_n`; when replicates are
written, these keys replace entries of `metadata` with those names, with
a warning. One draw is written, and read, at a time, so memory holds the
replicates and a single draw on top of them. The replicates of a subset
(`x[i, j]` shares the draws of `x` and records a `selection`, see
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md))
are written as the selected transcripts and samples only, and the
moments as the selected samples. `read_isoformic()` restores the list of
matrices with its `samp_type` and `moments` attributes, as they were
(including the storage mode of the matrices and of `n`); for a subset,
as plain matrices of the selected draws, without the `selection`. Only
these are kept: names of the draws, other attributes of the list or of
the matrices and names of their dimnames are dropped, with a warning.
Replicates that are not numeric matrices with transcript and sample
names are not written, with a warning. Directories written by earlier
versions have no `inf_reps/` and are read without replicates. Use
`add_uncertainty(keep_draws = FALSE)`, or remove `metadata$inf_reps`, to
save an object without its replicates.

With `copy_annotation = TRUE` the annotation directory is copied inside
the experiment directory so that the result is self-contained and
portable. The cached ORF table (`annot_data_orfs.parquet`) is left out
of the copy and computed again when needed. An annotation that lies
inside the directory being replaced is always copied, and so is an
annotation directory under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) (such as the private
copy made by
[`example_isoformic()`](https://luciorq.github.io/isoformic/reference/example_isoformic.md))
when `path` is not under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) itself: the
temporary directory is removed when the R session ends. A message says
so. Otherwise `annot_path` is stored as an absolute path and must still
exist when the object is read back (it can be overridden with the
`annot_path` argument of `read_isoformic()`). Without that argument,
`read_isoformic()` uses an `annotation/` directory inside `path`
whenever there is one.

The new directory is built next to `path` and moved in place at the end,
so a failed write leaves an existing experiment as it was. With
`overwrite = TRUE`, a non-empty `path` is replaced only when it holds an
`experiment.parquet` file.

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

  Not used; reserved for methods.

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
identical(tx_assay(iso), tx_assay(iso2))
#> [1] TRUE
```
