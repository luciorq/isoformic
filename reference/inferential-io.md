# Read the inferential uncertainty that Salmon writes

Salmon quantifies how well each transcript abundance is determined by
the reads, and writes that information next to `quant.sf`:

- `read_inferential_replicates()` reads the Gibbs (`--numGibbsSamples`)
  or bootstrap (`--numBootstraps`) replicates from `aux_info/bootstrap/`
  (the format is identical in Salmon 1.x and in the Rust rewrite, Salmon
  2). It returns a list with one transcript x sample matrix of counts
  per replicate, the layout used by
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  and
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md).

- `equivalence_support()` reads the equivalence classes (`--dumpEq`,
  `aux_info/eq_classes.txt.gz`) and reports, for every transcript, how
  many of its fragments are compatible with no other transcript.

- `salmon_run_qc()` reads `aux_info/meta_info.json` and reports the
  mapping rate, the library type, whether the optimiser converged and
  how many replicates were drawn.

- `inferential_replicates_table()` and `as_inferential_replicates()`
  convert the list of replicates to and from a long table (one row per
  transcript, sample and replicate), which can be stored as Parquet.

## Usage

``` r
read_inferential_replicates(
  paths,
  transcripts = NULL,
  n_draws = NULL,
  ignore_after_bar = TRUE
)

equivalence_support(paths, transcripts = NULL, ignore_after_bar = TRUE)

salmon_run_qc(paths, min_percent_mapped = 50)

inferential_replicates_table(reps)

as_inferential_replicates(table, samp_type = "gibbs")
```

## Arguments

- paths:

  Character vector of Salmon output directories (containing `quant.sf`)
  or paths to `quant.sf` files. Names, when present, are the sample IDs;
  otherwise the directory name is used.

- transcripts:

  Optional transcript IDs to keep. Reading every transcript of a human
  annotation for many samples takes gigabytes of memory, so
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  only reads the expressed ones.

- n_draws:

  Number of replicates to keep (the first ones). `NULL` keeps all of
  them.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- min_percent_mapped:

  Samples below this mapping rate (percent) get the `low_mapping_rate`
  flag in `salmon_run_qc()`.

- reps:

  A list of transcript x sample count matrices, one per replicate.

- table:

  A long table with `transcript_id`, `sample_id`, `draw` and `count`
  columns, as returned by `inferential_replicates_table()`.

- samp_type:

  Replicate type recorded on the result, `"gibbs"` or `"bootstrap"`.

## Value

- `read_inferential_replicates()`: a list of numeric matrices
  (transcripts by samples, counts) with the replicate type (`"gibbs"` or
  `"bootstrap"`) in the `samp_type` attribute.

- `inferential_replicates_table()`: a tibble with `transcript_id`,
  `sample_id`, `draw` and `count`. `as_inferential_replicates()` is its
  inverse.

- `equivalence_support()`: a tibble with one row per transcript:
  `transcript_id`, `n_classes` (mean number of equivalence classes per
  sample), `n_fragments` and `n_unique_fragments` (summed over samples),
  `unique_fraction` (share of the fragments of the transcript that map
  to it alone) and `n_samples`.

- `salmon_run_qc()`: a tibble with one row per sample.

## Examples

``` r
quant_dirs <- list.files(
  system.file("extdata", "mini_quants", package = "isoformic"),
  full.names = TRUE
)
names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
# The bundled quantifications were made without replicates
salmon_run_qc(quant_dirs)[, c("sample_id", "percent_mapped", "qc_flag")]
#> # A tibble: 6 × 3
#>   sample_id   percent_mapped qc_flag      
#>   <chr>                <dbl> <chr>        
#> 1 SRR11498039            100 no_replicates
#> 2 SRR11498040            100 no_replicates
#> 3 SRR11498041            100 no_replicates
#> 4 SRR11498079            100 no_replicates
#> 5 SRR11498080            100 no_replicates
#> 6 SRR11498081            100 no_replicates
```
