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
  The samples are read one at a time, so memory holds the result and the
  replicates of a single sample. When the samples have different numbers
  of replicates, the smallest number is kept. The mean and variance over
  the kept replicates of every transcript, including those left out by
  `transcripts`, are recorded in the `moments` attribute:
  [`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  estimates its prior from them, so it does not depend on which
  transcripts are kept.

- `equivalence_support()` reads the equivalence classes (`--dumpEq`,
  `aux_info/eq_classes.txt.gz`) and reports, for every transcript, how
  many of its fragments are compatible with no other transcript.

- `salmon_run_qc()` reads `aux_info/meta_info.json` and reports the
  mapping rate, the library type, whether the optimiser converged and
  how many replicates were drawn.

- `inferential_replicates_table()` and `as_inferential_replicates()`
  convert the list of replicates to and from a long table (one row per
  transcript, sample and replicate), which can be stored as Parquet with
  [`write_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  and read back with
  [`read_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  for example to exchange them with other tools.
  `inferential_moments_table()` converts the `moments` attribute to a
  table in the same way (one row per transcript and sample); pass it to
  the `moments` argument of `as_inferential_replicates()` to restore it.
  The long table takes several times the memory of the replicates: to
  save the replicates of an object, use
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md),
  which writes one wide file per draw. Both apply the `selection` of the
  replicates of a subset `x[i, j]` (see
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)):
  the replicate table holds the kept transcripts and samples only, the
  moment table the kept samples (the moments keep every transcript). The
  selection is an attribute of the unmodified `x@metadata$inf_reps`: a
  list operation such as
  [`lapply()`](https://rdrr.io/r/base/lapply.html) or `[` drops it, and
  the list then holds the draws of the parent. Pass the object as `x`
  (alone, or with such a list as `reps`) to read the draws of its
  transcripts and samples in every case.

The functions also read Oarfish output (see
[`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)):
its bootstrap replicates (`.infreps.pq`), its unique and ambiguous read
counts (`.ambig_info.tsv`, `n_classes` is then `NA`) and its run
metadata (`.meta_info.json`). One call reads the output of one tool.

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

inferential_replicates_table(reps = NULL, x = NULL)

inferential_moments_table(reps = NULL, x = NULL)

as_inferential_replicates(table, samp_type = "gibbs", moments = NULL)
```

## Arguments

- paths:

  Character vector of Salmon output directories (containing `quant.sf`)
  or paths to `quant.sf` files, or Oarfish output prefixes. Names, when
  present, are the sample IDs; otherwise the directory name (or the
  prefix) is used.

- transcripts:

  Optional transcript IDs to keep. Reading every transcript of a human
  annotation for many samples takes gigabytes of memory, so
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  only reads the expressed ones.

- n_draws:

  Number of replicates to keep (the first ones). `NULL` keeps all of
  them; a number above the replicates available keeps all of them, with
  a message.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- min_percent_mapped:

  Samples below this mapping rate (percent) get the `low_mapping_rate`
  flag in `salmon_run_qc()`.

- reps:

  A list of transcript x sample count matrices, one per replicate. In
  `inferential_replicates_table()` and `inferential_moments_table()`, it
  defaults to `x@metadata$inf_reps` when `x` is given.

- x:

  Optional
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).
  When given, the replicates are read as those of its transcripts and
  samples, matched by ID, even when `reps` lost the `selection` of a
  subset.

- table:

  A long table with `transcript_id`, `sample_id`, `draw` and `count`
  columns, as returned by `inferential_replicates_table()`.

- samp_type:

  Replicate type recorded on the result, `"gibbs"` or `"bootstrap"`.

- moments:

  Optional table of `inferential_moments_table()` (columns
  `transcript_id`, `sample_id`, `n_draws`, `mean` and `var`) to record
  as the `moments` attribute of the result.

## Value

- `read_inferential_replicates()`: a list of numeric matrices
  (transcripts by samples, counts) with the replicate type (`"gibbs"` or
  `"bootstrap"`) in the `samp_type` attribute and the moments of every
  transcript in the `moments` attribute (a list with the `mean` and
  `var` matrices, transcripts by samples, and the number of replicates
  `n`).

- `inferential_replicates_table()`: a tibble with `transcript_id`,
  `sample_id`, `draw` and `count`. `as_inferential_replicates()` is its
  inverse. `inferential_moments_table()`: a tibble with `transcript_id`,
  `sample_id`, `n_draws`, `mean` and `var`.

- `equivalence_support()`: a tibble with one row per transcript:
  `transcript_id`, `n_classes` (mean number of equivalence classes that
  contain the transcript, per sample), `n_fragments` and
  `n_unique_fragments` (summed over samples), `unique_fraction` (share
  of the fragments of the transcript that map to it alone) and
  `n_samples` (samples with equivalence classes).

- `salmon_run_qc()`: a tibble with one row per sample: `sample_id`,
  `salmon_version`, `library_type`, `detected_library_type`,
  `num_processed`, `num_mapped`, `percent_mapped`, `replicate_type`,
  `num_replicates`, `optimiser`, `em_converged`, `seq_bias_correct`,
  `gc_bias_correct`, `pos_bias_correct`, `frag_length_mean`,
  `frag_length_sd` and `qc_flag` (`low_mapping_rate`, `no_replicates`
  and `em_not_converged`, separated by `;`, or `NA`). For Oarfish output
  the columns are `sample_id`, `num_processed`, `num_mapped`,
  `percent_mapped`, `replicate_type`, `num_replicates`, `optimiser`,
  `alignment_source`, `prob_model` and `qc_flag`, and the `tool`
  attribute is `"oarfish"`.

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
#> 1 SRR11498039          0.324 low_mapping_rate
#> 2 SRR11498040          0.422 low_mapping_rate
#> 3 SRR11498041          0.486 low_mapping_rate
#> 4 SRR11498077          0.449 low_mapping_rate
#> 5 SRR11498078          0.551 low_mapping_rate
#> 6 SRR11498079          0.507 low_mapping_rate
```
