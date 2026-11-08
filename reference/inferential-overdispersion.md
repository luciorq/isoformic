# Put the quantification uncertainty into the counts

Reads that map to several isoforms make the counts of an isoform more
variable than a Poisson count, by an amount that differs between
transcripts. Baldoni et al. (2024) estimate this overdispersion for
every transcript from the Salmon replicates and divide it out of the
counts, so that tests built for counts (`edgeR`, `limma`, `DESeq2`) keep
their error rate without a model of their own for the read ambiguity.

- `inferential_overdispersion()` computes the overdispersion of every
  transcript from the replicates stored in the object, as
  [`edgeR::catchSalmon()`](https://rdrr.io/pkg/edgeR/man/catchSalmon.html)
  does: the Pearson statistic of a quasi-Poisson fit to the replicates
  of each sample, pooled over the samples, moderated towards the median
  with a prior of `prior_df` degrees of freedom and floored at one. The
  estimates use the transcripts whose replicates are in the object
  (those with a mean count of at least `min_count` in
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)),
  so the prior differs from the one `catchSalmon()` computes over every
  transcript of the index; transcripts without replicates get the prior.

- `scale_counts()` divides an assay by the overdispersion and stores the
  result as a new assay, ready for a count-based test.

The overdispersion complements the confidence classes of
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md):
it changes the test, while the classes describe how far each call can be
trusted.

## Usage

``` r
inferential_overdispersion(self, prior_df = 3, ...)

scale_counts(self, assay = "counts", name = paste0(assay, "_scaled"), ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  with inferential replicates (see
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)).

- prior_df:

  Degrees of freedom of the prior (3 in `edgeR`).

- ...:

  Unused.

- assay:

  Name of the counts assay to scale.

- name:

  Name of the new assay.

## Value

`inferential_overdispersion()`: a tibble with `transcript_id`,
`overdispersion` and `df` (degrees of freedom of the estimate), with the
prior in the `prior` attribute. `scale_counts()`: the object with the
scaled assay added.

## References

Baldoni PL, Chen Y, Hediyeh-Zadeh S, Liao Y, Dong X, Ritchie ME, Shi W,
Smyth GK (2024). Dividing out quantification uncertainty allows
efficient assessment of differential transcript expression with edgeR.
*Nucleic Acids Research* 52(3):e13.

## See also

[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)

## Examples

``` r
extdata <- system.file("extdata", package = "isoformic")
quant_dirs <- list.files(file.path(extdata, "mini_quants"), full.names = TRUE)
names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
iso <- read_salmon(
  quant_dirs,
  annot_path = file.path(extdata, "mini_annotation"),
  col_data = data.frame(
    sample_id = names(quant_dirs),
    condition = rep(c("control", "treatment"), each = 3)
  )
)
iso@metadata$inf_reps <- as_inferential_replicates(read_parquet(
  file.path(extdata, "mini_draws.parquet")
))
iso <- scale_counts(iso)
assay_names(iso)
#> [1] "abundance"     "counts"        "length"        "counts_scaled"
head(inferential_overdispersion(iso))
#> # A tibble: 6 × 3
#>   transcript_id      overdispersion    df
#>   <chr>                       <dbl> <dbl>
#> 1 ENST00000367468.10           2.96   114
#> 2 ENST00000681605.1            5.30   114
#> 3 ENST00000680451.1            3.93   114
#> 4 ENST00000490885.6            1      114
#> 5 ENST00000559627.1            1.58   114
#> 6 ENST00000466691.1            1      114
```
