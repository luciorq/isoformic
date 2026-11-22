# Put the quantification uncertainty into the counts

Reads that map to several isoforms make the counts of an isoform more
variable than a Poisson count, by an amount that differs between
transcripts. Baldoni et al. (2024) estimate this overdispersion for
every transcript from the Salmon replicates and divide it out of the
counts, with the aim that tests built for counts (`edgeR`, `limma`) keep
their error rate without a model of their own for the read ambiguity.
The scaled counts are not integers. In a genome-wide null of control
placentas (see the article "Is quantification uncertainty worth
checking?") dividing the counts removed the small-sample concentration
of nominal false calls in uncertain isoforms (3 and 5 samples per
group), but it did not lower the false discoveries at q \< 0.05 (on
average 1.21 against 0.20 per 3 v 3 split, and 0.80 against 0.095 per 7
v 7 split, most of them from one split). With known changes added at 7
samples per group it found fewer of them (21,349 against 31,886 true
calls), so it costs power.

- `inferential_overdispersion()` computes the overdispersion of every
  transcript from the replicates stored in the object, with the
  estimator of
  [`edgeR::catchSalmon()`](https://rdrr.io/pkg/edgeR/man/catchSalmon.html)
  (checked against edgeR 4.8): the Pearson statistic of a quasi-Poisson
  fit to the replicates of each sample, pooled over the samples where
  the transcript has a non-zero mean, moderated towards the median with
  a prior of `prior_df` degrees of freedom and floored at one. The prior
  comes from every transcript of the index, as in `catchSalmon()`, also
  when
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  keeps only the replicates of the transcripts with a mean count of at
  least `min_count`:
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  records the mean and variance of the replicates of every transcript,
  and the overdispersion of a transcript without stored replicates is
  estimated from them as well. Replicates without these moments (for
  example rebuilt with
  [`as_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  without its `moments` argument) give estimates from the stored
  transcripts only; transcripts without replicates then get the prior.

- `scale_counts()` divides an assay by the overdispersion and stores the
  result as a new assay, ready for a count-based test (the counts
  divided by the `Overdispersion` column of `catchSalmon()` in Baldoni
  et al. 2024).
  [`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
  divides the counts in the same way when the object holds replicates.

The overdispersion complements the stability columns of
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md):
it changes the test, while those columns describe the quantification
behind each call.

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

  Not used; reserved for methods.

- assay:

  Name of the counts assay to scale.

- name:

  Name of the new assay.

## Value

`inferential_overdispersion()`: a tibble with `transcript_id`,
`overdispersion` and `df` (degrees of freedom of the estimate, the
number of samples where the transcript has a non-zero mean times the
number of replicates minus one), with the prior in the `prior`
attribute. `scale_counts()`: the object with the scaled assay added.

## References

Baldoni PL, Chen Y, Hediyeh-Zadeh S, Liao Y, Dong X, Ritchie ME, Shi W,
Smyth GK (2024). Dividing out quantification uncertainty allows
efficient assessment of differential transcript expression with edgeR.
*Nucleic Acids Research* 52(3):e13.

## See also

[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)

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
    condition = rep(c("preeclampsia", "control"), each = 3)
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
#> 1 ENST00000367468.10           4.05   114
#> 2 ENST00000681605.1            3.30   114
#> 3 ENST00000680451.1            3.94   114
#> 4 ENST00000490885.6            3.00   114
#> 5 ENST00000559627.1            2.42   114
#> 6 ENST00000466691.1            1       95
```
