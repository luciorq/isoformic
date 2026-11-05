# Test differential transcript usage

`test_usage()` tests whether the share of a gene's expression carried by
each isoform differs between groups of samples, with one of three
established frameworks:

- `"diffsplice"`: `edgeR` quasi-likelihood fits and `diffSplice`, the
  pipeline of the 2025 *Nucleic Acids Research* paper below. Fast, and
  the one designed for counts divided by the read-ambiguity
  overdispersion.

- `"satuRn"`: quasi-binomial models of transcript usage (Gilis et al.,
  2021), suited to many samples and to single-cell data.

- `"DEXSeq"`: negative binomial models of exon (here transcript) usage
  (Anders et al. 2012), the slowest and often the most conservative.

Recent benchmarks of usage tests (2024 and 2025) rank these among the
best; none dominates, so the backend is an argument. With inferential
replicates in the object (see
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)),
the counts are divided by the overdispersion that the read ambiguity
adds (see
[`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md))
before the test (`scale = TRUE`).

Transcripts are kept when at least `min_count` reads and a share of
`min_prop` of the gene fall on them in as many samples as the smallest
group, and genes need two such transcripts.

The gene-level p-value combines the transcript p-values of a gene with
Simes' method, and both levels are adjusted with the Benjamini-Hochberg
procedure within each contrast, so that the three backends report the
same columns.

With more than two groups, every contrast of `contrasts` is tested and
the results are stacked, one block per contrast. `"diffsplice"` and
`"satuRn"` fit all the groups once and test each contrast on that fit,
so the dispersions use every sample; `"DEXSeq"` compares a model with
and without a group effect, so it is fitted on the samples of the two
groups of each contrast.

`contrasts = "omnibus"` gives one p-value per transcript for a
difference among any of the groups, with `contrast = "omnibus"` and no
fold-change (`usage_log2FC` is `NA`; look at the contrasts for the
direction). `"DEXSeq"` compares its models with and without the group
effect on all the samples. `edgeR` and `satuRn` test one contrast at a
time, so for them the p-value is the Simes combination of the contrasts
of every level against the first: the groups can only differ if some
level differs from the first, and those contrasts share the reference,
so they are positively correlated and Simes' test keeps its error rate.

## Usage

``` r
test_usage(
  self,
  group,
  levels = NULL,
  contrasts = NULL,
  method = "diffsplice",
  assay = "counts",
  scale = NULL,
  min_count = 10,
  min_prop = 0.1,
  workers = 1L,
  gene_name = NULL,
  ...
)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- group:

  Column of
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  with the groups to compare.

- levels:

  Levels of `group` to use, reference first. Samples of other groups are
  ignored. `NULL` takes all the levels of `group`.

- contrasts:

  The comparisons to test: `NULL` (default) compares every level with
  the first one; `"pairwise"` tests every pair of levels (the later
  level against the earlier one); `"omnibus"` tests whether usage
  differs among any of the groups (see Details); otherwise a list of
  `c(numerator, denominator)` pairs of levels, or one such pair. With
  two levels the default is the single contrast second against first.

- method:

  `"diffsplice"`, `"satuRn"` or `"DEXSeq"`. The backend must be
  installed.

- assay:

  Counts assay to test.

- scale:

  Logical, divide the counts by the overdispersion of the replicates.
  `NULL` (default) does it when the object has replicates. In a
  read-split experiment (technical replicates of the same libraries, see
  the uncertainty article) this removed most false positives of
  `"diffsplice"` and `"DEXSeq"` and made their calls reproducible, so
  keep it on.

- min_count, min_prop:

  Filters described above.

- workers:

  Number of workers for `satuRn` and `DEXSeq` (parallel on Unix-like
  systems).

- gene_name:

  Gene names or IDs to keep in the result. The filter is applied after
  the test.

- ...:

  Unused.

## Value

A tibble with one row per tested transcript and contrast: `gene_id`,
`gene_name`, `transcript_id`, `transcript_name`, `transcript_type`,
`contrast` (`"numerator vs denominator"`), `usage_log2FC` (the change in
usage on the log2 scale: for `"diffsplice"` the transcript against the
other transcripts of the gene, for `"satuRn"` the log2 odds ratio, for
`"DEXSeq"` the usage fold-change), `pvalue`, `qvalue`, `gene_pvalue`,
`gene_qvalue`, `method` and `scaled`.

## References

Dividing out quantification uncertainty enables assessment of
differential transcript usage with limma and edgeR (2025). *Nucleic
Acids Research*, doi:10.1093/nar/gkaf1305.

Gilis J et al. (2021). satuRn: scalable analysis of differential
transcript usage for bulk and single-cell RNA-sequencing applications.
*F1000Research* 10:374.

Anders S, Reyes A, Huber W (2012). Detecting differential usage of exons
from RNA-seq data. *Genome Research* 22:2008-2017.

## See also

[`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md),
[`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md),
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)

## Examples

``` r
iso <- example_isoformic()
if (requireNamespace("edgeR", quietly = TRUE)) {
  res <- test_usage(iso, group = "condition", method = "diffsplice")
  head(res)
  table(res$contrast)
}
#> 
#> control vs preeclampsia 
#>                     164 
```
