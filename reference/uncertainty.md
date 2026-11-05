# Quantify and report the uncertainty of isoform-level results

An isoform-level result is only as good as the quantification behind it.
Reads that map to several isoforms of a gene can be shared out between
them in more than one way, so the estimate of each isoform carries an
*inferential* uncertainty that is unrelated to biological variability.
Salmon measures it by drawing replicates from the posterior (Gibbs
sampling, `--numGibbsSamples`) and by writing the equivalence classes
(`--dumpEq`). These functions turn that information into columns that
can be filtered on and reported.

Every isoform gets four kinds of evidence:

- **Quantification uncertainty**: the inferential relative variance
  (`mean_inf_rv`, Zhu et al. 2019, as in `fishpond`), the variance of
  the replicates in excess of what counting noise alone would give,
  relative to the expression. Isoforms with a high value are not well
  determined by the reads. `uncertain` flags
  `mean_inf_rv >= inf_rv_cutoff`.

- **Identifiability**: the correlation between the replicates of
  isoforms of the same gene. A strong negative correlation means the
  reads cannot tell the two isoforms apart and every read taken from one
  is given to the other. Pairs at or below `cor_cutoff` are *ambiguous*
  and are joined into `ambiguity_group`s that can be tested as one
  feature.

- **Read support**: the share of the fragments of an isoform that map to
  it alone (`unique_fraction`, from the equivalence classes).

- **Stability of the fold-change**: the contrast between two groups is
  repeated on every replicate. `lfc_lo`/`lfc_hi` is the 95% interval of
  the log2 fold-change over replicates, `sign_agreement` the share of
  replicates that agree on its direction and `frac_draws_sig` the share
  in which the isoform is called differentially expressed.

`de_confidence()` combines the last three into one class per isoform:
`high` when the fold-change is stable
(`sign_agreement >= min_agreement`), the isoform is not `uncertain` and
has no ambiguous sibling; `moderate` when exactly one of these fails;
`low` when two or more fail. The class says how far a call on the
isoform can be trusted, so read it for the isoforms flagged by the
primary test (`n_det_sig`, `confident_de`): an isoform without a real
change is `moderate` or `low` simply because its fold-change is not
stable.

- `add_uncertainty(x, paths)` reads the replicates, equivalence classes
  and run metadata from the Salmon output directories `paths` (named by
  sample ID, in any order), computes the tables below and stores them in
  the object, so that they are written by
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md).
  When `group` is given it also stores the fold-change stability of that
  contrast.

- `tx_uncertainty(x)`: per-transcript table, for the transcripts whose
  draws are in the object (the expressed ones).

- `isoform_ambiguity(x)`: per-pair table of replicate correlations.

- `de_confidence(x, group)`: per-feature table of fold-change stability
  and confidence class; with `collapse_groups = TRUE`, ambiguous
  isoforms are summed into one feature before testing.

- `uncertainty_report(x)`: a list of summary tables for a methods or
  results section.

- `uncertainty_methods(x)`: a ready-to-edit paragraph that states the
  settings used.

The replicates themselves are kept in memory only (in
`metadata(x)$inf_reps`) and are not written by
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md);
the tables derived from them are. Draws are kept for the transcripts
with an estimated mean count of at least `min_count`; the others are
treated as unexpressed.

## Usage

``` r
add_uncertainty(
  self,
  paths = NULL,
  group = NULL,
  levels = NULL,
  n_draws = NULL,
  min_count = 1,
  inf_rv_cutoff = 0.5,
  cor_cutoff = -0.5,
  min_agreement = 0.95,
  equivalence = TRUE,
  run_qc = TRUE,
  keep_draws = TRUE,
  collapse_groups = FALSE,
  groups = NULL,
  ...
)

tx_uncertainty(self, gene_name = NULL, ...)

isoform_ambiguity(self, gene_name = NULL, ambiguous_only = FALSE, ...)

de_confidence(
  self,
  group = NULL,
  levels = NULL,
  collapse_groups = FALSE,
  gene_name = NULL,
  pseudocount = 5,
  inf_rv_cutoff = NULL,
  min_agreement = NULL,
  ...
)

uncertainty_report(self, ...)

uncertainty_methods(self, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- paths:

  Salmon output directories (see
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)),
  named by sample ID. `NULL` uses the replicates already stored in the
  object.

- group:

  Name of a
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  column with the sample groups.

- levels:

  The two group levels to contrast, reference first. `NULL` uses the
  levels of `group`, which must then be exactly two.

- n_draws:

  Number of replicates to read; `NULL` reads all of them.

- min_count:

  Minimum mean estimated count for a transcript to be assessed.

- inf_rv_cutoff:

  Inferential relative variance from which an isoform is called
  uncertain. `NULL` uses the value the tables were computed with.

- cor_cutoff:

  Mean replicate correlation at or below which two isoforms of a gene
  are called ambiguous. `NULL` uses the stored value.

- min_agreement:

  Share of replicates that must agree on the direction of a fold-change
  for it to be called stable. `NULL` uses the stored value.

- equivalence:

  Logical, read the equivalence classes when present, or the table
  returned by
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
  to attach read support when the draws are supplied without `paths`.

- run_qc:

  Logical, add the Salmon run metadata (see
  [`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md))
  to
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md).

- keep_draws:

  Logical, keep the replicates in the object (needed to recompute the
  tables and for
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  uncertainty).

- collapse_groups:

  Logical, sum ambiguous isoforms into one feature before testing the
  contrast.

- groups:

  External groups of transcripts that the reads cannot tell apart,
  replacing the groups built from the correlation of the replicates: a
  table with `transcript_id` and `group_id`, or the path to Terminus
  output (see
  [`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md)),
  or the cut of TreeTerminus trees (see
  [`cut_transcript_trees()`](https://luciorq.github.io/isoformic/reference/treeterminus.md)).
  The confidence classes and `collapse_groups` then follow these groups;
  the pairs of `isoform_ambiguity()` keep describing the correlation.

- ...:

  Unused.

- gene_name:

  Optional gene names or IDs to keep.

- ambiguous_only:

  Logical, `isoform_ambiguity()` returns only the ambiguous pairs.

- pseudocount:

  Pseudocount of the log2 fold-change (5, as in `swish`).

## Value

- `add_uncertainty()`: the object with `dea$uncertainty`,
  `dea$ambiguity` and, when `group` is given, `dea$confidence` filled
  in.

- `tx_uncertainty()`: a tibble with `transcript_id`, `gene_id`, names
  and type, `mean_count`, `mean_inf_rv`, `max_inf_rv`, `mean_inf_cv`,
  `uncertain`, `min_sibling_cor`, `n_ambiguous_partners`,
  `ambiguity_group` and, when equivalence classes were read,
  `n_classes`, `n_unique_fragments` and `unique_fraction`.

- `isoform_ambiguity()`: a tibble of isoform pairs of the same gene with
  `mean_cor`, `min_cor`, `n_samples` and `ambiguous`.

- `de_confidence()`: a tibble with one row per feature: `feature_id`,
  `feature_type` (`"transcript"` or `"group"`), gene and transcript
  names, the interval and stability of the fold-change, `mean_inf_rv`,
  `ambiguous`, `confidence` (an ordered factor), `reason`, the existing
  differential expression result (`det_log2FC`, `det_qvalue`,
  `n_det_sig`) and `confident_de`, which is `TRUE` when a significant
  isoform has `high` confidence.

- `uncertainty_report()`: a list with `settings`, `run`, `transcripts`
  (by transcript type), `genes`, `de` and `switches` tibbles.

- `uncertainty_methods()`: a character string.

## Details

The thresholds are heuristics, not universal constants: `inf_rv_cutoff`
(0.5), `cor_cutoff` (-0.5) and `min_agreement` (0.95) are conservative
starting points. With 30 draws, sampling noise alone puts the
inferential relative variance of a well-determined, highly expressed
transcript near 0.1, so cutoffs much below 0.5 flag noise; poorly
determined isoforms typically reach values of several units. Report the
numbers next to the classes, and rerun with other cutoffs if the data
call for it. The replicates measure the uncertainty of the
quantification only; they do not replace biological replicates and
cannot make a noisy contrast significant.

Fold-changes are computed on counts scaled to a common library size, as
`log2(mean_b + pseudocount) - log2(mean_a + pseudocount)` with a
pseudocount of 5 (the `swish` convention), and the per-replicate test is
a Welch t-test on the log scale with Benjamini-Hochberg adjustment
across features. This is a screening statistic for stability, not a
replacement for the primary test (`swish`, `DESeq2`, ...), whose results
are shown next to it. When the sample table (`col_data`) has a
`library_size` column, with the total estimated count of the full
quantification, it is used for the scaling. Set it when the object holds
only a subset of the transcripts, so that the subset does not change the
scaling.

## References

Zhu, A., Srivastava, A., Ibrahim, J.G., Patro, R., Love, M.I. (2019).
Nonparametric expression analysis using inferential replicate counts.
*Nucleic Acids Research* 47(18), e105.

## See also

[`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
[`plot_de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md),
[`isoform_switches()`](https://luciorq.github.io/isoformic/reference/diversity.md)

## Examples

``` r
quant_dirs <- list.files(
  system.file("extdata", "mini_quants", package = "isoformic"),
  full.names = TRUE
)
names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
annot <- system.file("extdata", "mini_annotation", package = "isoformic")
col_data <- data.frame(
  sample_id = names(quant_dirs),
  condition = rep(c("control", "treatment"), each = 3)
)
iso <- read_salmon(quant_dirs, annot_path = annot, col_data = col_data)
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion

# With your own data, `read_salmon(dirs, annot, inferential = TRUE)` or
# `add_uncertainty(iso, paths = dirs)` read the replicates that
# `salmon quant --numGibbsSamples 30` wrote; here bundled draws are used
draws <- as_inferential_replicates(read_parquet(
  system.file("extdata", "mini_draws.parquet", package = "isoformic")
))
iso@metadata$inf_reps <- draws
iso <- add_uncertainty(iso, group = "condition")
head(tx_uncertainty(iso))
#> # A tibble: 6 × 14
#>   transcript_id     gene_id gene_name transcript_name transcript_type mean_count
#>   <chr>             <chr>   <chr>     <chr>           <chr>                <dbl>
#> 1 ENST00000367468.… ENSG00… PTGS2     PTGS2-201       protein_coding       22.2 
#> 2 ENST00000681605.1 ENSG00… PTGS2     PTGS2-206       nonsense_media…       6.69
#> 3 ENST00000680451.1 ENSG00… PTGS2     PTGS2-205       protein_coding        4.26
#> 4 ENST00000490885.6 ENSG00… PTGS2     PTGS2-203       retained_intron       2.87
#> 5 ENST00000559627.1 ENSG00… PTGS2     PTGS2-204       nonsense_media…       1.54
#> 6 ENST00000466691.1 ENSG00… PTGS2     PTGS2-202       retained_intron       1.01
#> # ℹ 8 more variables: assessed <lgl>, mean_inf_rv <dbl>, max_inf_rv <dbl>,
#> #   mean_inf_cv <dbl>, uncertain <lgl>, min_sibling_cor <dbl>,
#> #   n_ambiguous_partners <int>, ambiguity_group <chr>
isoform_ambiguity(iso, ambiguous_only = TRUE)
#> # A tibble: 2 × 10
#>   gene_id   gene_name transcript_id_1 transcript_id_2 mean_cor min_cor n_samples
#>   <chr>     <chr>     <chr>           <chr>              <dbl>   <dbl>     <int>
#> 1 ENSG0000… PTGS2     ENST0000036746… ENST0000068160…   -0.618  -0.873         6
#> 2 ENSG0000… NFKBIA    ENST0000069795… ENST0000069796…   -0.528  -0.787         5
#> # ℹ 3 more variables: ambiguous <lgl>, transcript_name_1 <chr>,
#> #   transcript_name_2 <chr>
head(de_confidence(iso))
#> # A tibble: 6 × 23
#>   feature_id       feature_type gene_id gene_name transcript_ids transcript_name
#>   <chr>            <chr>        <chr>   <chr>     <chr>          <chr>          
#> 1 ENST00000263341… transcript   ENSG00… IL1B      ENST000002633… IL1B-201       
#> 2 ENST00000307407… transcript   ENSG00… CXCL8     ENST000003074… CXCL8-201      
#> 3 ENST00000697966… transcript   ENSG00… NFKBIA    ENST000006979… NFKBIA-221     
#> 4 ENST00000496280… transcript   ENSG00… IL1B      ENST000004962… IL1B-208       
#> 5 ENST00000651122… transcript   ENSG00… CCL5      ENST000006511… CCL5-203       
#> 6 ENST00000415035… transcript   ENSG00… STAT1     ENST000004150… STAT1-205      
#> # ℹ 17 more variables: transcript_type <chr>, mean_count <dbl>,
#> #   lfc_median <dbl>, lfc_lo <dbl>, lfc_hi <dbl>, lfc_sd <dbl>,
#> #   sign_agreement <dbl>, frac_draws_effect <dbl>, frac_draws_sig <dbl>,
#> #   mean_inf_rv <dbl>, ambiguous <lgl>, confidence <ord>, reason <chr>,
#> #   det_log2FC <dbl>, det_qvalue <dbl>, n_det_sig <int>, confident_de <lgl>
uncertainty_report(iso)$transcripts
#> # A tibble: 4 × 7
#>   transcript_type n_assessed n_uncertain pct_uncertain median_inf_rv n_ambiguous
#>   <chr>                <int>       <int>         <dbl>         <dbl>       <int>
#> 1 protein_coding          47          21          44.7         0.373           1
#> 2 retained_intron         27          12          44.4         0.324           2
#> 3 nonsense_media…         17          13          76.5         1.36            1
#> 4 protein_coding…          6           5          83.3         1.36            0
#> # ℹ 1 more variable: median_unique_fraction <dbl>
```
