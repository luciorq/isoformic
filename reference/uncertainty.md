# Quantify and report the uncertainty of isoform-level results

An isoform-level result is only as good as the quantification behind it.
Reads that map to several isoforms of a gene can be shared out between
them in more than one way, so the estimate of each isoform carries an
*inferential* uncertainty that is unrelated to biological variability.
Salmon measures it by drawing replicates from the posterior (Gibbs
sampling, `--numGibbsSamples`) and by writing the equivalence classes
(`--dumpEq`). These functions turn that information into columns that
describe each isoform and each call, to be reported next to the results.

Every isoform gets four kinds of evidence:

- **Quantification uncertainty**: the inferential relative variance
  (`mean_inf_rv`, Zhu et al. 2019, as in `fishpond`), the variance of
  the replicates in excess of what counting noise alone would give,
  relative to the expression. Isoforms with a high value are not well
  determined by the reads. `uncertain` flags
  `mean_inf_rv >= inf_rv_cutoff`.

- **Ambiguity between isoforms**: the correlation between the replicates
  of isoforms of the same gene. A strong negative correlation means the
  reads cannot tell the two isoforms apart and every read taken from one
  is given to the other. Pairs at or below `cor_cutoff` are *ambiguous*
  and are joined into `ambiguity_group`s that can be tested as one
  feature. This is the idea of `mmcollapse` (Turro et al. 2014), with a
  fixed cutoff.

- **Read support**: the share of the fragments of an isoform that map to
  it alone (`unique_fraction`, from the equivalence classes).

- **Stability of the fold-change**: the contrast between two groups is
  repeated on every replicate. `lfc_lo`/`lfc_hi` is the 95% interval of
  the log2 fold-change over replicates, `sign_agreement` the share of
  replicates that agree on its direction and `frac_draws_sig` the share
  in which the isoform is called differentially expressed.

`de_confidence()` repeats the contrast on every replicate. `unstable`
flags a fold-change whose direction fewer than `min_agreement` (95%) of
the replicates agree on, and `unstable_de` a call of the primary test
with such a direction. The other columns describe the quantification of
the isoform, and `flags` lists those that apply. None of them is a
verdict on a call: report them next to it.

- `add_uncertainty(x, paths)` reads the replicates, equivalence classes
  and run metadata from the Salmon output directories (or Oarfish output
  prefixes, see
  [`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md))
  `paths` (named by sample ID, in any order), computes the tables below
  and stores them in the object, so that they are written by
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md).
  When `group` is given it also stores the fold-change stability of that
  contrast.

- `tx_uncertainty(x)`: per-transcript table, for the transcripts whose
  draws are in the object (the expressed ones).

- `isoform_ambiguity(x)`: per-pair table of replicate correlations.

- `de_confidence(x, group)`: per-feature table of fold-change stability,
  recomputed from the replicates; with `collapse_groups = TRUE`,
  ambiguous isoforms are summed into one feature before testing. Without
  `group` it returns the table that `add_uncertainty(group = )` stored.

- `uncertainty_report(x)`: a list of summary tables for a methods or
  results section.

- `uncertainty_methods(x)`: a ready-to-edit paragraph that states the
  settings used.

The replicates themselves are kept in `x@metadata$inf_reps`, a list with
one transcript x sample matrix per replicate.
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
writes them, one Parquet file per replicate, with the tables derived
from them, and
[`read_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
restores them; with `keep_draws = FALSE` only the tables are kept and
written. `x[i, j]` keeps the replicates of the kept transcripts and
samples as a view: it shares the draws of `x` and records the selection,
so the draws of `x` stay in memory while the subset exists (see the
Subsetting section of
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md));
to contrast two groups of a larger object, pass `levels` instead of
subsetting. Draws are read only for the transcripts of `x`; with a
`counts` assay, only for those with a mean count of at least
`min_count`, and the others are treated as unexpressed.

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
  identifiability = FALSE,
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

  Salmon output directories or Oarfish output prefixes (see
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)),
  named by sample ID. `NULL` uses the replicates already stored in the
  object.

- group:

  Name of a
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  column with the sample groups. In `de_confidence()`, `NULL` returns
  the stored table.

- levels:

  The two group levels to contrast, reference first. Samples of other
  groups are left out of the contrast. `NULL` uses the levels of
  `group`, which must then be exactly two.

- n_draws:

  Number of replicates to read; `NULL` reads all of them.

- min_count:

  Minimum mean estimated count for a transcript to be assessed (and,
  with a `counts` assay, for its replicates to be read).

- inf_rv_cutoff:

  Inferential relative variance from which an isoform is called
  uncertain. In `de_confidence()`, `NULL` uses the value the tables were
  computed with; without `group`, a value recomputes the "uncertain
  quantification" flag of the stored table. It describes the isoform and
  is not used to judge calls.

- cor_cutoff:

  Mean replicate correlation at or below which two isoforms of a gene
  are called ambiguous: a number, or `"adaptive"` to derive it from the
  data as `mmcollapse` does (minus the 97.5th percentile of the largest
  mean correlation of each isoform, which treats the positive
  correlations as noise). The value used is stored and stated by
  `uncertainty_methods()`.

- min_agreement:

  Share of replicates that must agree on the direction of a fold-change;
  below it the fold-change is flagged `unstable`. In `de_confidence()`,
  `NULL` uses the stored value; without `group`, a value recomputes
  `unstable`, `flags` and `unstable_de` of the stored table from its
  `sign_agreement`, without the replicates.

- equivalence:

  Logical, read the equivalence classes when present, or the table
  returned by
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
  to attach read support when the draws are supplied without `paths`.

- identifiability:

  Logical, also compute which transcripts the equivalence classes can
  tell apart
  ([`equivalence_identifiability()`](https://luciorq.github.io/isoformic/reference/equivalence_identifiability.md);
  minutes per sample, needs `paths` and, for human samples, the `Matrix`
  package), or the table that function returned. `de_confidence()` then
  reports `identifiable`.

- run_qc:

  Logical, add the run metadata of Salmon or Oarfish (see
  [`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md))
  to
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  with the column names prefixed by `salmon_` or `oarfish_`.

- keep_draws:

  Logical, keep the replicates in the object (needed to recompute the
  tables, for `de_confidence(group = )`,
  [`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  and the uncertainty of
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)).
  Kept replicates are also written by
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md).

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
  `ambiguous` and `collapse_groups` then follow these groups; the pairs
  of `isoform_ambiguity()` keep describing the correlation.

- ...:

  Not used; reserved for methods.

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

- `tx_uncertainty()`: a tibble with `transcript_id`, `gene_id`,
  `gene_name`, `transcript_name`, `transcript_type`, `mean_count`,
  `assessed` (`mean_count >= min_count`; the statistics are `NA` for the
  others), `mean_inf_rv`, `max_inf_rv`, `mean_inf_cv`, `uncertain`,
  `min_sibling_cor`, `n_ambiguous_partners`, `ambiguity_group` and, when
  equivalence classes were read, `n_classes`, `n_fragments`,
  `n_unique_fragments`, `unique_fraction` and `low_information` (`TRUE`
  for an uncertain isoform without a single fragment of its own, the
  isoforms `mmcollapse` selects for collapsing). With `identifiability`,
  also `identifiable_fraction` and `identifiable`.

- `isoform_ambiguity()`: a tibble of isoform pairs of the same gene,
  sorted by `mean_cor`, with `gene_id`, `gene_name`, `transcript_id_1`,
  `transcript_id_2`, `mean_cor`, `min_cor`, `n_samples`,
  `collapse_gain`, `ambiguous`, `transcript_name_1` and
  `transcript_name_2`. `collapse_gain` is the score Terminus groups on
  (Sarkar et al. 2020): the inferential relative variance of the sum of
  the two isoforms minus the mean of their own, averaged over samples.
  The more negative, the more is gained by testing the pair as one
  feature.

- `de_confidence()`: a tibble with one row per feature, sorted by the
  absolute median fold-change: `feature_id`, `feature_type`
  (`"transcript"` or `"group"`), `gene_id`, `gene_name`,
  `transcript_ids` (the members, separated by `;`), `transcript_name`,
  `transcript_type`, `mean_count`, the fold-change over the replicates
  (`lfc_median`, `lfc_lo`, `lfc_hi`, `lfc_sd`), `sign_agreement`,
  `frac_draws_effect` (share of replicates that agree on the direction
  with an absolute fold-change of at least the cutoff of
  [`de_cutoffs()`](https://luciorq.github.io/isoformic/reference/dea.md)),
  `frac_draws_sig`, `mean_inf_rv` (`NA` for a feature that is not
  assessed, which is then never flagged as uncertain), `ambiguous`,
  `identifiable` (only when the identifiability was computed; `NA` for
  groups), `flags` (the descriptive flags that apply, separated by
  `"; "`: "unstable direction", "uncertain quantification", "ambiguous
  with a sibling isoform", "not identifiable from the equivalence
  classes"; `NA` when none applies), the current differential expression
  result (`det_log2FC`, `det_qvalue`, `n_det_sig`, the number of
  significant members; read from
  [`de_tx()`](https://luciorq.github.io/isoformic/reference/dea.md) each
  time, so they follow `de_tx<-` and
  [`de_cutoffs()`](https://luciorq.github.io/isoformic/reference/dea.md))
  and `unstable_de`, `TRUE` for a feature with a significant member and
  an unstable direction (`NA` without stored results). `unstable` (after
  `sign_agreement`) is `sign_agreement < min_agreement`, `NA` for a
  feature that is not assessed. The sign agreement is measured against
  the median sign over the replicates, not against `det_log2FC`.

- `uncertainty_report()`: a list with `settings` (a list of the
  cutoffs), and the tibbles `replicates` (`replicate_type`, `n_draws`,
  `resampling_df`, the degrees of freedom over the study, and
  `median_overdispersion`, `NA` unless the draws are in the object and
  close to 1 for long reads), `run` (the run metadata in
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)),
  `transcripts` (by transcript type, with `median_inf_rv`, `q25_inf_rv`
  and `q75_inf_rv`), `genes`, `de` (significant features per flag; the
  flags describe the calls and are not a verdict: one row per flag,
  `unstable` first, with `flag`, `rule`, `n_significant`, `n_assessed`,
  `n_flagged` and `pct_flagged`; a collapsed group with a significant
  member counts once, and is never `ambiguous`, since its members were
  summed) and `switches` (genes with isoform switches, all and without
  unstable directions: `n_genes`, `n_opposite_direction` and
  `n_det_not_deg`). `run` and `de` are `NULL` when their input is
  missing; `switches` needs only the transcript-level table
  ([`de_tx()`](https://luciorq.github.io/isoformic/reference/dea.md))
  and is `NULL` without it, its "without unstable directions" row is
  `NA` without the fold-change stability of `add_uncertainty(group = )`,
  and its `n_det_not_deg` is `NA` without a gene-level table
  ([`de_gene()`](https://luciorq.github.io/isoformic/reference/dea.md)).

- `uncertainty_methods()`: a character string.

## Details

What the evidence supports. In a genome-wide evaluation on a placenta
cohort (GENCODE 50, 30 Gibbs draws; see the article "Is quantification
uncertainty worth checking?"), calls with an unstable direction
replicated less often in held-out halves of the cohort than the other
calls (35% against 64% pooled; odds ratio 0.36, 95% interval 0.31 to
0.41 with genes resampled as clusters, adjusted for count, test
statistic, number of isoforms of the gene and transcript type; the
interval is for that cohort and does not include the choice of samples).
Their raw replication rate was lower in all 60 directed comparisons (30
half-splits of the same samples, each used in both directions); against
flags of the same size made of the lowest-count or lowest-\|z\| calls,
the gap was larger in all 60. The flag marked 6% of the calls and only
about 10% of the calls that failed to replicate, so a stable direction
is not evidence that a call is right: in comparisons of control samples,
93% of the false calls had a stable direction. The flag was a secondary
measure of that evaluation, chosen as the headline after the results;
its threshold of 0.95 was fixed in advance and no other was tried, and
other primary tests (DESeq2, swish) and long-read bootstraps were not
evaluated. Ambiguity predicted replication little (odds ratio 0.87);
isoforms not identifiable from the equivalence classes, or in the top
third of InfRV, replicated somewhat less (0.56 and 0.63). An InfRV of at
least 0.5 marked 87% of the expressed isoforms.

The thresholds are heuristics, not universal constants.
`min_agreement = 0.95` is the threshold of that evaluation. It is a
share of the replicates: with fewer than 20 replicates no replicate may
disagree, with 30 one may, and with 100 five may. `inf_rv_cutoff = 0.5`
marks inferential relative variance above the sampling noise of about 30
draws: with 30 draws, sampling noise alone puts the value of a
well-determined, highly expressed transcript near 0.1, while poorly
determined isoforms typically reach several units. It describes the
isoform and is not used to judge calls. `cor_cutoff` (-0.5) defines the
ambiguity groups and the features of `collapse_groups`. Rerun with other
cutoffs if the data call for it. The replicates measure the uncertainty
of the quantification only; they do not replace biological replicates
and cannot make a noisy contrast significant.

`lfc_lo`-`lfc_hi` is the spread of the screening statistic over the
replicates, not a confidence interval for the point estimate of the
primary test; for low counts that estimate can lie outside it. The sign
agreement is measured against the sign of the median fold-change over
the replicates (`lfc_median`), not against `det_log2FC`.

What is not measured. The replicates describe how much an estimate moves
under the quantification model. They do not find the transcripts whose
abundance is not identifiable, that is, where several abundance vectors
explain the reads equally well (Zheng et al. 2022 report that only about
half of the transcripts flagged by replicate variance are also
non-identifiable in that sense), and they do not find an estimate that
is wrong because the annotation lacks the transcript the reads came from
(Ma and Kingsford 2019). An ambiguous pair is a statement about the
replicates, not a proof of non-identifiability.

How many replicates. The inferential relative variance and the
correlations are computed per sample, from the replicates of that
sample. The overdispersion used to scale counts (see
[`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md))
pools all samples, and its precision depends on the total
`n_samples * (n_replicates - 1)`: Baldoni, Chen and Smyth (2024)
recommend at least 60 and prefer 200 to 300, and find that Gibbs samples
resolve low counts better than bootstrap replicates. `add_uncertainty()`
says so when the total is below 60.

Fold-changes are computed on counts scaled to a common library size, as
the difference between the groups of the mean of
`log2(scaled count + pseudocount)`, with a pseudocount of 5 (the `swish`
convention): a log ratio of geometric means, on the scale of the
per-replicate test, a Welch t-test on that log scale with
Benjamini-Hochberg adjustment across features. It is not the log ratio
of the arithmetic group means, and the two can differ in sign when one
sample dominates a group. This is a screening statistic for stability,
not a replacement for the primary test (`swish`, `DESeq2`, ...), whose
results are shown next to it. When the sample table (`col_data`) has a
`library_size` column, with the total estimated count of the full
quantification, it is used for the scaling. Set it when the object holds
only a subset of the transcripts, so that the subset does not change the
scaling.

## References

Zhu, A., Srivastava, A., Ibrahim, J.G., Patro, R., Love, M.I. (2019).
Nonparametric expression analysis using inferential replicate counts.
*Nucleic Acids Research* 47(18), e105.

Sarkar, H., Srivastava, A., Bravo, H.C., Love, M.I., Patro, R. (2020).
Terminus enables the discovery of data-driven, robust transcript groups
from RNA-seq data. *Bioinformatics* 36(Suppl 1), i102-i110.

Turro, E., Astle, W.J., Tavare, S. (2014). Flexible analysis of RNA-seq
data using mixed effects models. *Bioinformatics* 30(2), 180-188.

Zheng, H., Ma, C., Kingsford, C. (2022). Deriving ranges of optimal
estimated transcript expression due to nonidentifiability. *Journal of
Computational Biology* 29(2), 121-139.

Ma, C., Kingsford, C. (2019). Detecting, categorizing, and correcting
coverage anomalies of RNA-seq quantification. *Cell Systems* 9(6),
589-599.

Baldoni, P.L., Chen, L., Smyth, G.K. (2024). Faster and more accurate
assessment of differential transcript expression with Gibbs sampling and
edgeR v4. *NAR Genomics and Bioinformatics* 6(4), lqae151.

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
  condition = rep(c("preeclampsia", "control"), each = 3)
)
iso <- read_salmon(quant_dirs, annot_path = annot, col_data = col_data)

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
#> 1 ENST00000367468.… ENSG00… PTGS2     PTGS2-201       protein_coding        9.93
#> 2 ENST00000681605.1 ENSG00… PTGS2     PTGS2-206       nonsense_media…       5.11
#> 3 ENST00000680451.1 ENSG00… PTGS2     PTGS2-205       protein_coding        9.70
#> 4 ENST00000490885.6 ENSG00… PTGS2     PTGS2-203       retained_intron       8.98
#> 5 ENST00000559627.1 ENSG00… PTGS2     PTGS2-204       nonsense_media…       6.61
#> 6 ENST00000466691.1 ENSG00… PTGS2     PTGS2-202       retained_intron       1.58
#> # ℹ 8 more variables: assessed <lgl>, mean_inf_rv <dbl>, max_inf_rv <dbl>,
#> #   mean_inf_cv <dbl>, uncertain <lgl>, min_sibling_cor <dbl>,
#> #   n_ambiguous_partners <int>, ambiguity_group <chr>
isoform_ambiguity(iso, ambiguous_only = TRUE)
#> # A tibble: 6 × 11
#>   gene_id   gene_name transcript_id_1 transcript_id_2 mean_cor min_cor n_samples
#>   <chr>     <chr>     <chr>           <chr>              <dbl>   <dbl>     <int>
#> 1 ENSG0000… NFKBIA    ENST0000113342… ENST0000021679…   -0.879  -0.937         6
#> 2 ENSG0000… CCL5      ENST0000060514… ENST0000100676…   -0.878  -0.983         6
#> 3 ENSG0000… CCL2      ENST0000113714… ENST0000022583…   -0.801  -0.962         6
#> 4 ENSG0000… NFKBIA    ENST0000069795… ENST0000069796…   -0.707  -0.787         6
#> 5 ENSG0000… CXCL10    ENST0000030660… ENST0000113500…   -0.707  -0.932         6
#> 6 ENSG0000… ICAM1     ENST0000113664… ENST0000026483…   -0.624  -0.704         6
#> # ℹ 4 more variables: collapse_gain <dbl>, ambiguous <lgl>,
#> #   transcript_name_1 <chr>, transcript_name_2 <chr>
head(de_confidence(iso))
#> # A tibble: 6 × 23
#>   feature_id       feature_type gene_id gene_name transcript_ids transcript_name
#>   <chr>            <chr>        <chr>   <chr>     <chr>          <chr>          
#> 1 ENST00000907726… transcript   ENSG00… SOCS3     ENST000009077… SOCS3-203      
#> 2 ENST00000540176… transcript   ENSG00… STAT1     ENST000005401… STAT1-212      
#> 3 ENST00000624362… transcript   ENSG00… CCL2      ENST000006243… CCL2-203       
#> 4 ENST00000580907… transcript   ENSG00… CCL2      ENST000005809… CCL2-202       
#> 5 ENST00000697958… transcript   ENSG00… NFKBIA    ENST000006979… NFKBIA-215     
#> 6 ENST00000225831… transcript   ENSG00… CCL2      ENST000002258… CCL2-201       
#> # ℹ 17 more variables: transcript_type <chr>, mean_count <dbl>,
#> #   lfc_median <dbl>, lfc_lo <dbl>, lfc_hi <dbl>, lfc_sd <dbl>,
#> #   sign_agreement <dbl>, unstable <lgl>, frac_draws_effect <dbl>,
#> #   frac_draws_sig <dbl>, mean_inf_rv <dbl>, ambiguous <lgl>, flags <chr>,
#> #   det_log2FC <dbl>, det_qvalue <dbl>, n_det_sig <int>, unstable_de <lgl>
uncertainty_report(iso)$transcripts
#> # A tibble: 4 × 9
#>   transcript_type  n_assessed n_uncertain pct_uncertain median_inf_rv q25_inf_rv
#>   <chr>                 <int>       <int>         <dbl>         <dbl>      <dbl>
#> 1 protein_coding           86          66          76.7          2.96      0.643
#> 2 nonsense_mediat…         28          27          96.4          3.22      1.94 
#> 3 retained_intron          27          16          59.3          1.27      0.154
#> 4 protein_coding_…          7           5          71.4          2.67      1.08 
#> # ℹ 3 more variables: q75_inf_rv <dbl>, n_ambiguous <int>,
#> #   median_unique_fraction <dbl>
```
