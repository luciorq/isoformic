# Quantification uncertainty: literature notes

Written 2026-11-11. Sources: Consensus, scite and Undermind searches run that
day. Purpose: ground the definitions in `R/uncertainty.R`,
`R/usage.R` and the uncertainty article, and list what the package does not
cover.

**Where to start.** Section 1 is the conclusion. Section 4 maps each
definition to the function that implements it. Section 5 is the backlog:
ten items (U1 to U10) with the function to change, the effect, the tests and
the effort, and its status (5.1). Section 2 holds
the formulas and numbers behind them. Section 8 is the pre-registered
genome-wide evaluation of 2026-11-18 (does uncertainty change the results?),
and 8.5 the package changes it proposed (U11 to U13), applied on 2026-11-20.

Evidence level per paper is marked: `[full]` the full text was read (through
Undermind PDF readers or scite excerpts), `[abs]` abstract only. Statements
from `[abs]` papers are what the abstract says, nothing more. All DOIs in the
reference list were checked with scite `verify_papers`: all resolved, no
retraction, correction or concern notices.

The Undermind deep search (91 ranked papers) is in the workspace
"Transcript quantification uncertainty: definitions and methods":
<https://app.undermind.ai/projects/84b80218-7f49-4444-bc06-a5aafe355fba?path=/Definitions%20of%20transcript%20quantification%20uncertainty>

## 1. Main conclusion

"Uncertainty" of a transcript abundance has two meanings in the literature,
and they are not equivalent:

1. **Inferential variability given the model.** How much the estimate moves
   when the reads are resampled (bootstrap) or the posterior is sampled
   (Gibbs). This is what Salmon draws measure, and what InfRV, the edgeR
   overdispersion, sleuth and Swish use.
2. **Non-identifiability.** Several abundance vectors have the same
   likelihood, so the optimum is a set, not a point. This is a property of
   the likelihood (annotation topology plus read design), not of sampling
   noise. Kingsford group.

A third, separate thing is **model misfit**: the annotation is incomplete or
the bias model is wrong, so the estimate is confidently wrong. Neither (1)
nor (2) under a complete reference detects it.

isoformic implements (1). It has no measure of (2) or of misfit, and since
2026-11-11 its documentation says so. See section 4. Since 2026-11-12
`equivalence_identifiability()` measures the structural part of (2) from
equivalence class membership (section 5.4); ranges of optima and misfit
are still not measured.

## 2. Definitions found

### 2.1 Read-to-transcript ambiguity (RTA) overdispersion (Smyth group)

Baldoni et al. 2024, NAR 52:e13 (online 2023) `[full]`, preprint title used
"mapping ambiguity overdispersion".

- Terminology. The paper lists the synonyms in use: "assignment ambiguity,
  mapping ambiguity, quantification uncertainty or inferential uncertainty",
  and picks "read to transcript ambiguity (RTA)" to avoid confusion with
  other sources of uncertainty.
- Model. The Poisson model for sequencing variability is replaced by a
  quasi-Poisson model fitted to the resampled counts of each transcript.
  Bootstrap resampling acts on equivalence classes, where reads map
  unambiguously, so without ambiguity the variation would be Poisson.
- Estimator, transcript `t`, samples `i = 1..n`, resamples `b = 1..B`,
  resampled count `u_tib`, `lambda_ti = mean_b(u_tib)`:

  ```
  sigma2_t  = (1 / d_t) * sum_i sum_b (u_tib - lambda_ti)^2 / lambda_ti
  d_t       = n * (B - 1)
  sigma2_0  = max(1, median_t(sigma2_t) / median(F(d_med, d_0)))   d_0 = 3
  sigma2~_t = max(1, (d_0 * sigma2_0 + d_t * sigma2_t) / (d_0 + d_t))
  ```

  `d_med` is the median `d_t` over expressed transcripts. The shrinkage is
  minor because `d_t` is usually 2 to 3 orders of magnitude above `d_0`
  (example in the paper: `d_t = 990` for 10 samples and 100 bootstraps), so
  the estimates are treated as known constants.
- Use. Scaled counts `z_ti = y_ti / sigma2~_t` go into a gene-level pipeline
  (edgeR, limma-voom). No rounding needed.
- The paper gives no rule for terms with `lambda_ti = 0` (zero denominator).
- Values. Only 11% of expressed transcripts of single-transcript genes have
  overdispersion above 1.11, against 90% for multi-transcript genes. Mean
  6.75 for genes with at least 10 annotated transcripts. Some values reach
  about 10^2.5. Predictors: transcript overlap (number of transcripts per
  gene), shorter reads, single-end reads (50 bp single-end: about 60% higher
  mean than 150 bp paired-end). Expression level is not named as a
  predictor. ONT long reads: estimates mostly close to 1.
- The paper does not relate its overdispersion to InfRV or to sleuth's
  inferential variance numerically. It describes sleuth as a decomposition
  into biological and inferential variance and Swish as a different testing
  approach.

Baldoni, Chen and Smyth 2024, NAR Genomics and Bioinformatics 6:lqae151
`[full]`.

- Replaces bootstrap with Gibbs resampling and edgeR v3 with v4. Describes
  the estimate as "the squared coefficient of variation of the resampled
  transcript counts, bounded below by 1, and with very light empirical Bayes
  moderation between transcripts".
- Why Gibbs is better: "The Bayesian framework gives non-zero probability to
  unobserved ECs and is able to obtain higher resolution resampling counts
  that are not limited to single EC scenarios or to granular resampling of
  reads in low count situations." Example: transcript HABP2-201 has
  bootstrap RTA-dispersion 1.00 and Gibbs 1.42; BBC3-202 is filtered out
  after bootstrap scaling and detected as DE with Gibbs.
- Number of resamples. Accuracy depends on the total degrees of freedom
  `N * (B - 1)` over the whole study (`N` = all samples, not per group):
  at least 60; 200 to 300 preferable for Gibbs; `B = 2` is enough when
  `N >= 200`. Bootstrap is already good at about 50.
- Speed: about 50% less Salmon time with Gibbs than bootstrap for 50 million
  read pairs and 30 resamples.
- What it does not capture: the model separates RTA dispersion (annotation
  topology) from other technical variation (for example PCR duplication) and
  from biological variation; those are handled by the count model.
- No mention of Salmon thinning. Salmon run with default options otherwise.

Baldoni et al. 2025, NAR, doi 10.1093/nar/gkaf1305 `[abs]`: the same divided
counts feed `diffSplice` (limma and edgeR) for differential transcript usage;
claims better power and FDR control than specialised DTU methods, especially
with few samples.

Chen et al. 2025, edgeR v4, NAR `[abs]`: fractional counts, new
quasi-likelihood treatment that is more accurate for small counts, and DTE /
DTU functionality.

### 2.2 Inferential relative variance, InfRV (Love and Patro groups)

Zhu et al. 2019, Swish, NAR 47:e105 `[abs]`: nonparametric test over
inferential replicates (extends SAMseq); better FDR control for transcripts
with high inferential uncertainty.

Sarkar et al. 2020, Terminus, Bioinformatics 36:i102 `[full]`:

```
InfRV(g_i) = max(0, var(g_i) - mean(g_i)) / (mean(g_i) + p) + s
p = 5 (pseudocount), s = 0.01 (shift)
score(i, j) = InfRV(g_i + g_j) - (InfRV(g_i) + InfRV(g_j)) / 2
```

- A pair is a collapse candidate only if (a) the two transcripts co-occur in
  an equivalence class, (b) at least one passes the screen
  `mean(g) >= 1 and (max(g) - min(g)) / mean(g) > 0.1`, and (c) the score is
  at or below a data-derived threshold: the 2.5th percentile of the scores of
  random expressed pairs.
- Before that, transcripts with the same equivalence-class labels and
  near-identical conditional probabilities are collapsed (the data cannot
  distinguish them).
- Greedy: collapse the edge with the lowest score, update, repeat.
- Consensus over samples: union graph of per-sample groups, edges weighted by
  the number of samples, edges below the consensus threshold dropped
  (default `m / 4` samples), connected components are the groups.
- Terminus does **not** use a correlation cutoff. It credits the
  anti-correlation idea to mmcollapse.
- 100 Gibbs samples in all experiments.

`fishpond::computeInfRV()` (version 2.16.0, checked in the installed
package):

```r
InfRV <- pmax(infVar - infMean, 0) / (infMean + pc) + shift   # pc = 5, shift = 0.01
```

Singh, Love and Patro 2023, TreeTerminus, iScience `[abs]`: transcripts in a
tree where inferential uncertainty decreases, on average, towards the root.
Singh et al. 2025, mehenDi, Genome Research `[abs]`: differential testing on
that tree, selected nodes chosen from the data, no two in an
ancestor/descendant relation.

### 2.3 Anti-correlation of posterior samples (mmcollapse)

Turro, Astle and Tavare 2014, Bioinformatics 30:180 `[full]`.

- For a set `S`, aggregate expression in sample `i` is
  `phi_iS = sum_{t in S} mu_it`.
- The correlation is between the aggregates of two sets over the joint MMSEQ
  posterior samples, computed per biological sample, then averaged over
  samples. It is not a correlation of point estimates across samples.
- Iteration: merge the two sets with the lowest mean posterior correlation;
  stop when the minimum mean correlation exceeds a threshold.
- Threshold is data-adaptive: "We use the right tail of the empirical
  distribution of maximum mean correlations as a control and set the stopping
  threshold to minus the 97.5th percentile by default." No universal numeric
  cutoff.
- Initialisation: identical sequences are collapsed first; the procedure then
  runs on "low-information transcripts" (no uniquely mapping reads and
  posterior SD above the maximum for transcripts with at least one unique
  read). Transcripts with zero alignments in all samples are excluded.
- Rationale: the reads may identify the total of a set while not
  distinguishing its members, so the set is the unit of inference.

### 2.4 Variance decomposition and measurement error

- Pimentel et al. 2017, sleuth, Nature Methods `[abs]`: response-error linear
  model that separates biological variance from inferential variance
  estimated from kallisto bootstraps.
- Young, Van Buren and Rashid 2023, CompDTU / CompDTUme, Biostatistics
  `[abs]`: compositional regression for DTU; the `me` variant adds
  quantification uncertainty as measurement error and helps for genes with
  high uncertainty when the sample size is sufficient.
- From the Undermind list, not read: BANDITS (Tiberi and Robinson 2019), RATs
  (Froussios et al. 2019), SEESAW (Wu et al. 2023), BitSeq (Glaus et al.
  2012), RSEM credibility intervals (Li and Dewey 2011), MMSEQ (Turro et al.
  2011), Polee (Jones and Ruzzo 2021), IsoEM2 bootstrap intervals (Mandric et
  al. 2017).

### 2.5 Non-identifiability (Kingsford group)

Zheng, Ma and Kingsford 2022, J Comput Biol 29:121 (bioRxiv 2019) `[full]`.

- **Range of optima** of transcript `T_j`: minimum and maximum abundance over
  all estimates that are equally optimal under the model. Every value in
  between is attained by some optimal estimate.
- Complete reference: keep the totals of the sufficient statistics fixed
  (sum of abundances of the transcripts containing each junction; for
  fragments spanning several junctions, each "phasing path"), keep abundances
  non-negative, and solve two linear programs per transcript (minimise,
  maximise). A transcript is non-identifiable when the two differ.
- Incomplete reference (graph quantification): find the optimal edge flow on
  the splice graph; every decomposition of the flow into paths is equally
  optimal; bound the weight of the target path with max-flow (AND-Quant).
- Hybrid: `[l^k, u^k] = [k*l^1 + (1-k)*l^0, k*u^1 + (1-k)*u^0]`, `k` =
  assumed share of expression from unannotated paths (reference completeness
  `1 - k`; shown as lambda in the figures).
- Human Body Map, 16 samples: 35 to 50% of transcripts have a non-singleton
  range; for 20 to 47% the range is wide enough that the ranking against at
  least one sibling isoform cannot be decided.
- Relation to draws. The authors treat bootstrap variance as uncertainty from
  finite sample size, a different source: "Only half of transcripts for which
  the uncertainty estimates are caused by finite sample size are common to
  the transcripts under the model nonidentifiability case" (prostate sample,
  k = 70%). And: "the difference between Terminus and our approach for
  identifying uncertain expression estimates lies in the source of error and
  whether incomplete reference is considered."
- DE reliability flag: overlap = size of the intersection of the two group
  mean ranges divided by the size of the smaller range; flagged when above
  25%. Result: 5 of 257 DE transcripts (MCF10) and 19 of 3,152 (CD8 T cells)
  flagged. So most DE calls survive.
- Code: <https://github.com/Kingsford-Group/subgraphquant>. The paper does
  not say that plain Salmon output is sufficient input.

Ma, Zheng and Kingsford 2021, Algorithms Mol Biol 16:5 `[abs]`: graph
quantification with variable fragment lengths, where every full path of the
splice graph is a possible transcript.

Related, from the Undermind list, not read: Hiller et al. 2009
(identifiability of isoform deconvolution); Ferrer-Bonsoms et al. 2022
(identifiability and fragment length choice; DOI not verified, do not cite
without checking); Zhang et al. 2026, Brief Bioinform, doi 10.1093/bib/bbag134
`[abs]` (condition number of a gene-specific matrix as a bound on
quantification error; decreases with read length).

### 2.6 Model misfit and evaluation without truth

Ma and Kingsford 2019, SAD, Cell Systems 9:589 `[full, bioRxiv version]`.

- Anomaly: observed fragment-start coverage of a transcript against Salmon's
  expected coverage (with its sequence, GC and positional bias models).
  Over-expression score of a region `[a, b]` is
  `max(sum_{i=a..b}(obs[i] - exp[i]), 0)`, under-expression the reverse; the
  transcript score is the maximum over contiguous regions, tested for
  significance.
- Reads are reassigned among isoforms by linear programming. If the anomaly
  disappears it is "adjustable" (read assignment mistake); if not,
  "unadjustable" (likely a missing transcript).
- 103 genes have unadjustable anomalies common to all 30 GEUVADIS and 16
  Human Body Map samples. 50 to 60% of unadjustable anomalies have a matching
  isoform assembled by StringTie or Scallop. Adjusted quantification removes
  2.4 to 6% of false positive DE calls. No overall fraction of transcripts
  affected is reported.
- On draws: "statistical techniques such as bootstrapping and Gibbs sampling
  can associate confidence intervals to expression estimates, but these
  techniques provide little insight into the causes of low confidence or
  misquantification, and detect a subset of misquantifications."

Fan, Chan and Patro 2022, Perplexity, Algorithms Mol Biol 17:6 `[full]`:
`PP(F_heldout, theta) = exp(-(1/N) * log P(F_heldout | theta))`, computed
over range-factorized equivalence classes, with smoothing for fragments of
zero likelihood. Lower is better. Point for us: uncertainty of an estimate is
not an evaluation of its accuracy.

From the Undermind list, not read: Soneson et al. 2018, junction coverage
compatibility (JCC) score; Papastamoulis et al. 2014, variational Bayes
underestimates posterior variance.

### 2.7 Long reads

- RTA overdispersion is close to 1 for ONT reads (Baldoni et al. 2024).
- Wissel et al. 2026, Genome Biology `[abs + excerpt]`: PacBio Kinnex against
  Illumina. Illumina detects more transcripts but many are unstable or
  ambiguous in complex genes; after correcting Illumina for inferential
  variability the two platforms are highly concordant. Oarfish nominated as
  the most efficient tool for Kinnex.
- Jousheghani, Singh and Patro 2025, Oarfish, Bioinformatics `[abs]`:
  coverage score in the fragment assignment probability; error metrics nearly
  halved with it, larger gain for lower-fidelity data.
- Undermind summary (abstract level): long-read tools are less explicit about
  per-isoform uncertainty; DELongSeq (Hu et al. 2023) uses EM information
  matrices; miniQuant (Li et al. 2025) not read.

## 3. Terms to use

| Term | Meaning | Source |
|---|---|---|
| RTA overdispersion | Variance inflation of a transcript count over Poisson, from draws, pooled over samples, at least 1 | Baldoni 2024 |
| InfRV | Excess of draw variance over the mean, relative to the mean, per sample, averaged | Zhu 2019, Sarkar 2020 |
| Posterior anti-correlation | Negative correlation of draws of two transcripts (or sets) within a sample | Turro 2014 |
| Range of optima | Interval of abundance over equally optimal solutions | Zheng 2022 |
| Non-identifiability | The range of optima is not a single point | Zheng 2022 |
| Coverage anomaly | Observed coverage contradicts the model expectation | Ma and Kingsford 2019 |

"Identifiability" must be kept for the fourth and fifth rows.

## 4. Comparison with isoformic

Checked against the code on 2026-11-11, before the fixes of section 5 (U1,
U2, U3, U9 are now applied; the rows below describe what was found).
Function names are given instead of line numbers.

| Item | Where in isoformic | Literature | Verdict |
|---|---|---|---|
| Overdispersion | `overdispersion_from_moments()` through `object_overdispersion()` (`R/uncertainty.R`), used by `inferential_overdispersion()`, `scale_counts()`, `overdispersion_for()` and `test_usage()` | Baldoni 2024 estimator (2.1) | Same: Pearson statistic over draws, `df = n_expressed * (B - 1)`, prior `max(1, median(phi) / qf(0.5, median(df), 3))`, moderated, floored at 1. Terms with zero mean are left out of sum and df, a case the paper does not specify. The moments come from every transcript of the index (`moments` attribute of `read_inferential_replicates()`), so the prior is the one of `catchSalmon()` whatever `min_count` keeps (checked on GENCODE 50 in the evaluation of section 8: 654,828 placenta targets, largest relative difference 1.35e-11; the earlier GENCODE 49 check gave 4e-12 on 517,038 transcripts). Replicates without moments fall back to the stored transcripts |
| Number of draws | `add_uncertainty(n_draws = )`; 30 Gibbs draws in examples and articles; stored in `metadata$uncertainty_n_draws` | `N * (B - 1) >= 60`, better 200 to 300 | Met from 3 samples (87); 7 samples reach 203. Not checked or reported by the package |
| Replicate type | `read_meta_info()` reads `samp_type`; stored in `metadata$uncertainty_replicate_type` | Gibbs preferred over bootstrap (2.1) | Recorded, no advice given |
| InfRV | `inferential_rv()`: `pmax(0.01, var - mean) / (mean + 5)` | `pmax(var - mean, 0) / (mean + 5) + 0.01` (fishpond, Terminus) | **Differs.** Ours is 0.01 lower when `var - mean > 0.01`, and equals `0.01 / (mean + 5)` instead of `0.01` when `var <= mean`. The roxygen says "as in fishpond" |
| Ambiguity rule | `inferential_correlation()` and `assign_ambiguity_groups()`: mean per-sample correlation of draws between isoforms of a gene; ambiguous at or below `cor_cutoff = -0.5`; connected components | mmcollapse (2.3): same statistic on sets, greedy, data-adaptive threshold | Same idea, fixed cutoff, pairs only (no re-evaluation of merged sets). Turro et al. 2014 is not cited |
| Wording | `?add_uncertainty` calls the correlation rule "Identifiability" | Identifiability is a property of the likelihood (2.5) | Loose. Zheng 2022 reports about half overlap between the two notions |
| Read support | `equivalence_support()`: `unique_fraction` | mmcollapse selects "low-information transcripts" by "no uniquely mapping reads" | Consistent precedent; can be cited |
| Equivalence classes | `read_eq_classes()` keeps members and counts, **drops the weights** | Ranges of optima and Terminus both need the conditional probabilities | The membership part of U5 is done without them; ranges and a weighted variant need them kept |
| Fold-change stability | `compute_de_confidence()`: Welch t-test per draw, `sign_agreement`, `frac_draws_sig` | Swish tests over draws; Zheng 2022 flags DE by range overlap above 25% | Ours is a screening statistic, as documented |
| Usage tests | `usage_diffsplice()` on counts divided by `overdispersion_for()` | Baldoni 2025 | Same pipeline, already cited in `R/usage.R` |
| Groups and trees | `read_terminus_groups()`, `run_terminus()`, `read_treeterminus()`, `cut_transcript_trees()` | Sarkar 2020, Singh 2023 | Covered. mehenDi (Singh 2025) testing on the tree is not |
| Long reads | `read_oarfish()`, Parquet replicates | RTA overdispersion near 1 for ONT | Covered; the expectation "near 1" is not stated in the docs |
| Non-identifiability | none | Ranges of optima (2.5) | `equivalence_identifiability()`: structural, from class membership (U5, section 5.4); ranges not implemented |
| Misfit | none | SAD, JCC, perplexity (2.6) | Not implemented, see section 5.5 |

## 5. Backlog: what is relevant and implementable

Status on 2026-11-12: U1, U2, U3 and U9 are applied (see 5.3), U4, U6 and
U8 are applied and U7 was tested and not adopted (see 5.4); the first step
of U5 (structural identifiability) is applied (see 5.4); U10 is open.
Effort: S = under half a day, M = one to two days, L = more. "Derivation" marks reasoning of these notes that is not in a
paper and must be validated before it is trusted.

### 5.1 Summary

| ID | Item | Effort | Needs new input | Changes numbers | Source | Status |
|---|---|---|---|---|---|---|
| U1 | InfRV in the fishpond form | S | no | yes, slightly | 2.2 | done |
| U2 | Rename "Identifiability", cite mmcollapse, state what is not measured | S | no | no | 2.3, 2.5 | done |
| U3 | Report and check `N * (B - 1)` and the replicate type | S | no | no | 2.1 | done |
| U4 | Terminus collapse score per pair in `isoform_ambiguity()` | S-M | no | adds a column | 2.2 | done |
| U5 | Range of optima from equivalence classes | L | `--dumpEqWeights` | adds columns | 2.5 | first step done; ranges open |
| U6 | Data-adaptive `cor_cutoff` | M | no | changes groups (opt-in) | 2.3 | done |
| U7 | DE reliability flag by interval overlap between groups | S-M | no | adds a column | 2.5 | not adopted |
| U8 | Low-information flag | S | equivalence classes | adds a column | 2.3 | done |
| U9 | State long-read expectation in docs and report | S | no | no | 2.1, 2.7 | done |
| U10 | Gene-level separability score from the annotation | L | read the paper first | new function | 2.5 | open |

Open items: U10, and the ranges of U5, which need a decision on the model
first (see 5.4).

U11 to U13 come from the genome-wide evaluation and are listed in 8.5
(applied on 2026-11-20: `unstable` reported as the flag of a call, the
confidence class removed, `inf_rv_cutoff = 0.5` kept as a description).

### 5.2 Items

**U1. InfRV in the fishpond form.**
- Where: `inferential_rv()` in `R/uncertainty.R`. One line:
  `pmax(moments$var - moments$mean, 0) / (moments$mean + pseudocount) + shift`.
- Effect: `mean_inf_rv` and `max_inf_rv` rise by about 0.01; transcripts with
  `var <= mean` get exactly 0.01 instead of a value near 0.
- Also touch: the sentence in `?add_uncertainty` details ("near 0.1 with 30
  draws", cutoff 0.5), the formula in `ARCHITECTURE.md` section
  "Quantification uncertainty" (it states the current form), the inlined
  copy of the same formula in `cut_transcript_trees()` (`R/treeterminus.R`,
  `rv[, s] <- pmax(defaults$shift, va - mu) / ...`; make it call
  `inferential_rv()`), and the numbers in
  `isoformic-uncertainty.qmd` and `isoformic-quantification-stability.qmd`
  (CLAUDE.md rule: check every number the articles state). Bundled tables
  under `inst/extdata/knockdown/` hold precomputed `mean_inf_rv`; rebuild
  with `data-raw/build_knockdown_example.R` or note the offset.
- Tests: `tests/testthat/test-uncertainty.R` uses inequalities (`< 0.5`,
  `> 1`, `> 3`) that should survive; add one equality test against
  `fishpond::computeInfRV()` on a small array (`skip_if_not_installed`).
- Alternative with no number change: keep the code, drop "as in fishpond"
  and state the formula.

**U2. Wording and citation.**
- Where: roxygen of `add_uncertainty` (bullet "Identifiability"),
  `uncertainty_methods()` text, `isoformic-uncertainty.qmd`,
  `vignettes/articles/library.bib` (add Turro 2014, Zheng 2022, Baldoni 2024
  NARGAB, Ma and Kingsford 2019).
- Text to carry: the correlation rule is the mmcollapse idea (Turro et al.
  2014) with a fixed cutoff; it shows that the draws of two isoforms trade
  reads; it is not the non-identifiability of Zheng et al. 2022, which the
  package does not measure; and none of the columns detects a transcript
  that is confidently wrong because the annotation is incomplete (Ma and
  Kingsford 2019).
- Use the terms of section 3.

**U3. Report and check the resampling degrees of freedom.**
- Where: `add_uncertainty()` (after the draws are read),
  `uncertainty_report()` (settings table), `uncertainty_methods()`
  (sentence), `salmon_run_qc()` (already returns `replicate_type` and
  `num_replicates`).
- What: compute `N * (B - 1)`; inform when below 60 (Baldoni, Chen and Smyth
  2024), and mention that Gibbs draws resolve low counts better when
  `replicate_type` is `"bootstrap"`. A message, not an error.
- Scope note: the rule is about the overdispersion estimate. InfRV and the
  correlation are per sample and depend on `B` itself; the stability article
  (PLAN workstream L, 30 against 100 draws) is the evidence for those.
- Tests: message appears for 2 samples x 10 draws, not for 6 x 30.

**U4. Terminus collapse score per pair.**
- Where: `inferential_correlation()` / `isoform_ambiguity()`; the draws array
  from `draws_array()` is already in hand.
- What: for each sibling pair, per sample,
  `score = InfRV(g_i + g_j) - (InfRV(g_i) + InfRV(g_j)) / 2`, averaged over
  samples; new column `collapse_gain` (negative = the sum is better
  determined than the parts).
- Why: correlation says the two trade reads; the score says how much
  variance the sum removes. Terminus decides on the score, not on the
  correlation. It needs no Terminus run.
- Check: on the knockdown data compare with the groups of `run_terminus()`
  (`data-raw/knockdown/compare_terminus.R`); expect pairs inside a Terminus
  group to have the most negative scores.
- Cost: one more pass over the pairs; pairs are already generated 50,000 at
  a time.

**U5. Range of optima from equivalence classes.** Derivation.
- Idea. Under the equivalence-class likelihood the data enter only through
  the class totals. With `A[c, t]` = conditional probability of class `c`
  given transcript `t` (the weights of `--dumpEqWeights`; 1 for members
  when weights are absent) and `x` = estimated counts, any `x' >= 0` with
  `A x' = A x` has the same likelihood. The range of transcript `t` is
  `[min x'_t, max x'_t]` over that set: two linear programs per transcript,
  as in Zheng et al. 2022 section 2.2, with class totals in the place of
  their junction and phasing-path totals.
- Not established: that Salmon classes are an adequate stand-in for the
  paper's sufficient statistics, and how effective-length and bias terms
  enter. Validate against `subgraphquant` on a few genes before exposing
  anything.
- Cheap first step, no solver: per connected component of the class graph,
  the null space of `A` (base `qr()` or `svd()`). A transcript whose
  coordinate is zero in every null vector is identifiable; report
  `identifiable` (logical) and the component nullity. The range itself needs
  the non-negativity constraints, hence a solver.
- Where: extend `read_eq_classes()` to keep the weights (a line holds
  `2 * size + 2` fields when weights are present); new internal
  `equivalence_matrix()`; new exported `abundance_range()` or columns
  `range_lo`, `range_hi`, `identifiable` in `tx_uncertainty()`; store per
  sample like `equivalence_support()` does.
- Dependencies: a solver in Suggests (`highs` or `lpSolve`), behind
  `rlang::check_installed()`. Components are small (one gene, sometimes a
  family), so problems are tiny; the cost is the count: transcripts x
  samples x 2. Restrict to expressed transcripts in components with
  nullity above zero.
- Limits to document: complete-reference case only (`k = 0`); the
  incomplete-reference range needs graph quantification and is out of
  reach from Salmon output.
- Payoff: the only item that measures meaning (2) of section 1. Expect a
  minority of transcripts flagged (35 to 50% in the paper, under different
  sufficient statistics).

**U6. Data-adaptive correlation cutoff.**
- Where: `uncertainty_defaults()`, `add_uncertainty(cor_cutoff = )`.
- What: `cor_cutoff = "adaptive"` sets the cutoff to minus the 97.5th
  percentile of the per-transcript maximum mean correlation, the mmcollapse
  default.
- Risk: changes the ambiguity groups and so the confidence classes of every
  analysis; keep -0.5 as the default and offer the adaptive value as an
  option. Report the value used in `uncertainty_methods()`.

**U7. DE reliability flag by interval overlap.** Derivation.
- Zheng et al. 2022 flag a DE call when the mean ranges of the two groups
  overlap by more than 25% of the smaller range. With draws in the place of
  ranges: per transcript, the interval of the group mean over draws (2.5 to
  97.5%) for each group, overlap divided by the shorter interval.
- Where: `compute_de_confidence()`; new column `interval_overlap`.
- This is not the paper's quantity (draw intervals are not ranges of
  optima); name it accordingly. With U5 done, the same code gives the
  paper's flag from `range_lo` / `range_hi`.
- Overlaps with `sign_agreement`; check on the knockdown data whether it adds
  anything before keeping it.

**U8. Low-information flag.**
- Where: `build_uncertainty()`, from columns already present.
- What: `low_information = unique_fraction == 0 & mean_inf_rv >= cutoff`,
  the mmcollapse selection (no unique reads and high posterior spread).
  Needs equivalence classes; `NA` without them.

**U9. Long-read expectation.**
- Where: `isoformic-uncertainty.qmd` (Oarfish section),
  `uncertainty_report()`.
- What: state that overdispersion close to 1 is expected for long reads
  (Baldoni 2024) and that Wissel et al. 2026 found short-read and Kinnex
  quantification concordant once inferential variability is corrected.
  Report the median overdispersion in `uncertainty_report()` so that the
  expectation can be checked.

**U10. Gene-level separability from the annotation.**
- Zhang et al. 2026 bound the quantification error with the condition
  number of a gene-specific matrix built from the isoform structures and the
  read-length distribution; Hiller et al. 2009 and Ferrer-Bonsoms et al.
  2022 give identifiability conditions. Only abstracts were read.
- Would give a design-time score ("can 2 x 100 bp reads separate the
  isoforms of this gene at all") from `exons` alone, next to the features of
  the stability article (PLAN workstream L, `isoform_features.R`).
- Read the full texts before estimating the work.

### 5.3 Applied on 2026-11-11

The effect sizes below were measured on the GENCODE 49 quantifications of
the knockdown data, before the bundled data moved to GENCODE 50 (commit
f3a1a87); they were not rerun and describe the change, not the current
tables. The current numbers of the articles are checked against the
GENCODE 50 bundles. Where the evaluation of section 8 (GENCODE 50) measured
the same thing, its number is given in the item.

- **U1.** `inferential_rv()` is now
  `pmax(var - mean, 0) / (mean + pseudocount) + shift`; `cut_transcript_trees()`
  calls it instead of its inline copy; `data-raw/knockdown/isoform_features.R`
  has the same form. A test compares it with `fishpond::computeInfRV()`;
  on GENCODE 50 the snapshot used in section 8 matches `computeInfRV()` on
  the placenta data to a largest relative difference of 1.29e-11.
  Effect on the six-sample knockdown data (GENCODE 49) (116,229 assessed transcripts):
  median InfRV of the groups up by 0.01; among the large changes 2 isoforms
  moved from high to moderate confidence (203 to 201); one isoform changed
  InfRV bin in the replication table; uncertain isoforms left alone by the
  correlation rule 78,892 to 79,208 (316 more are uncertain); TreeTerminus
  cut: 5,543 to 5,496 groups. Rebuilt with the dev package:
  `inst/extdata/knockdown/stability/*` (`isoform_features.R`; R2 of the log
  InfRV model 0.629 to 0.634, of the disagreement model 0.079 to 0.062) and
  `read_split/terminus_comparison`, `replication_rv`, `grouping_comparison`
  (`compare_terminus.R`, `validate_uncertainty.R`, `compare_treeterminus.R`;
  `null_calls`, `null_stable`, `replication_class` and `isar_comparison`
  came out identical). Article text changed where a stated number moved.
  The old read-split cache is kept as `read_split_contrasts.old_infrv.rds`
  in the knockdown work directory.
- **U2.** `?add_uncertainty`: bullet renamed "Ambiguity between isoforms",
  mmcollapse credited, a "What is not measured" paragraph, four references
  added. Same in `uncertainty_methods()` and the uncertainty article;
  `library.bib` has `turro2014`, `zheng2022`, `ma2019`, `baldoni2024gibbs`.
- **U3.** `resampling_df()` and `inform_resampling_df()`; `add_uncertainty()`
  stores `metadata$uncertainty_resampling_df` and informs below 60 (class
  `isoformic_few_inferential_replicates`); `uncertainty_methods()` states the
  number; `uncertainty_report()$replicates` is new. Section 8 checked the
  sufficiency of 30 draws on GENCODE 50: InfRV thirds from 15 and 30 draws
  agree for 93.5% of the expressed placenta transcripts, the default class
  flag (class not `high`) for 98.0% of the replication calls, and the overdispersion prior moves
  by +2.0% from 30 to 100 draws. Whether 100 draws change the calls was not
  tested.
- **U9.** The article and `?add_uncertainty` state the long-read expectation;
  `uncertainty_report()$replicates$median_overdispersion` lets it be checked
  (needs the draws in the object).

### 5.4 Applied or tested on 2026-11-12

Checked on the six GSE129540 samples (`quants_sub`, GENCODE 49, 116,229
assessed transcripts, 639,909 sibling pairs; Terminus groups from
`terminus_all/out`). The script was a scratch file and is not in the
repository; the numbers come from one run. They are GENCODE 49 numbers and
were not rerun on GENCODE 50; the identifiability items (U5) have GENCODE 50
counterparts from section 8, given at the end of U5.

- **U4, done.** `isoform_ambiguity()$collapse_gain`, computed in
  `inferential_correlation()` from the two variances and the correlation, so
  no second pass over the draws. Pairs inside one Terminus group (32,142):
  median -1.41, 58% negative. Other pairs (607,767): median 0.86, 13%
  negative. As a ranking of "same Terminus group" it is no better than the
  correlation (AUC 0.67 against 0.68), and the two are only moderately
  related (Spearman 0.37), so it is a second view and not a replacement.
  4,122 of the 4,131 ambiguous pairs have a negative gain; 95,477 pairs have
  a negative gain without being ambiguous.
- **U6, done as an option.** `add_uncertainty(cor_cutoff = "adaptive")`,
  resolved by `resolve_cor_cutoff()`; default unchanged. Adaptive value
  -0.246: ambiguous isoforms 8,076 to 36,910, groups 3,945 to 13,329. Of
  the pairs inside a Terminus group, the fixed cutoff calls 2,716 ambiguous
  and the adaptive one 8,918; outside, 1,415 and 15,275. When no positive
  percentile exists the default is used with a message
  (`isoformic_no_adaptive_cutoff`).
- **U8, done.** `tx_uncertainty()$low_information`; 47,116 of the 87,251
  uncertain isoforms have no unique fragment. Only 4,621 of them are in an
  ambiguity group of the correlation rule.
- **U7, not adopted.** Interval overlap of the group means over draws
  (2.5 to 97.5%, overlap over the shorter interval, flag above 25%) against
  `sign_agreement >= 0.95`: Spearman -0.82; the two agree on 110,523 of
  116,229 isoforms. 5,410 stable isoforms overlap by more than 25%; they
  have a higher InfRV (median 2.06 against 1.26), and the interval of the
  fold-change excludes zero for 68% of them against 98%. `lfc_lo`/`lfc_hi`
  already show that, so the column would add a second threshold for the
  same information. No truth was used; the read-split replication
  (`validate_uncertainty.R`) could test whether these 5,410 replicate less.
  The paper's own quantity (ranges of optima) is not available from
  Salmon output; see U5 below.
- **U5, first step applied.** Structural identifiability from
  the equivalence classes of one knockdown sample (`GSE129540_SRR8870780`,
  12 M read pairs, GENCODE 49; prototype first, then the package function).
  - Method. `M` = class by transcript matrix; a transcript is identifiable
    when its unit vector lies in the row space of `M`. (1) Peel: a class
    with a single unresolved member resolves it; repeat. (2) Merge
    transcripts with identical columns (not identifiable from each other)
    and peel again. (3) Split the rest into connected components; dense
    `svd()` up to 400 columns, and for larger ones the projection of two
    random vectors on the null space by CGLS with `Matrix` (a transcript
    with a non-zero coordinate is not identifiable). CGLS agreed with the
    clean gap expected (no value between 1e-10 and 1e-4).
  - Size. 470,395 classes, 315,803 transcripts in a class. The class graph
    has one component of 118,498 transcripts; peeling resolves 74,133
    transcripts in 17 rounds and merging marks 77,589; the largest block
    left has 30,898 columns (18,000 CGLS iterations per vector). About
    4 minutes per sample, 2.5 of them in that block.
  - Result without weights: 130,468 of 315,803 transcripts (41%) are not
    identifiable; 25,956 of the 115,655 transcripts of the six-sample
    uncertainty table that are in a class of this sample (22%).
  - With weights (`salmon quant --dumpEqWeights`, same reads, rerun in
    `quants_eqw/` of the knockdown work directory; 527,324 classes): of the
    193,096 unresolved columns in components up to 400 columns, 87,708 are
    still not identifiable, and 51,514 singular values are below 1e-14 of
    the largest. So the weights do not remove it: isoforms that share all
    their fragments have proportional columns (the ratio is the ratio of
    effective lengths). Another 12,000 singular values lie between 1e-12
    and 1e-3 without a gap, so near-dependence is common as well and any
    rank tolerance is a choice.
  - It is not what the draws measure. Among those transcripts the median
    InfRV is 2.31 for the non-identifiable against 2.25 for the others, the
    share called uncertain 73% against 76%, and 998 of 25,956 are in an
    ambiguity group against 7,076 of 89,699. `unique_fraction` is zero for
    83% of them against 37%.
  - Why the range of optima does not follow (derivation, not validated).
    The classes in the file are the observed ones; the likelihood also
    depends on the normalisation `sum(theta * effective_length)`. A move
    along the null space of `M` that lowers this sum lets every class
    probability be scaled up, so it raises the likelihood. The maximum is
    therefore the point of the set `{M theta = b, theta >= 0}` with the
    smallest length-weighted total: a linear program whose optimum is in
    general one vertex. Under the full model the range of a structurally
    non-identifiable transcript is then a single point, and the estimate is
    decided by effective length (the shorter isoform takes the fragments),
    not by reads. Zheng et al. 2022 get ranges because their sufficient
    statistics are junction and path totals without that term; whether that
    is their exact model was not checked against `subgraphquant`.
  - Replication test (decides the item). Each of the six GSE129540 samples
    has its reads split in two halves (`quants_half/`); error = absolute
    log2 ratio of the estimated counts (plus 1) of the two halves;
    identifiability from the classes of the 12 M pair run of the same
    sample; transcripts of the six-sample uncertainty table, binned by mean
    count. Median error, non-identifiable against identifiable: 2.32
    against 1.94 at 1 to 5 reads, 3.96 against 1.82 at 5 to 20, 2.41
    against 0.70 at 20 to 100 (16,101 against 78,820 sample-transcripts),
    0.30 against 0.20 at 100 to 1,000. Per sample at 20 to 100 reads the
    medians agree to the second decimal. Among transcripts the draws call
    uncertain: 2.44 against 0.83 at 20 to 100 reads. Among those the draws
    call stable: 3.52 against 0.50 at 5 to 20 reads (1,512 against 18,203);
    at higher counts almost no stable transcript is non-identifiable (75 at
    20 to 100). So the flag carries information the draws do not: it splits
    the "uncertain" class into a part that replicates and a part that does
    not.
  - Across samples: 13,790 of the 116,389 table transcripts are
    non-identifiable in all six samples and 39,681 in at least one; of the
    13,790, 73% are uncertain by draws (75% of the others) and all have
    `unique_fraction == 0`.
  - Decision: adopted on 2026-11-12 as `equivalence_identifiability()`
    (`R/identifiability.R`), a separate function and not a column of
    `add_uncertainty()`: 30 to 185 s per sample here, and `Matrix` (now in
    `Suggests`) for blocks above 400 columns. The function reproduces the
    prototype (130,468 for the first sample). Membership only; the name
    says "identifiability", not "range of optima".
  - Used in the confidence since 2026-11-12 (decision of the maintainer):
    `add_uncertainty(identifiability = TRUE)` stores `identifiable` and
    `identifiable_fraction`, and `de_confidence()` counts a
    non-identifiable isoform as a fourth failed condition. Off by default
    (cost; the article numbers do not change). How many confident calls of
    GSE129540 this removes was not measured; section 8 measured the flag on
    its own instead (next item).
  - GENCODE 50, placenta (section 8; 25 samples, 2026-11-18). 22,181 of the
    57,294 expressed transcripts (38.7%) are not identifiable in at least
    one sample, in 10,631 of 14,844 genes; the function took 203 minutes
    for the 25 samples (about 8 minutes per sample, against 30 to 185 s on
    the 12 M pair knockdown samples). In the held-out replication,
    non-identifiable calls (20.4% of the calls) replicate in 51.6% against
    65.5% (adjusted OR 0.56, 0.54-0.58, planned bootstrap). In the
    7 v 7 null they do not attract more false calls (adjusted RR 0.83,
    0.43-1.57, sensitivity model). So the flag predicts non-replication more
    than the class does, but less than `sign_agreement` (OR 0.36).
  - Open: a weighted variant;
    whether the flip between halves is the effective-length mechanism
    above. The GENCODE 49 numbers above are one run on one dataset.
  - `read_eq_classes()` reads files with weights correctly (it takes the
    first `size` fields as members and the last as the count) and ignores
    the weights.

### 5.5 Relevant but not implementable from the inputs the package reads

- **Coverage anomalies (SAD)** and **junction coverage compatibility
  (JCC)**: need per-base or per-junction coverage from alignments. Salmon
  output has none. Document as a limit; a user with BAM files can run SAD or
  JCC outside and join the result through `add_annotation_columns()`.
- **Perplexity**: needs held-out reads and a second quantification.
- **Incomplete-reference ranges** (graph quantification): need the splice
  graph flow, not available from Salmon.
- **mehenDi**: testing on the TreeTerminus tree; available in the
  `fishpond`/`mehenDi` ecosystem, could be wrapped like `run_terminus()` if
  wanted. Not assessed.
- **sleuth / CompDTUme**: alternative ways to propagate uncertainty; the
  package already has Swish (external) and divided counts. No gap.

## 6. Search log

- Consensus: 4 searches (inferential uncertainty and DTE; non-identifiability
  and Kingsford; edgeR divided counts; long-read uncertainty), 15 to 30
  results each.
- scite: 3 topic searches, 3 targeted full-text excerpt reads (Baldoni 2023
  preprint, Baldoni 2024 NARGAB, Zheng 2022), `verify_papers` on 20 DOIs.
  One DOI guessed for Ferrer-Bonsoms et al. 2022 resolved to an unrelated
  paper and was dropped.
- Undermind: 4 semantic searches, 1 deep search (91 papers), 7 PDFs read by
  sub-agents (Baldoni 2023, Baldoni 2024, Zheng 2019/2022, Sarkar 2020, Turro
  2014, Ma and Kingsford 2019 bioRxiv, Fan 2022).
- Not found or not checked: a published calibration study of Salmon Gibbs
  against bootstrap variance versus true error other than Baldoni 2024; the
  Zhu 2019 full text (InfRV formula taken from Terminus and from the fishpond
  source).

## 7. References

Verified DOIs.

- Baldoni PL et al. Dividing out quantification uncertainty allows efficient
  assessment of differential transcript expression with edgeR. Nucleic Acids
  Res 2024;52:e13. doi:10.1093/nar/gkad1167 (preprint
  doi:10.1101/2023.04.02.535231)
- Baldoni PL, Chen L, Smyth GK. Faster and more accurate assessment of
  differential transcript expression with Gibbs sampling and edgeR v4. NAR
  Genom Bioinform 2024;6:lqae151. doi:10.1093/nargab/lqae151
- Baldoni PL et al. Dividing out quantification uncertainty enables
  assessment of differential transcript usage with limma and edgeR. Nucleic
  Acids Res 2025. doi:10.1093/nar/gkaf1305
- Chen Y et al. edgeR v4. Nucleic Acids Res 2025. doi:10.1093/nar/gkaf018
- Zhu A et al. Nonparametric expression analysis using inferential replicate
  counts. Nucleic Acids Res 2019. doi:10.1093/nar/gkz622
- Pimentel H et al. Differential analysis of RNA-seq incorporating
  quantification uncertainty. Nat Methods 2017. doi:10.1038/nmeth.4324
- Sarkar H et al. Terminus enables the discovery of data-driven, robust
  transcript groups from RNA-seq data. Bioinformatics 2020.
  doi:10.1093/bioinformatics/btaa448
- Singh N, Love MI, Patro R. TreeTerminus. iScience 2023.
  doi:10.1016/j.isci.2023.106961
- Singh N et al. Tree-based differential testing using inferential
  uncertainty for RNA-seq. Genome Res 2025. doi:10.1101/gr.279981.124
- Turro E, Astle WJ, Tavare S. Flexible analysis of RNA-seq data using mixed
  effects models. Bioinformatics 2014;30:180. doi:10.1093/bioinformatics/btt624
- Zheng H, Ma C, Kingsford C. Deriving ranges of optimal estimated transcript
  expression due to nonidentifiability. J Comput Biol 2022;29:121.
  doi:10.1089/cmb.2021.0444
- Ma C, Zheng H, Kingsford C. Exact transcript quantification over splice
  graphs. Algorithms Mol Biol 2021. doi:10.1186/s13015-021-00184-7
- Ma C, Kingsford C. Detecting, categorizing, and correcting coverage
  anomalies of RNA-seq quantification. Cell Syst 2019.
  doi:10.1016/j.cels.2019.10.005
- Fan J, Chan S, Patro R. Perplexity: evaluating transcript abundance
  estimation in the absence of ground truth. Algorithms Mol Biol 2022.
  doi:10.1186/s13015-022-00214-y
- Young A, Van Buren S, Rashid NU. Differential transcript usage analysis
  incorporating quantification uncertainty via compositional measurement
  error regression modeling. Biostatistics 2023.
  doi:10.1093/biostatistics/kxad008
- Zhang B et al. Quantifying transcript complexity via the condition number
  of gene-specific random matrix. Brief Bioinform 2026. doi:10.1093/bib/bbag134
- Jousheghani ZZ, Singh N, Patro R. Oarfish. Bioinformatics 2025.
  doi:10.1093/bioinformatics/btaf240
- Wissel D et al. A systematic benchmark of high-accuracy PacBio long-read
  RNA sequencing for transcript-level quantification. Genome Biol 2026.
  doi:10.1186/s13059-026-03988-1
- Patro R et al. Salmon provides fast and bias-aware quantification of
  transcript expression. Nat Methods 2017. doi:10.1038/nmeth.4197

## 8. Genome-wide evaluation (2026-11-18)

A pre-registered evaluation of whether quantification uncertainty changes
transcript-level results, and whether the flags of `add_uncertainty()` and
`de_confidence()` predict which calls fail to replicate. Negative results
are kept. The article
`vignettes/articles/isoformic-uncertainty-evidence.qmd` presents it; every
number there is computed inline from the bundled tables.

### 8.1 Where things are

- Scripts: `data-raw/uncertainty_evidence/` (numbered in run order; the
  README gives the run order, the inputs and the environment variables).
  `prereg.yaml` is the frozen plan (seeds, methods, thresholds, rules R-A
  and R-B); `prereg_amendments.yaml` lists every change, A1 to A24; A5 to A9
  and A11 to A24 were made or recorded after the results were seen and say
  so (`post_hoc: true`; A10 changes no analysis). A18 records the F4 tables
  recomputed with the repository code instead of the snapshot, and A19 that
  the mouse null was not pre-registered. Neither file is in an external registry, and both
  were committed together with the results, so the git history does not
  show that the plan came first: that rests on `written` in `prereg.yaml`
  and the snapshot sha256 it records (A16). The A4 pilot wrote the null
  results of its 9 splits before A4 fixed the design; they are kept in the
  reported results (the only DESeq2 rows at 5 v 5 and 3 v 3), and nothing
  records whether they were looked at (A16).
- Code: every script loads a frozen `git archive` snapshot of commit
  4a6e983 (its sha256 is in `prereg.yaml`), not the live `R/`, and stamps
  every table with the git sha, the snapshot sha256, the seed and the
  script.
  `92_package_flag_check.R` (post hoc, A24) is the link to the live code:
  it recomputes the sign agreement of the discovery half of all 60
  directed pairs with the repository's `add_uncertainty()` and
  `de_confidence()` and counts the calls whose `unstable` differs from the
  evaluated flag (keys `B.pkgcheck.*`; 0 of 167,318 on 2026-11-20). Rerun
  it after any change to `compute_de_confidence()` and the evidence
  numbers still describe the package only if the count stays 0.
- Outputs: everything is written under `$ISOFORMIC_UE_OUT` (not in the
  repository): `S0/`, `store/`, `null_placenta/`, `placenta_audit/`,
  `null_mouse/`, `rhalf/`, `wpc_summary/`, `od_prior/`, `article_tables/`,
  and a summary report `REPORT.md`.
- Article data: `90_article_tables.R` copies every stated number with its
  CI, CI type and source table into `article_tables/`;
  `91_bundle_article_data.R` copies those tables into
  `vignettes/articles/data/uncertainty_evidence/`, pruned to the keys the
  articles use (about 19 KB).

### 8.2 Design

- Data (all GENCODE 50, Salmon, 30 Gibbs draws): placenta GSE148241 /
  PRJNA622166, 8 preeclampsia and 17 controls (15 control units: two
  labelled A/B pairs, probably twins, keep one member each); mouse adBNST
  PRJNA984605, 10 control males (GENCODE vM39); random read halves of the
  two knockdown series GSE129540 and GSE75489.
- Universe: 57,294 transcripts of 14,844 genes expressed in the
  preeclampsia samples; InfRV thirds from those samples only (cut points
  2.17 and 7.59).
- WP-A, placenta null: control-only splits 3 v 3, 5 v 5 (100 each) and
  7 v 7 (210 leave-one-unit-out); methods M1 (edgeR QL, point counts), M2
  (divided counts), M6 (swish on draws), M7 (swish on tied point counts, a
  control), M9a (M1 kept in class `high`), M10a (M1 after dropping as many
  low-expression transcripts), M11 (gene level), plus M3/M4 (limma-voom,
  7 v 7 only) and M5 (DESeq2 Wald, 20 7 v 7 splits; 3 pilot splits at 5 v 5
  and 3 v 3); DTU (diffSplice, satuRn, DEXSeq) and NMD; a gene-thinning
  power arm on 30 7 v 7 and 30 5 v 5 splits (5% of the filtered
  multi-isoform genes, sampled within InfRV thirds, plus as many
  single-isoform genes, each with |log2FC| uniform in [0.58, 2]). Rule R-A:
  adjusted RR (false calls p < 0.01, top against bottom InfRV third,
  quasi-Poisson with expression decile and isoform class) above 1, removed
  by an uncertainty-aware method, absent at gene level. prereg.yaml names
  no group size; the analysis judges it at 7 v 7, the design A4 gave the
  most splits (A13, recorded after the results).
- WP-B, placenta real contrast: prevalence of fragile calls; held-out
  replication over 30 half-splits in both directions (60 directed pairs,
  167,318 calls at q < 0.10, a script rule with a top-500 fallback that
  prereg.yaml does not contain, A14; replicated = same sign at one-sided
  p < 0.05 in the other half); flags from the discovery half only, except
  `unique_fraction == 0` (all PE samples). Rule R-B: the default flag
  (class not `high`) has an adjusted OR CI below 1 and beats the
  count-matched flag. prereg.yaml describes that control as having the same
  mean-count distribution; it was implemented as C4, the same number of
  lowest-count calls per pair, a different control (A12, recorded after the
  results). prereg.yaml registers only the class flag; sign agreement and
  the other flags are secondary, and making sign agreement the headline was
  a post hoc choice (its 0.95 threshold was fixed).
- WP-C: mouse null (all 126 5 v 5 splits; not in prereg.yaml, which names
  no mouse dataset or design: it applies the R-A rule by analogy, A19); read
  halves (technical).
- WP-D: the overdispersion prior genome-wide against the bundled 39-gene
  panel prior.
- Intervals: placenta null, delete-one-unit jackknife with a Monte Carlo
  correction where valid (A5; the raw jackknife is the conservative
  interval and is reported where the correction removes most of the
  variance); mouse, exact jackknife over all splits per animal (A8); WP-B,
  the planned two-level bootstrap (half-splits, then genes within each
  pair; prereg.yaml says transcripts, A15), plus gene-clustered and
  split x gene bootstraps (A11);
  WP-D, bootstrap over genes.

### 8.3 Verdicts

- **R-A fails** at 7 v 7, the size the analysis judged it at (A13): M1
  adjusted RR 1.00 (raw jackknife 0.56-1.76; the MC-corrected 0.92-1.08
  keeps only about 2% of the jackknife variance, within the noise of the MC
  term, and is not reliable; same for 5 v 5, 1.81, raw 1.17-2.80). The raw
  interval does not rule out an excess of up to about 1.8. edgeR QL on point
  counts is calibrated there (inflation at p < 0.01 0.83; 0.095 false
  calls per split at q < 0.05), as is limma-voom (0.076); DESeq2 Wald is
  not (43.9 false calls per split, raw jackknife 0-115, inflation 1.82,
  20 splits). swish makes more false calls than M1 in the placenta too
  (5 v 5: 2.39, 1.71-3.07, against 0.11; 7 v 7: 0.67 against 0.095).
  Nothing above 7 per group was tested.
- **Small n (secondary sizes, not judged by R-A)**: 5 v 5 1.81 (raw jackknife
  1.17-2.80), 3 v 3 2.44 (1.65-3.60); overall inflation stays near 1 (1.13
  at 3 v 3), so false calls are redistributed towards uncertain isoforms,
  not added. Mouse 5 v 5: 2.70 (1.68-4.33), passing R-A only under the A8
  deviations (the mouse design was not pre-registered, A19); the mouse controls hold a hidden library batch, and the
  excess survives on the batch-balanced splits (3.86, 1.31-11.4).
- **Divided counts** remove the small-n concentration (ratio of RRs M2/M1
  0.45 at 5 v 5, 0.28 at 3 v 3) but do not reduce false discoveries (3 v 3:
  0.20 against 1.21 per split at q < 0.05, 0.39 without M2's worst split)
  and cost power (thinning: 21,349 against 31,886 true calls; 2,507 against
  7,254 in the top InfRV third).
- **R-B passes formally but is fragile**: class-not-high OR 0.68
  (0.647-0.726); it beats C4 by 0.018 (0.004-0.031) under the planned
  bootstrap, but C4 is not the control prereg.yaml describes (A12), and the
  CI includes 0 with genes as clusters (-0.018 to 0.053)
  or split x gene (-0.019 to 0.057). Closer to the planned control, the
  class keeps a matched deficit within count x |z| cells of 0.099
  (0.087-0.110), and the class permuted within those cells has -0.017. A
  same-size flag from count and |z|
  (C4b, post hoc) beats the class (deficit 0.200 against 0.085). In the
  null, class `high` holds 8.4% of tested transcripts but 18.5% of false
  calls (enrichment 2.22), because 93% of null calls are stable across
  draws.
- **`sign_agreement < 0.95` predicts non-replication** (secondary flag,
  made the headline after the results; strongest by OR in the placenta
  half-splits but narrow, it catches about 10% of failures; in the read
  halves other flags do as well or better): 5.7% of the calls of the
  discovery halves (q < 0.10), 2.1% of the 6,636 full-cohort calls at
  q < 0.05;
  replication 34.7% against 64.3% (pooled half-split rates; per pair
  24.0-51.1% against 51.3-85.9%, lower in all 60 pairs); deficit matched
  within count x |z| cells 0.249 (0.232-0.266); adjusted OR 0.36,
  planned-bootstrap CI 0.33-0.39, gene-clustered 0.31-0.41; pair fixed effects
  0.35; splines on count and |z| 0.36; sign-concordance outcome 0.29. PPV
  for non-replication 65% (base 37%), catches 10% of failures, removes 3.2%
  of replicating calls. Read halves (technical consistency only): OR 0.069
  (0.030-0.130) GSE129540 and 0.154 (0.079-0.268) GSE75489, the latter on
  fallback calls (29 and 45 calls at q < 0.10, so the top 500 by p).
- **InfRV >= 0.5** marks 87.1% of the expressed placenta transcripts; the
  top InfRV third has OR 0.63; an ambiguous partner 0.87; no unique reads
  0.69; not identifiable 0.56.
- **Point against draws** (post hoc): 23.0% of the 6,636 full-contrast DTE
  calls (q < 0.05) have a point log2FC outside the 2.5-97.5% range of their
  30 per-draw refits, 2.0% by at least FLT1-214's 1.53 (FLT1-214: 5.28
  against 2.51-3.75; LEP-203, 1.51, is just short of it and not among
  the 134); 1.8% of class `high`, 27.0% `moderate`, 39.4% `low`.
  In 1,239 of the 1,529 the point estimate is further from 0 than the whole
  range: the per-draw refits are pulled towards 0 relative to the EM point
  estimate, a systematic EM against Gibbs difference; neither is ground
  truth, and the draw range is not a CI for the point estimate. 6 of the 8
  tested FLT1 isoforms lie outside their range (5 further from 0; FLT1-204
  by 0.01 only). The refits are edgeR refits per draw; `de_confidence()`'s
  `lfc_median`/`lfc_lo`/`lfc_hi` (difference of the group means of
  log2(scaled count + 5), a log ratio of geometric means) were
  not compared with the point estimate.
- **Story genes** (prereg `story_genes`; `story_genes.parquet`): ENG
  (gene log2FC 1.88, replicated 59 of 59) and LEP (6.79, 60 of 60); ENG-202
  (12,749 reads) and ENG-216 (3,927) are class `low` yet replicate 42 of 42
  and 60 of 60. The prereg fragile rule selected no gene (0 candidates); a
  substitute rule of 32_story.R (F2/F4, not in prereg.yaml) picked LAS1L
  (of 120): gene point log2FC 0.40 while only 3% of the per-draw refits
  share its sign, yet the call replicates in 14 of 16 pairs; LAS1L-202
  (point 4.74, refits -0.96 to 1.61) replicates 40 of 40. Replication uses
  the same EM point counts, so it cannot catch a distortion repeated in all
  samples. The robust rule picked TBC1D26 (of 411; one isoform, class
  `high`, 60 of 60).
- **DTU**: scaling helps diffSplice at 5 v 5 (3.49 -> 0.17 false genes per
  split, intervals separated), satuRn inflation 0.71 (0.15-3.48) -> 1.25
  (0.13-11.68) cannot be told apart, DEXSeq is anti-conservative either
  way.
- **Prior**: genome-wide 1.477 (1.447-1.504); the bundled panel prior is
  2.94 times the genome-wide one for GSE129540 but swapping them changes no
  DTE call.

### 8.4 Corrections made to the articles

`isoformic-uncertainty.qmd` (2026-11-19): the class gradient, the "all 72
below 0.5" scale, the "about three times chance" excess on raw counts and
the usage gains are now qualified as technical-half results with the
biological numbers next to them; the reporting advice no longer implies
filtering on the class and names `sign_agreement` as the column to report;
`confident_only` is described as a sensitivity analysis; the cut-off
caution says which default was supported. The stability article sentence
"uncertainty is not a problem of low counts" is only partly supported (the
flags keep an effect after adjustment for count and |z|, but a count and
|z| flag beats the class and the largest distortions are in low-count
isoforms); on 2026-11-20 it was narrowed to the InfRV result ("A high
InfRV is not a sign of low counts"), with a pointer to the low-count
exceptions of the evidence article.

`isoformic-uncertainty-evidence.qmd` (2026-11-19, on review): the
pre-registration is described as it is (two rules; sign agreement a
secondary flag made the headline after the results; files in the
repository, not a registry); R-A's group size, R-B's control and the
discovery rule are described as implemented, with amendments A12 to A14;
LEP-203 is "close to" FLT1-214, not among the 134; the mouse batch-aligned
false calls are attributed to the one fully confounded split; the thinning
design and the recommendation (7 per group, edgeR QL only) are stated as
tested; a supplementary section gives every null method (limma-voom, DESeq2,
Wilcoxon, M9/M10 variants, gene level), the 5 v 5 thinning arm and the
robustness variants. `isoformic-uncertainty.qmd`: sign agreement is the
best, not the only, predictor; C4 is described as the lowest-count calls.
`REPORT.md` in the output directory is not tracked; its stale numbers
(enrichment CI 1.72-2.85, from before A6; "about 28" reads instead of 27;
"satuRn gets worse") were corrected on 2026-11-19, but
`workflow_results.json` still has the old enrichment CI: quote the bundled
tables, which are the source of record.

`isoformic-uncertainty-evidence.qmd` (2026-11-19, second review): sign
agreement is described as a secondary flag of the analysis scripts, with
only its threshold in the plan; the plan/pilot/results order is stated as
the files show it (A16); the replication bootstrap is the "planned
bootstrap" with its gene unit (A15); every story gene of the plan and both
rules are reported; the point-against-draws text no longer says the
uncertainty only changes the effect size, and the advice names the edgeR
refits that were evaluated and says `lfc_median` was not; the satuRn,
thinning FDR, thinning design, robustness, draws and placenta-size
sentences match the tables. `isoformic-uncertainty.qmd`: the satuRn
sentence and the `lfc_median` advice follow, and an HTML comment says where
each typed Terminus/TreeTerminus/IsoformSwitchAnalyzeR number comes from.

`isoformic-uncertainty-evidence.qmd` (2026-11-19, third review): the
read-half check now shows every flag (table `tbl-rhalf-flags`): there
sign agreement is not the best flag, and InfRV, "uncertain" and the
ambiguous partner have adjusted ORs as low or lower (the ambiguous partner
only after adjustment), so the read halves support "the draws carry
information", not sign agreement in particular; "best" is limited to the
placenta half-splits, by OR and replication gap (by delta-AUC over count
and |z| the non-identifiable flag is level with it). The plan's two other
R-B controls are reported: the placebo as implemented (the class shuffled
within count x |z| cells) and the label-permutation baseline C1 (deficit
0.096 under permuted labels against 0.085 observed; not a null with 4 PE
samples per discovery half), with amendment A17. WP-D is in the text:
divided counts made 7,887 against 6,235 calls on the common tested set
(about as many with their own filters, 7,889 against 7,975), of unknown
truth; the panel prior (2.9 and 2.4 times the genome-wide one) changed no
knockdown call. The sign-agreement headline and advice give its
sensitivity (10% of failures) and the 93% of null false calls stable over
the draws; the small-sample excess is "at nominal p < 0.01" with the q <
0.05 false calls (0.20 and 0.11 per split) next to it. Smaller fixes: the
fragile and robust rules use "validation halves" (not held out for a cohort
call), the substitute rule counts 120 candidate results (117 genes), the
TBC1D26 numbers are given for the gene and for TBC1D26-204, the FLT1 gap
sentence is hedged, the tutorial pointer says "two knockdown series", and
Limits says the prior was the only quantity computed from 100 draws
(+2.0%). `isoformic-uncertainty.qmd`: "best predicted" is limited to the
placenta half-splits with its 10% sensitivity, "thirty replicates ...
change little" is replaced by what was and was not tested (and the
stability article's 100-draw advice), the class is no longer described as
how far a call "can be trusted", and the panel-prior and real-contrast
WP-D results are stated. The note that `R/uncertainty.R` still called the
class a measure of trust was closed by U11 to U13 (8.5, 2026-11-20).

Articles aligned with U11 to U13 (2026-11-20). `isoformic-uncertainty.qmd`
reports `unstable` / `unstable_de` / `flags` instead of the class, states
the measured effect (OR 0.36, gene-clustered 0.31-0.41, 60 of 60 pairs) with
its 10% sensitivity, that the flag was chosen after the results, and that
R-A failed at 7 v 7 (the 3 v 3 excess is nominal, p < 0.01, only); it
describes `uncertain` as a description of the isoform (87% of expressed
placenta isoforms), says to read `unstable` on the calls only (null
transcripts have a median log2FC near 0), and presents `drop_unstable` as
a sensitivity analysis. The knockdown scripts count by flag instead of
class and were rerun with the working tree on GENCODE 50
(`validate_uncertainty.R`, `compare_terminus.R`, `compare_treeterminus.R`):
`replication_class` became `replication_flag` (per flag: unstable 15 of
110 half-1 calls, same direction 0.80 against 0.96, none called again);
`null_stable` counts uncertain and ambiguous isoforms (1,189 of 1,192 stable
large null changes are uncertain); `terminus_comparison` and
`grouping_comparison` count the large changes with an unstable direction
(1,698 of 5,110, the same under every grouping) and with an ambiguous
sibling (1,003 correlation, 2,867 Terminus, 1,061 TreeTerminus); every
other value of these tables, `null_calls` and `replication_rv` came out
identical. `isar_comparison` lost its `high_confidence` column (the rest is
unchanged; IsoformSwitchAnalyzeR was not rerun). The new read-split cache
replaces the old one in the GENCODE 50 work directory; the old one is kept
as `read_split_contrasts.pre_unstable.rds`. `isoformic-uncertainty-evidence.qmd`
calls the class "former default" or "of the evaluated snapshot", replaces
the "defaults predate this evaluation" note with the change, and adds to
Limits that the package flag reproduces the evaluated flag on all 60 pairs
and is a secondary, post hoc flag whose effect elsewhere is untested. The
stability, introduction, preprocess and enrichment articles, the README and
`_pkgdown.yml` no longer mention a confidence of the calls. The untracked
evidence `REPORT.md` was synced with the article (read halves are a
technical check, not two independent series; gene-clustered intervals; R-A
fails; small-n excess nominal only).

Fourth review (2026-11-20). The 5.7% share of unstable calls is that of
the discovery halves (q < 0.10); the evidence article and the uncertainty
tutorial now give the full-cohort share next to it (2.1% of 6,636 calls at
q < 0.05, key `B.prevalence.F1.pct_sign_agree_lt_0.95`, added to
`90_article_tables.R` from `summary_family.parquet`). The 2.22-fold excess of
class `high` among null calls is no longer read as the class working the
wrong way: a stable direction alone is over-represented 3.00 times (31% of
tested transcripts, 93% of null calls), since selecting calls at p < 0.01
selects stable signs; the evidence against the class in the null is the
equal-removal comparison (M9a minus M10a -0.097, -0.24 to 0.05). The usage
section of the tutorial says that scaling raised satuRn's null calls at
q < 0.05 from 0 to 6 and 2, about as many as its 11 extra real calls.
`REPORT.md`: the FLT1-214 example no longer says the uncertainty changes
only the size of the effect or that LEP-203 is among the 134, and the
Baldoni pattern is reported as holding only in part (no power gain from
divided counts).

### 8.5 Package changes (applied 2026-11-20)

Status: U11, U12 and U13 applied on 2026-11-20 (spec of that day, kept
outside the repository). The articles, the README and the knockdown
read-split tables were aligned the same day (8.4, last paragraph).

- **U11 applied.** `de_confidence()` reports `unstable`
  (`sign_agreement < min_agreement`, `NA` when not assessed) right after
  `sign_agreement`, and `unstable_de` (a significant call with an unstable
  direction); `reason` is renamed `flags` ("unstable direction" first).
  `uncertainty_report()$de` has one row per flag (`unstable`, `uncertain`,
  `ambiguous`, `not_identifiable` when computed) with its rule and counts;
  `uncertainty_methods()` describes the sign agreement and makes no trust
  claim; `plot_de_confidence()` colours by the direction over the
  replicates. The help text states the measured effect (OR 0.36,
  gene-clustered 0.31-0.41, 60 of 60 pairs, about 10% of failures caught,
  93% of null false calls stable) and that the headline was chosen after the
  results.
- **U12 applied: the class was removed**, not redefined (it was never
  released: v0.1.2 exports none of these functions). `confidence` and
  `confident_de` are gone; `isoform_switches(confident_only = )` is defunct
  (error of class `isoformic_defunct_argument`; the formal goes after the
  0.1.3 release), `drop_unstable = TRUE` leaves out the unstable calls and
  keeps the not-assessed ones, `n_det_unstable` replaces `n_det_confident`,
  `uncertainty_report()$switches` compares all calls with the calls without
  unstable directions, and `confidence_palette()` is removed. Identifiability
  stays a separate column, never combined with `unstable` (combining weak
  flags is what diluted the class OR from 0.36 to 0.68). Tables stored by
  earlier development builds are upgraded by `read_isoformic()` (one message,
  class `isoformic_legacy_confidence`) and silently in memory.
- **U13: kept 0.5, no longer a gate; quantile cut not adopted**
  (object-dependent, changes with panels and `x[i, ]`; top-third OR 0.63 is
  weaker than `unstable` and non-identifiability; `cut_transcript_trees()`
  relies on the absolute 0.5). The docs now say it marks InfRV above the
  sampling noise of about 30 draws and describes the isoform without judging
  calls; `uncertainty_report()$transcripts` gains `q25_inf_rv` and
  `q75_inf_rv`.
- The per-draw range (`lfc_lo`-`lfc_hi`) is documented as the spread of the
  screening statistic, not a confidence interval for the point estimate of
  the primary test, and the sign agreement as measured against the median
  draw sign, not `det_log2FC`.
- The fold-change of `de_confidence()` (`lfc_median`, `lfc_lo`, `lfc_hi`,
  and so `sign_agreement`) is documented as what the code and the evaluated
  rule (`m9_flags()` of `21_null_dte.R`) compute: the difference between
  the groups of the mean of log2(scaled count + 5), a log ratio of
  geometric means on the scale of the Welch test. Earlier text called it a
  log ratio of group means; the two can differ in sign. A unit test pins
  `welch_matrix()` to it.

The original proposals:

| ID | Change | Evidence (8.3) | Effort |
|---|---|---|---|
| U11 | Report `sign_agreement < min_agreement` as the primary flag: a count in `uncertainty_report()$de`, a sentence in `uncertainty_methods()`, highlighted in `plot_de_confidence()` | OR 0.36 (gene-clustered 0.31-0.41), robust to power adjustment and pair fixed effects; read halves 0.07-0.15 | S |
| U12 | Change the role of the class in `confident_de` (`R/uncertainty.R`, `n_det_sig > 0 & confidence == "high"`) and in `isoform_switches(confident_only = TRUE)` / `uncertainty_report()$switches`: base it on sign agreement (optionally identifiability), or keep the class as a description only and drop the filter wording | class flags 82.5% of calls, OR 0.68, beats count-matched C4 by 0.018 with a CI that includes 0 under gene clustering, beaten by a count and |z| flag; `high` enriched 2.22 times for null calls; at equal removal no better than dropping low-expression transcripts | M (changes numbers in articles and tests) |
| U13 | Revisit `inf_rv_cutoff = 0.5`: a data-relative cut (for example the top InfRV third of the expressed transcripts), or keep InfRV out of the class | 87.1% of expressed placenta transcripts are above 0.5; the bottom third ends at 2.17; top-third flag OR 0.63; `cor_cutoff = -0.5` (ambiguous partner) OR 0.87 | S-M |

Also proposed, and applied with U11: state in `?de_confidence` that the
per-draw fold-change range is not a confidence interval for the point
estimate, and that the point estimate of low-count isoforms can lie outside
it (8.3, point against draws). Not tested here: long reads, `delta_z`, `lfc_z`,
`dominant_inf_rv`, and more than 30 draws for the calls.
