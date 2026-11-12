# Quantification uncertainty: literature notes

Written 2026-11-11. Sources: Consensus, scite and Undermind searches run that
day. Purpose: ground the definitions in `R/uncertainty.R`,
`R/usage.R` and the uncertainty article, and list what the package does not
cover.

**Where to start.** Section 1 is the conclusion. Section 4 maps each
definition to the function that implements it. Section 5 is the backlog:
ten items (U1 to U10) with the function to change, the effect, the tests and
the effort; U1 to U3 and U9 are small and can go in one pass. Section 2 holds
the formulas and numbers behind them.

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
2026-11-11 its documentation says so. See section 4.

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
| Overdispersion | `compute_overdispersion()` (`R/uncertainty.R`), used by `inferential_overdispersion()`, `scale_counts()`, `overdispersion_for()` and `test_usage()` | Baldoni 2024 estimator (2.1) | Same: Pearson statistic over draws, `df = n_expressed * (B - 1)`, prior `max(1, median(phi) / qf(0.5, median(df), 3))`, moderated, floored at 1. Terms with zero mean are left out of sum and df, a case the paper does not specify |
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
| Misfit | none | SAD, JCC, perplexity (2.6) | Not implemented, see section 5.3 |

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

- **U1.** `inferential_rv()` is now
  `pmax(var - mean, 0) / (mean + pseudocount) + shift`; `cut_transcript_trees()`
  calls it instead of its inline copy; `data-raw/knockdown/isoform_features.R`
  has the same form. A test compares it with `fishpond::computeInfRV()`.
  Effect on the six-sample knockdown data (116,229 assessed transcripts):
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
  number; `uncertainty_report()$replicates` is new.
- **U9.** The article and `?add_uncertainty` state the long-read expectation;
  `uncertainty_report()$replicates$median_overdispersion` lets it be checked
  (needs the draws in the object).

### 5.4 Applied or tested on 2026-11-12

Checked on the six GSE129540 samples (`quants_sub`, GENCODE 49, 116,229
assessed transcripts, 639,909 sibling pairs; Terminus groups from
`terminus_all/out`). The script was a scratch file and is not in the
repository; the numbers come from one run.

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
  - Open: use the flag in `transcript_confidence()` (a non-identifiable
    transcript should probably not be "confident"); a weighted variant;
    whether the flip between halves is the effective-length mechanism
    above. The numbers above are one run on one dataset.
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
