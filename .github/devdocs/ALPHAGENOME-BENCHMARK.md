# AlphaGenome pilot benchmark and GENCODE 50 bump: notes

Written 2026-11-11. Status (2026-11-18): pilot finished on 500 genes and
four conditions (2026-11-12), and rerun on 2026-11-17 with five seeds and a
bootstrap over genes after the coding model was retrained (section 4, step
3, "Rerun of 2026-11-17"). On the Karousis truth, the only one that does not
come from the annotation, AlphaGenome gives no detectable gain in pooled AUC
(95% interval -0.016 to +0.013). On the retained-intron biotype it gives a
small gain whose interval excludes zero (+0.007, 0.004 to 0.010). The stop
rule as written ("no improvement of the `retained_intron` or
`protein_coding` rows") is therefore not met for the retained-intron row.
Nothing is integrated, but that is now a judgement (reasons in section 4,
step 3), not the outcome of the rule; the maintainer should confirm it. A6
stays closed on that basis; A5 stays open without AlphaGenome.

**Where to start.** Section 4 "To resume" has the commands that rerun the
pilot (an API key is needed). Section 7 is the backlog (A1 to A9) with the
status of each item; two defects were found on the way (A1, A2; fixed, with
A7, on 2026-11-11).

## 1. Question

Do sequence-based splicing predictions of AlphaGenome add information on
productive versus unproductive isoforms beyond the structure and sequence
rules of `transcript_productivity()`?

## 2. What AlphaGenome can and cannot give

- Input: DNA sequence, 16 kb to 1 Mb (`SEQUENCE_LENGTH_1MB = 1048576`).
  Organisms: human and mouse.
- Outputs relevant here: `SPLICE_SITES`, `SPLICE_SITE_USAGE`,
  `SPLICE_JUNCTIONS`, `RNA_SEQ`, `CAGE`. Others exist (ATAC, DNase, ChIP,
  contact maps, methylation, single-cell outputs).
- It has no notion of reading frame, translation or decay. It cannot score
  coding potential or NMD directly. The coding potential model
  (`inst/extdata/coding_model/`) works on the mature transcript, which
  AlphaGenome never sees.
- Useful as indirect evidence: splice-site strength and junction support
  (then listed as not implemented in `ARCHITECTURE.md`; a first-order
  splice-site score is implemented since, A4), intron retention
  propensity from predicted intronic coverage, TSS / TES support.
- Circularity: it is trained on steady-state RNA-seq, so junctions of
  NMD-degraded isoforms are already under-represented in predicted usage. Low
  predicted usage mixes "weakly spliced" with "spliced then degraded". Use as
  supporting evidence, not as an NMD caller.

## 3. Client facts (alphagenome 2.0.0, inspected in the installed package)

- The current API is `alphagenome.models.v2.dna_client` ("AlphaGenome 2").
  `alphagenome.models` exports nothing at top level; a v1 module also exists.
- `dna_client.create(api_key)` connects to `gdmscience.googleapis.com`. An
  API key is required. No local weights were tried (the workstation GPU is a
  Quadro P4000, 8 GB).
- `predict_interval(interval, requested_outputs=..., condition=...)` takes
  **one condition per request**: a `ModelCondition` or a composite curie
  string built from CLO, UBERON, CL, EFO or NTR terms.
  `supported_conditions()` lists them with the imputation mode of each
  output type (observed or imputed). Request count = genes x conditions.
- `SPLICE_JUNCTIONS` returns `JunctionData`: the junctions the model predicts
  (`genome.Junction`, 0-based start and end, stranded) with one value per
  track. Annotated junctions that are absent are below detection.
- Track outputs are `TrackData` (`values[positions, tracks]`, `metadata` with
  `name` and `strand`, `resolution`).
- Verified with a key on 2026-11-12 (first genes of the pilot):
  - `supported_conditions()` lists 2,799 conditions; 227 have observed
    `splice_sites`, `splice_site_usage` and `splice_junctions`, and all 227
    also have observed `rna_seq`.
  - Within an observed condition some tracks are still imputed (for HepG2
    the GTEx/recount3 polyA track is `imputed_exploratory`): filter on the
    `imputation_mode` column of the track metadata, not on the condition.
  - Splice-site signal: with `pos` the first or last base of the intron
    (1-based), the peak is on the neighbouring exonic base in the frame of
    `02_score.py`: offset -1 at the left end of the intron and +1 at the
    right end, for both strands. `03_features_alphagenome.R` asserts it.
  - Junctions: `intron_start = junction_start + 1`, `intron_end =
    junction_end`. Every canonical intron of the first 22 genes is among
    the predicted junctions, and 40% of the others.
  - One request of 1 Mb with four outputs takes about 5 s when the server
    answers at once; the full run averaged closer to 50 s per gene of four
    conditions. No quota error was met.
- From the documentation of 2.0.0 (released 2026-11-11, read 2026-11-12:
  alphagenomedocs.com, user guides "AG1 vs AG2" and "imputation", CHANGELOG):
  - AG2 takes one `condition=` per call where AG1 took `ontology_terms`;
    `output_metadata()` became `supported_conditions()` and
    `observed_metadata()`; `ontology_curie` became `biosample_curie`; every
    track has `imputation_mode`. `TrackData` and `JunctionData` have
    `filter_to_observed()`.
  - `observed`: "This exact channel was experimentally measured in this
    condition and used for training." `imputed_supported`: inferred in a
    benchmarked direction. `imputed_exploratory`: "not benchmarked and can
    be poor or systematically misleading." Splicing outputs are always
    exploratory when not measured; the guide says to restrict splicing
    scoring to conditions where the splicing tracks or bulk RNA-seq were
    measured.
  - `SPLICE_SITES` is now condition-specific (five channels per condition:
    donor and acceptor on each strand, and padding); AG1 had four
    condition-agnostic tracks. `PROCAP` was removed.
  - `dna_client.create(api_key, model_version=)`: `ModelVersion` has
    `ALL_FOLDS` and `FOLD_1`; the default is not documented. With all folds
    the pilot genes were in the training data, so the predictions at
    annotated sites are partly recall of measured RNA-seq. Which regions
    `FOLD_1` holds out was not found.
  - Terms: free for non-commercial use; outputs "should not be used for the
    training of other machine learning models"; nothing found on
    redistributing predictions (the terms page itself was not read). This
    weighs against shipping predictions as frozen data (A6).
  - Scorers `SpliceJunctionScorer` and `GeneMaskSplicingScorer` exist for
    variants; not used here.
- Not verified: the behaviour of `condition = None`, quotas, and the mouse
  assembly.

## 4. Pilot design

Files: `data-raw/alphagenome/` (the scripts, `pixi.toml` and `pixi.lock`
are tracked by the `.gitignore` exceptions for `data-raw`; outputs are not).
Own pixi project (`pixi.toml`, `pixi.lock`; python 3.12, pyarrow, pandas,
alphagenome from PyPI). Work directory: `ISOFORMIC_ALPHAGENOME_DIR`, on a
data disk; every output named below is relative to it.

### Step 1, `01_targets.R` (run on GENCODE v50)

- Eligible genes: protein coding, chr1-22 and X, no haplotype, length at most
  900 kb, a canonical transcript, at least 2 protein-coding transcripts.
- Stratum `karousis` (250 genes): genes with a validated NMD / non-NMD
  transcript pair of Karousis et al. 2021 (`nmd_exons` frozen dataset) where
  both Ensembl IDs are in the annotation. Many IDs of that table are
  `MSTRG.*` (StringTie) and cannot be used.
- Stratum `biotype` (250 genes): random genes with at least one
  `nonsense_mediated_decay` and one `retained_intron` transcript.
- Seed 1. Window: 1 Mb centred on the gene, shifted to stay inside the
  chromosome (lengths from `annot_seq_regions.parquet`; the first run used
  UCSC `hg38.chrom.sizes` because of defect A1, and the targets are identical).
- Output in `<work_dir>/targets/`: `genes`, `transcripts`, `exons`,
  `introns`, `sites`, `regions`, `karousis_pairs` (Parquet).
  `ISOFORMIC_ALPHAGENOME_TARGETS` writes to another sub-directory, to
  rebuild the targets with newer code and compare them with the scored ones.

Result: 500 genes, 18,378 transcripts, 19,226 unique introns, 252 NMD and 257
non-NMD labelled transcripts; 259 Karousis pairs can be evaluated.

Rules-only baseline on the sample (`use_biotype = FALSE`,
`use_evidence = FALSE`), with the GENCODE 50 coding model (targets rebuilt
2026-11-17 in `targets_2026-11-18/`; the first targets of 2026-11-11, with
the GENCODE 49 model, in brackets where they differ):

| Biotype | n | unproductive | productive | altered_productive | uncertain |
|---|---|---|---|---|---|
| `nonsense_mediated_decay` | 5,042 | 5,031 | 0 | 10 | 1 |
| `retained_intron` | 1,781 | 1,110 [1,136] | 3 | 478 [463] | 190 [179] |
| `protein_coding` | 10,439 | 40 | 2,309 | 7,011 | 1,079 |
| `protein_coding_CDS_not_defined` | 1,106 | 818 [820] | 2 | 216 | 70 [68] |
| `non_stop_decay` | 5 | 5 | 0 | 0 | 0 |
| `protein_coding_LoF` | 2 | 2 | 0 | 0 | 0 |
| `TEC` | 3 | 2 | 0 | 1 | 0 |

So: NMD biotype recall 0.998, retained intron 0.62 (1,110 / 1,781; was
0.64), protein coding called productive or altered 0.89. Room for
AlphaGenome is in the retained-intron class and the 478 retained-intron
transcripts called `altered_productive`.

The rebuilt targets are identical to the scored ones (genes, transcripts,
exons, introns, sites, regions, Karousis pairs) except for the `rule_*`
columns: 52 calls in 41 genes changed, none of a Karousis-labelled
transcript, so `02_score.py` was not rerun. All 52 come from the coding
model: the rules with the current code and the GENCODE 49 model (ORF cache
off) give the old calls exactly, and with the GENCODE 50 model the new ones;
coding probabilities near the 0.364 cutoff crossed it in both directions.

### Step 2, `02_score.py` (done 2026-11-12, 500 genes, about 3 hours)

Conditions of the pilot: HepG2 (`EFO:0001187`), K562 (`EFO:0002067`),
frontal cortex (`UBERON:0001870`) and liver (`UBERON:0002107`); all four have
observed splicing and RNA-seq tracks. Smoke test passed. Log:
`score.log` in the work directory; the run resumes when restarted with the
same command.

One request per gene and condition for the four outputs. Keeps only values at
the targets, nothing aggregated over tracks:

- `scores/sites/<gene>.parquet`: splice-site and usage tracks at each
  annotated site, offsets -2..2 (the offset with the signal is settled in
  step 3).
- `scores/regions/<gene>.parquet`: mean RNA-seq coverage per exon and intron.
- `scores/junctions/<gene>.parquet`: predicted junctions inside the gene, on
  its strand, in AlphaGenome coordinates.
- `scores/tracks.parquet`: track metadata.

Resumable (the junction file marks a gene as done), retries with backoff,
`--list-conditions` and `--smoke` modes. Only `site_scores()`,
`region_scores()` and `metadata_table()` were checked, on synthetic tracks.

### Step 3, `03_evaluate.R` (done 2026-11-12; rerun 2026-11-17)

The script holds the baseline that the AlphaGenome features have to beat:
the calls of the rules and the splice-site strength of
`splice_site_strength()` (backlog A4), with the truths, the grouped
cross-validation and the paired tests already in place. AlphaGenome features
go in as further entries of `feature_sets`.

Design (since the rerun; the header of the script has the details):

- Logistic regression, 5 folds grouped by gene, pooled out-of-fold AUC.
  Seeds 1 to 5 of the fold assignment (`ISOFORMIC_ALPHAGENOME_SEEDS`); seed
  1 keeps the draws of the first single run.
- Uncertainty: bootstrap over genes (2,000 replicates, seed 20261117,
  `ISOFORMIC_ALPHAGENOME_BOOT`): genes drawn with replacement, all
  transcripts of a drawn gene kept, applied to the out-of-fold predictions.
  Each replicate averages the AUC over the five seeds, so the interval is
  that of the five-seed mean. Differences between feature sets are paired
  (the same draw for every set and seed of a replicate).
- Within a truth every feature set is fitted on the same transcripts: those
  with a value for every feature. Single-exon transcripts have no
  splice-site score and are left out (Karousis 0, NMD biotype 3, retained
  intron 32). The first version fitted each set on its own complete cases,
  so its rules row of the retained-intron truth had 12,157 transcripts
  against 12,125 for the others.
- `ISOFORMIC_ALPHAGENOME_SCORES=false` runs the baseline alone on all 500
  genes even when the AlphaGenome scores exist.
- Writes `features`, `auc`, `auc_summary`, `auc_contrasts`, `predictions`,
  `pairs` and `diagnostics` (Parquet) to `<work_dir>/evaluation/`, or to the
  sub-directory named by `ISOFORMIC_ALPHAGENOME_EVALUATION`.

Baseline on the v50 targets, all 500 genes (rerun of 2026-11-17 in
`evaluation_2026-11-18_baseline/`, on the rebuilt targets). Each cell: AUC
with seed 1; mean over the five seeds (bootstrap 95% interval of the mean).

| Truth | n (positive) | Rules | PWM | Intron number | Rules + PWM | Baseline |
|---|---|---|---|---|---|---|
| Karousis NMD against non-NMD | 509 (252) | 0.784; 0.780 (0.736-0.821) | 0.729; 0.735 (0.692-0.778) | 0.578; 0.578 (0.550-0.608) | 0.805; 0.807 (0.769-0.843) | 0.799; 0.800 (0.760-0.836) |
| NMD against protein coding biotype | 15,478 (5,042) | 0.996; 0.996 (0.994-0.998) | 0.599; 0.604 (0.588-0.622) | 0.517; 0.517 (0.494-0.540) | 0.997; 0.997 (0.995-0.998) | 0.997; 0.997 (0.995-0.998) |
| Retained intron against protein coding biotype | 12,188 (1,752) | 0.863; 0.865 (0.853-0.877) | 0.840; 0.839 (0.820-0.858) | 0.796; 0.798 (0.773-0.824) | 0.962; 0.962 (0.954-0.968) | 0.962; 0.963 (0.955-0.970) |

- Baseline: rules, PWM and intron number.
- Against the first single run (2026-11-12) seed 1 reproduces every value
  except two of the retained-intron row: rules 0.870 to 0.863 and rules +
  PWM 0.963 to 0.962. Recomputed with the same fold draws: the old calls on
  the old 12,220 transcripts give 0.8697 for the rules, the new calls on
  them 0.8648 (coding model), the old calls on the common 12,188 0.8682
  (transcript set), the new calls on the common set 0.8632. Rules + PWM
  was always fitted on 12,188: 0.9625 to 0.9619, the coding model alone.
- Differences (mean over seeds, 95% interval, share of replicates at or
  below zero): rules + PWM minus rules, Karousis +0.027 (0.012 to 0.045,
  0), NMD biotype +0.0006 (0.0000 to 0.0013, 0.029), retained intron
  +0.097 (0.086 to 0.107, 0). Baseline minus rules + PWM, Karousis -0.0077
  (-0.0124 to -0.0031): the intron-number term lowers the AUC there;
  retained intron +0.0015 (-0.0007 to 0.0039, 0.091).
- PWM features: weakest donor and acceptor over all introns, the same over
  the introns that the canonical transcript of the gene does not have
  ("novel"), and whether there is a novel intron. Canonical transcripts are
  `protein_coding`, the negative class of both biotype truths, so they have
  no novel intron by construction.
- Karousis pairs (259, same gene): the NMD transcript has the weaker novel
  donor in 142 of 233 pairs that differ (61%, sign test p = 0.001) and the
  weaker novel acceptor in 149 of 244 (61%, p = 0.0007). Over all introns the
  direction is the opposite for the donor (35 of 107, 33%): the restriction
  to novel introns is what carries the signal. Unchanged by the rerun.
- Reading. On the only truth that is independent of the annotation
  (Karousis) the sequence score adds about 0.03 AUC to the rules (interval
  0.012 to 0.045). The NMD biotype row is circular for the rules (the
  biotype follows the same 50 nt rule) and leaves no room. The
  retained-intron row is the one with room (0.863 to 0.962), but intron
  number alone gives 0.796 there: retained-intron models are short
  fragments with few introns, so most of that gain is transcript structure,
  not splice-site sequence.
- Consequence for the comparison below: AlphaGenome has to add to rules +
  PWM (0.805 on Karousis, 0.962 on retained introns). Requiring it to beat
  0.805 with the intron-number term in the model, as the first version of
  this note did, was not a fair bar: that term costs 0.008 on Karousis.

Result with the AlphaGenome features (`03_features_alphagenome.R`; observed
tracks only). 497 of the 500 genes have a predicted junction on their strand
and enter the models; the baseline is refitted on them. Rerun of 2026-11-17
(`evaluation_2026-11-18/`, log `evaluate_2026-11-18.log`); the first single
run is kept in `evaluation/` and `evaluate.log`. Each cell: AUC with seed 1;
mean over the five seeds (bootstrap 95% interval of the mean).

| Truth | n (positive) | Rules + PWM | Baseline | Baseline + contained intron | AlphaGenome alone | Rules + AlphaGenome | Baseline + AlphaGenome |
|---|---|---|---|---|---|---|---|
| Karousis NMD against non-NMD | 509 (252) | 0.805; 0.807 (0.769-0.843) | 0.799; 0.800 (0.760-0.836) | 0.796; 0.796 (0.753-0.835) | 0.731; 0.726 (0.678-0.772) | 0.802; 0.807 (0.766-0.845) | 0.802; 0.799 (0.755-0.838) |
| NMD against protein coding biotype | 15,405 (5,025) | 0.997; 0.997 (0.995-0.998) | 0.997; 0.997 (0.995-0.998) | 0.997; 0.997 (0.995-0.998) | 0.590; 0.596 (0.581-0.611) | 0.998; 0.997 (0.996-0.998) | 0.997; 0.997 (0.996-0.998) |
| Retained intron against protein coding biotype | 12,125 (1,745) | 0.962; 0.961 (0.953-0.969) | 0.965; 0.964 (0.957-0.971) | 0.969; 0.969 (0.961-0.976) | 0.858; 0.861 (0.844-0.877) | 0.967; 0.967 (0.961-0.973) | 0.971; 0.971 (0.964-0.978) |

Paired differences: mean over the five seeds [range over the seeds], 95%
interval, share of the 2,000 replicates at or below zero.

| Contrast | Karousis | NMD biotype | Retained intron |
|---|---|---|---|
| Baseline + AlphaGenome minus baseline | -0.0010 [-0.0075, 0.0029], -0.0158 to 0.0133, 0.54 | +0.0002 [-0.0003, 0.0006], -0.0002 to 0.0006, 0.13 | +0.0069 [0.0048, 0.0081], 0.0044 to 0.0096, 0 |
| Baseline + AlphaGenome minus baseline + contained intron | +0.0023 [-0.0115, 0.0109], -0.0095 to 0.0154, 0.34 | +0.0001 [-0.0001, 0.0005], -0.0003 to 0.0006, 0.23 | +0.0026 [0.0005, 0.0040], 0.0006 to 0.0047, 0.0055 |
| Baseline + contained intron minus baseline | -0.0034 [-0.0082, 0.0040], -0.0117 to 0.0032, 0.82 | +0.0001 [-0.0002, 0.0003], -0.0001 to 0.0002, 0.28 | +0.0043 [0.0036, 0.0049], 0.0026 to 0.0061, 0 |
| Rules + AlphaGenome minus rules | +0.0269 [0.0172, 0.0353], 0.0117 to 0.0427, 0 | +0.0010 [0.0006, 0.0014], 0.0004 to 0.0018, 0 | +0.1025 [0.0995, 0.1057], 0.0924 to 0.1133, 0 |
| Rules + PWM minus rules | +0.0271 [0.0211, 0.0289], 0.0116 to 0.0445, 0 | +0.0007 [0.0003, 0.0009], 0.0001 to 0.0013, 0.0035 | +0.0974 [0.0949, 0.1004], 0.0866 to 0.1086, 0 |
| Baseline minus rules + PWM | -0.0077 [-0.0153, -0.0038], -0.0124 to -0.0031, 0.999 | +0.0000 [-0.0003, 0.0004], -0.0003 to 0.0003, 0.47 | +0.0029 [0.0020, 0.0036], 0.0006 to 0.0052, 0.007 |

- Baseline: rules, PWM and intron number. "Contained intron" is the
  structural flag that a canonical intron lies inside an exon of the
  transcript, without any prediction. Baseline + AlphaGenome holds the
  baseline, the flag and the AlphaGenome features.
- A separate check, not in `03_evaluate.R` (other fold draws, 500
  replicates, so its absolute AUCs differ by about 0.005), left out the
  intron-number term: rules + PWM + AlphaGenome minus rules + PWM is
  +0.0078 (0.0047 to 0.0105) on retained introns and +0.0018 (-0.0097 to
  0.0138) on Karousis; with the contained-intron flag in both models,
  +0.0045 (0.0020 to 0.0069) on retained introns, 0 of 500 replicates at or
  below zero.
- AlphaGenome feature groups alone, Karousis / retained intron, seed 1
  (five-seed mean): junction 0.704 / 0.829 (0.704 / 0.837), site and usage
  0.675 / 0.774 (0.678 / 0.775), retention ratio 0.496 / 0.554 (0.495 /
  0.558; interval on retained introns 0.533 to 0.584); the flag alone 0.496
  / 0.631 (0.500 / 0.632).
- Diagnostics: 98.6% of the canonical introns and 37.5% of the others are
  among the predicted junctions; median site probability 0.893 against
  0.015, usage 0.663 against 0.006; 8 observed junction tracks.
- Karousis pairs: the NMD transcript has the lower AlphaGenome score more
  often than the lower PWM score: weakest junction 176 of 241 (73%), novel
  junction 182 of 245 (74%, p = 1.4e-14), novel site 163 of 241 (68%), novel
  usage 167 of 244 (68%), against 61% for the PWM. All sign tests p < 1e-5.
  In the 223 pairs where both the AlphaGenome novel junction and the PWM
  novel donor differ, AlphaGenome orders the pair right in 165 and the PWM
  in 137 (McNemar test p = 0.0036): within a gene AlphaGenome carries
  information that the PWM lacks, and the pooled logistic regression does
  not use it.
- Reading.
  - Karousis: no detectable gain in pooled AUC. The interval of baseline +
    AlphaGenome minus baseline (-0.016 to +0.013) does not rule out a gain
    of about 0.013 to 0.015. With 509 transcripts and up to 17 correlated
    features in an unpenalised `glm`, overfitting pulls the larger feature
    sets toward the null. The pair result above points the other way, so
    "no added value" on this truth needs a within-pair test (for example
    conditional logistic regression on the pairs, rules + PWM with and
    without AlphaGenome), which was not run.
  - Retained intron: a small gain with an interval that excludes zero,
    +0.007 over the baseline and +0.003 beyond the contained-intron flag
    (+0.0045 in the check without the intron-number term). The +0.0069 and
    +0.0045 contrasts have no replicate at or below zero; the +0.0026
    contrast (share 0.0055) would not survive a Bonferroni correction over
    the 18 contrast-truth rows.
  - NMD biotype: circular for the rules (AUC 0.996), so its small
    "significant" differences (rules + PWM minus rules +0.0007, rules +
    AlphaGenome minus rules +0.0010) are ceiling effects and count neither
    way.
  - The `protein_coding` row named in the stop rule has no truth here and
    was not evaluated.
  - The predicted retention ratio carries almost no signal for the
    retained-intron biotype (0.554).
- Decision. The stop rule ("if the `retained_intron` and `protein_coding`
  rows do not improve, do not integrate") is not met: the retained-intron
  row improves. No integration remains defensible as a judgement, not by
  the rule: the gain is small; the truth is a biotype label, not
  productivity (478 retained-intron transcripts are `altered_productive` by
  the rules, and the contained-intron flag is close to the definition of
  the label); the default model version was trained on all folds with steady-state RNA-seq,
  so it can recall measured RNA-seq and decay-lowered junction usage of
  these genes, which makes the positive results upper bounds; and the
  terms are non-commercial (section 3). A6 stays closed on that basis (to
  be confirmed by the maintainer); A5 cannot use the predicted coverage
  ratio.
- Limits: the bootstrap keeps the out-of-fold predictions fixed, so the
  intervals leave out the variance of refitting and are somewhat narrow;
  the share of replicates at or below zero is not a p-value; logistic
  regression on transcript-level minima wastes the pairwise signal (shown
  above); the GTEx isoform-fraction control was not run, and with a
  positive retained-intron result it is needed to tell splicing signal from
  an expression proxy; the retention feature is one design (canonical
  intron inside an exon over canonical exon coverage) and was not tuned.

#### Rerun of 2026-11-17

Why: the coding model was retrained on GENCODE 50 (A8) after the pilot,
which can move the rule calls, and the first evaluation had one seed and no
intervals. Run on 2026-11-17 and 2026-11-18 at commit 4a6e983, with no
uncommitted change under `R/` or `inst/`.

- Targets rebuilt (`targets_2026-11-18/`): identical apart from 52 rule
  calls, all from the coding model (step 1). `02_score.py` was not rerun.
- Evaluation with five seeds and 2,000 bootstrap replicates, on the rebuilt
  targets, with and without the AlphaGenome features (tables above).
  Seed 1 reproduces every Karousis AUC and every sign test of the first run.
  The retained-intron rows moved through the coding model and the common
  transcript set: rules 0.867 to 0.862, rules + PWM 0.963 to 0.962,
  baseline 0.966 to 0.965, baseline + AlphaGenome 0.972 to 0.971. Two NMD
  rows differ in the fifth decimal (rules, through the 3 transcripts left
  out; the retention ratio, cause not found); equal at three decimals.
- A7 rules benchmark rerun with `data-raw/benchmark_rules.R` (A7), A4
  tables recomputed with `splice_site_groups.R` (A4).
- Peak memory: targets 2.2 GB, evaluation 2.0 GB and 1.5 GB, A4 tables
  3.4 GB, A7 benchmark 4.2 GB (cached ORF table).

### To resume

```sh
export ALPHAGENOME_API_KEY=...        # from the git-ignored .env of the package root
cd data-raw/alphagenome
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --list-conditions
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --smoke
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --conditions <curies>
```

The R steps run from the package root in the pixi environment of the
package, for example (the names of the rerun):

```sh
ISOFORMIC_ALPHAGENOME_TARGETS=targets_2026-11-18 \
  pixi run --frozen Rscript data-raw/alphagenome/01_targets.R
ISOFORMIC_ALPHAGENOME_TARGETS=targets_2026-11-18 \
  ISOFORMIC_ALPHAGENOME_EVALUATION=evaluation_2026-11-18 \
  pixi run --frozen Rscript data-raw/alphagenome/03_evaluate.R
```

`ISOFORMIC_ALPHAGENOME_DIR` must be set for all of them. `test_v49/` in the
work directory is a disposable trial of step 1 on v49.

## 5. Integration, if the pilot is positive

Not at runtime. Follow the frozen-evidence pattern: score all junctions and
introns of the pinned GENCODE release once, write Parquet keyed by junction
and transcript, register it in `R/frozen-data.R`, join with an
`add_alphagenome_splicing()` in `R/productivity-evidence.R` through
`add_annotation_columns()`, and read the columns by name in
`transcript_productivity()` as an optional evidence tier, fetch gated by
`ISOFORMIC_FULL`. Check the licence for redistribution first.

## 6. GENCODE 50 / M39 bump (done the same day)

- Latest releases on the EBI server: human 50 (Ensembl 116, dated
  2026-04-08), mouse M39. GENCODE 49 is Ensembl 115.
- `download_reference()` defaults are now 50 and M39; roxygen, `man/`,
  `README.md`, `NEWS.md`, the generic lines of the preprocess article and
  `tests/testthat/test-download_reference.R` follow. That test file and
  `test-frozen.R` pass.
- GENCODE 50 `chr_patch_hapl_scaff`: 670,670 transcripts, 646,577 with
  sequence from the primary assembly genome (v49: 533,740 and 509,650).
  GFF3 is 166 MB; EBI served it at about 100 kB/s.
- **APPRIS naming.** In `e113v49` the `e113` is the Ensembl release and `v49`
  is the APPRIS annotation version. Ensembl 113 is GENCODE 47, so the old pin
  was labelled "GENCODE 49" by mistake. New pin: release `2026_08.v50`,
  directory `e116v50` (Ensembl 116 = GENCODE 50).
- Rebuilt frozen files (local cache only, not uploaded to the `frozen-data`
  release or Zenodo): `appris_e116v50.parquet`, 430,885 transcripts, MD5
  `c84750680dace8250ed99db78543f21c` (was 172,157 transcripts);
  `appris_domains_e116v50.parquet`, 612,203 rows over 286,888 transcripts,
  MD5 `f6f399a9779b1c01b01a08a49089573f` (was 179,111 rows). Same schema, no
  duplicate rows.
- **Source cache trap.** `fetch_frozen_data(build = "always")` reuses a
  cached source file by name ("already exists"), so changing only the URL in
  the registry rebuilt the old data with the old MD5. The old sources were
  moved to `sources/appris_e113v49` and `sources/appris_domains_e113v49` in
  the isoformic cache before the rebuild. A version-keyed source directory
  would remove the trap.
- Left on v49 on purpose: bundled example data and the mini annotation, the
  knockdown and case-study bundles (tied to existing quantifications), and
  the coding potential model (`data-raw/build_coding_model.R`,
  `release <- "49"`; retraining changes the coefficients). Later: the model
  was retrained on GENCODE 50 (A8), and the bundled example data, the mini
  annotation, the knockdown example and the case-study bundles were rebuilt
  on GENCODE 50 (mouse vM39) with new Salmon quantifications (f3a1a87; see
  `NEWS.md`).
- The benchmark table in `ARCHITECTURE.md` is from GENCODE v33; it was not
  recomputed on v50. Later: recomputed on v50 (A7).

## 7. Backlog: what is relevant and implementable

Status on 2026-11-18: A1, A2, A3, A4, A7 and A8 are done (A3, A4 and A7
rerun or rechecked on 2026-11-17 after A8), A6 is closed, A9 is partly done
and A5 is open. Effort: S = under half a day, M = one to two
days, L = more.

| ID | Item | Effort | Blocked by | Kind |
|---|---|---|---|---|
| A1 | `annot_seq_regions.parquet` keeps 1 of 706 sequence regions | S | done | defect |
| A2 | Frozen source cache is keyed by file name, not by version | S | done | defect |
| A3 | Pilot (done 2026-11-12; rerun 2026-11-17: no gain on Karousis, small gain on retained introns) | M | done | benchmark |
| A4 | Splice-site strength from the annotation itself (no AlphaGenome) | M | done | feature |
| A5 | Intron-retention evidence for retained-intron transcripts | M-L | A3 or coverage data | feature |
| A6 | `add_alphagenome_splicing()` evidence tier (closed on judgement, see A6) | L | closed | feature |
| A7 | Recompute the rules benchmark on GENCODE 50 | S | done | docs |
| A8 | Retrain the coding potential model on GENCODE 50 (done) | S-M | done | data |
| A9 | Upload the rebuilt APPRIS files (done); other datasets, Zenodo (open) | S | partly done | release |

**A1. Sequence regions. Done.** GENCODE GFF3 files put one
`##sequence-region` line before the features of each sequence (706 in v50
`chr_patch_hapl_scaff`), not all in the header, and only the header was
read: the annotation directories of v49, v50 and vM38 held one row (`chr1`).
Fix: `read_sequence_regions()` in `R/annotation.R` scans the whole file with
DuckDB; `prepare_isoformic_annotation()` merges the declared lengths with
the per-sequence largest feature end that the SQL already wrote. Cost: about
10 s on a full GENCODE file. Test in `tests/testthat/test-annotation.R`. The
three cached annotation directories on this machine were patched in place
(706, 706 and 61 rows); other existing directories need
`overwrite = TRUE`.

**A2. Source cache. Done.** `frozen_sources()` downloaded into
`get_isoformic_cache("sources", <name>)` and skipped files that exist; after
a version change in `frozen_registry()` the file name is the same
(`appris_data.appris.txt`), so the old source was rebuilt under the new name
with the old MD5. Fix: the default directory is
`get_isoformic_cache("sources", <name>, <version>)`. Test in
`tests/testthat/test-frozen.R`. The local cache was moved to that layout
(both APPRIS versions are kept, `e113v49` and `e116v50`).

**A3. Pilot. Done 2026-11-12; rerun 2026-11-17** (section 4, step 3). The
first run (one seed) was read as negative. The rerun (five seeds, bootstrap
over genes, rebuilt targets) gives no detectable gain in pooled AUC on
Karousis (baseline + AlphaGenome minus baseline -0.0010, 95% interval
-0.0158 to 0.0133) and a small gain on the retained-intron biotype (+0.0069,
0.0044 to 0.0096; +0.0026, 0.0006 to 0.0047, beyond the contained-intron
flag). Within Karousis pairs AlphaGenome orders more pairs right than the
PWM (165 against 137 of 223, McNemar p = 0.0036); a within-pair model test
is still open. The plan
was: steps in section 4. Decide the conditions after
`--list-conditions`: about four, chosen for having observed (not imputed)
`SPLICE_JUNCTIONS` and `RNA_SEQ` tracks. Budget: 500 genes x conditions requests of
1 Mb. Then add the AlphaGenome features to `03_evaluate.R` (its baseline is
in section 4) and record the result here and in `ARCHITECTURE.md`.

**A4. Done on 2026-11-12.** `R/splice-sites.R`: `extract_splice_sites()`,
`has_splice_sites()`, `splice_site_strength()`; table
`annot_data_splice_sites.parquet`, also written by
`extract_transcript_sequences(genome_fasta = )`. Not added to
`splicing_events()` / `poison_exons()`: it stays a separate table. Result
on GENCODE 50 (698,612 unique introns, matrices trained on 260,664 introns
of canonical transcripts; first computed by a scratch script, recomputed on
2026-11-17 by `data-raw/alphagenome/splice_site_groups.R` from the splice
sites of the cached annotation directory, with
`ISOFORMIC_SPLICE_EXTRACT=false`; every number below reproduced, and
`R/splice-sites.R` has had only documentation changes since it was added):

| Group | n | Median weakest donor | Median weakest acceptor |
|---|---|---|---|
| `protein_coding` transcripts | 276,953 | 4.72 | 4.08 |
| `nonsense_mediated_decay` | 91,818 | 4.43 | 3.64 |
| `retained_intron` | 33,567 | 6.25 | 6.13 |
| `lncRNA` | 184,403 | 6.23 | 5.95 |

The weakest site falls with the number of introns (median 8 to 9 for
coding and NMD, 2 to 3 for the others), so the groups are not comparable
without that. Per intron, on the main chromosomes: introns of canonical
transcripts have a median donor of 8.49 and acceptor of 8.10; introns seen
only in NMD transcripts 7.98 and 7.17. AUC of canonical against NMD-only
introns: donor 0.558, acceptor 0.587, the smaller of the two 0.616. AUC of
protein coding against NMD transcripts by weakest site: 0.548 (donor) and
0.550 (acceptor). Conclusion:
the first-order score is real but weak evidence, and the matrices were
trained on the canonical introns, which favours them. This is the baseline
that AlphaGenome must beat in A3; `03_evaluate.R` now holds it (numbers in
section 4, step 3). The v50
annotation directory in the local cache now holds the table.

Original note. `ARCHITECTURE.md` listed
splice-site strength as not implemented because position weight matrices are
missing. They can be built from the annotation: with the genome FASTA that
`extract_transcript_sequences()` already streams, take the donor (3 exonic +
6 intronic bases) and acceptor (20 intronic + 3 exonic) windows of the
junctions of canonical transcripts, build log-odds matrices against the
background, and score every junction. Where: next to
`sequences_from_genome()` in `R/sequences.R` (windows are extracted while a
chromosome is in memory) and a new column pair in `splicing_events()` /
`poison_exons()`. Gives `weakest_donor`, `weakest_acceptor` per transcript
at no external cost, and a baseline that AlphaGenome must beat in A3: add it
to the step 3 comparison. A first-order matrix is weaker than MaxEntScan;
say so in the docs.

**A5. Intron retention.** The rules call 478 of 1,781 retained-intron
transcripts `altered_productive` (the projected ORF ends inside a last
retained intron, so no downstream junction and no NMD by the 50-nt rule).
Evidence that the intron is retained in mature, exported RNA would settle
them; neither the rules nor the annotation have it. Candidates: AlphaGenome
intronic over exonic coverage (A3), or measured intron coverage from the
user's data, which the package does not read. 463 with the GENCODE 49
coding model. A3 (2026-11-12, rerun 2026-11-17): the predicted ratio barely
separates retained-intron transcripts (AUC 0.554 with seed 1, five-seed mean
0.558, 95% interval 0.533 to 0.584), so only measured intron coverage is
left; open.

**A6. Evidence tier. Closed 2026-11-12.** Section 5. The condition was
that A3 improves the `retained_intron` or `protein_coding` rows over the
rules plus A4. The first run was read as "it does not". The rerun of
2026-11-17 shows a small improvement of the `retained_intron` row with an
interval that excludes zero, so the condition is met for that row; the
`protein_coding` row was not evaluated (no truth for it). A6 stays closed
as a judgement, to be confirmed by the maintainer: the gain is small (0.007
AUC), the truth is a biotype label, not productivity, the model may have
seen these genes in training, and the non-commercial terms and the open
question of redistributing predictions (section 3) stand in the way.
Reopening it would need the within-pair test on Karousis and the GTEx
isoform-fraction control first (section 4, step 3).

**A7. Benchmark table. Done; rerun 2026-11-17.** `ARCHITECTURE.md`
"Benchmark of the sequence rules" holds the full GENCODE v50 run (646,577
transcripts with sequence) next to the v33 values: NMD biotype 0.998,
retained intron 0.588, protein coding 0.909, canonical 0.973, lncRNA 1.000.
The first v50 run (2026-11-11, GENCODE 49 coding model, `benchmark_v50.rds`)
came from an inline `Rscript -e` call; the rerun with the GENCODE 50 model
uses the tracked script `data-raw/benchmark_rules.R` (its header says how to
run it) and is kept in `benchmark_v50_2026-11-18.rds` (with the commit,
4a6e983, and the coefficients of the coding model). Only the retained-intron
row moved (0.598 to 0.588: unproductive 20,491 to 20,144, altered productive
9,254 to 9,424, uncertain 3,549 to 3,726); outside the table
`protein_coding_CDS_not_defined` (altered 4,908 to 4,970, uncertain 1,749 to
1,734, unproductive 19,836 to 19,789) and `TEC` (3 transcripts from
unproductive to uncertain) moved too. The `protein_coding` and
`protein_coding_LoF` break-downs are unchanged. The coding model is the
only input that changed; that it explains the moves is shown on the pilot
transcripts (step 1) and inferred for the full annotation, which was not
rerun with the old model. The rerun reused the ORF table cached in the
annotation directory on 2026-11-13: its key matches the current model and
package version, and the only code change in `R/orfs.R` since then
restructures an `if`. Classification took 49 s.
`protein_coding_LoF` (148 transcripts) was first left out of
`unproductive_biotypes()` and was added on 2026-11-12 (decision of the
maintainer): the calls describe the reference genome, in which a variant
stops the translation. The benchmark uses `use_biotype = FALSE`, so its
numbers do not change.

**A8. Coding model. Done 2026-11-12.** `data-raw/build_coding_model.R`
has `release <- "50"`; the shipped model is retrained (9,938 coding, 10,000
non-coding). Held-out GENCODE 50 transcripts (4,969 and 5,000): sensitivity
0.941, specificity 0.919, AUC 0.978; the v49 model on the same set: 0.928,
0.936, 0.977 (part of the set was in its training sample). Equivalent
models. The default test suite passes. Rerun with the new model on
2026-11-17: the A7 benchmark table and the pilot targets and evaluation
(section 4). No article, `README.md` or `NEWS.md` states productivity
counts on full GENCODE: the productivity article computes
`benchmark_productivity()` when it is rendered with `ISOFORMIC_FULL=true`.

**A9. Release. Done for APPRIS (2026-11-12).** The maintainer created the
`frozen-data` release of `luciorq/isoformic` (it did not exist before) with
`appris_e116v50.parquet` and `appris_domains_e116v50.parquet`; the
downloaded files match the MD5 of `frozen_registry()`. The other frozen
datasets are not in the release, so `fetch_frozen_data()` still rebuilds
them from their sources; no Zenodo deposit yet.

Redistribution terms of the frozen datasets, checked on 2026-11-12 (reading
of the published terms, not legal advice):

| Dataset | Terms found | Source of the terms | Upload |
|---|---|---|---|
| `nmd_exons` | Article CC BY 4.0; CC0 waiver for "the data made available in this article" | Permissions block of the full text (Europe PMC, PMC8361881) | yes, with citation |
| `gtex_transcript_tpm` | GTEx Portal Data License: redistribution allowed; acknowledge the Portal as source with the date; state that v8 is not the most current release; no implied endorsement | `gtexportal.org/home/license` (text taken from the portal application script) | yes, with these statements |
| `riboseq_orfs` | No licence on the files or READMEs. GENCODE: "All the GENCODE project data is open access". EMBL-EBI: "places no additional restrictions on the use or redistribution of the data", expects attribution | `gencodegenes.org/pages/data_access.html`, `ebi.ac.uk/about/terms-of-use` | probably; not explicit |
| `polyasite` | No licence on the home, about or atlas pages. Article is CC BY-NC (permissions link of the NAR page); the workflow code is Apache-2.0 | `polyasite.unibas.ch`, GitHub API | unknown; ask the Zavolan lab |
| `appris`, `appris_domains` | No licence on the data server. The code repository is CC BY-NC-ND 4.0 ("APPRIS by Jose Manuel Rodriguez is licensed under ..."); if that covers the data, a reformatted table is a derivative and may not be distributed | `github.com/appris/appris` `LICENSE.md`; `apprisws.bioinfo.cnio.es/pub/README.md` | uploaded on 2026-11-12 before this check; ask the authors or remove |

Not read: the APPRIS papers, the Ensembl legal page, the Mudge et al. 2022
paper.
