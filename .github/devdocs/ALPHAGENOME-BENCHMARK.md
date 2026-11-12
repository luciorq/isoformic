# AlphaGenome pilot benchmark and GENCODE 50 bump: notes

Written 2026-11-11. Status: targets built, scoring blocked on an API key,
evaluation not written.

**Where to start.** Section 4 "To resume" has the three commands that
continue the pilot once a key exists. Section 7 is the backlog (A1 to A9):
what is implementable now without AlphaGenome, what waits for the pilot, and
two defects found on the way (A1, A2; fixed, with A7, on 2026-11-11).

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
  (listed as not implemented in `ARCHITECTURE.md`), intron retention
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

Files: `data-raw/alphagenome/` (git-ignored by `data-raw/*`, like the other
`data-raw` sub-directories). Own pixi project (`pixi.toml`, `pixi.lock`;
python 3.12, pyarrow, pandas, alphagenome from PyPI). Work directory:
`ISOFORMIC_ALPHAGENOME_DIR` (used:
`/data/gamma/luciorq/workspaces/temp/isoformic-alphagenome`).

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

Result: 500 genes, 18,378 transcripts, 19,226 unique introns, 252 NMD and 257
non-NMD labelled transcripts.

Rules-only baseline on the sample (`use_biotype = FALSE`,
`use_evidence = FALSE`):

| Biotype | n | unproductive | productive | altered_productive | uncertain |
|---|---|---|---|---|---|
| `nonsense_mediated_decay` | 5,042 | 5,031 | 0 | 10 | 1 |
| `retained_intron` | 1,781 | 1,136 | 3 | 463 | 179 |
| `protein_coding` | 10,439 | 40 | 2,309 | 7,011 | 1,079 |
| `protein_coding_CDS_not_defined` | 1,106 | 820 | 2 | 216 | 68 |

So: NMD biotype recall 0.998, retained intron 0.64, protein coding called
productive or altered 0.89. Room for AlphaGenome is in the retained-intron
class and the 463 retained-intron transcripts called `altered_productive`.

### Step 2, `02_score.py` (running since 2026-11-12)

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

### Step 3, `03_evaluate.R` (baseline written and run; AlphaGenome part not written)

The script holds the baseline that the AlphaGenome features have to beat:
the calls of the rules and the splice-site strength of
`splice_site_strength()` (backlog A4), with the truths, the grouped
cross-validation and the paired test already in place. AlphaGenome features
go in as further entries of `feature_sets`. Writes
`evaluation/{features,auc,pairs}.parquet` in the work directory.

Baseline on the v50 targets (one run, seed 1, logistic regression, 5 folds
grouped by gene, pooled out-of-fold AUC):

| Truth | n (positive) | Rules | PWM | Intron number | Rules + PWM |
|---|---|---|---|---|---|
| Karousis NMD against non-NMD | 509 (252) | 0.784 | 0.729 | 0.578 | 0.805 |
| NMD against protein coding biotype | 15,478 (5,042) | 0.996 | 0.599 | 0.517 | 0.997 |
| Retained intron against protein coding biotype | 12,188 (1,752) | 0.870 | 0.840 | 0.796 | 0.963 |

- PWM features: weakest donor and acceptor over all introns, the same over
  the introns that the canonical transcript of the gene does not have
  ("novel"), and whether there is a novel intron.
- Karousis pairs (259, same gene): the NMD transcript has the weaker novel
  donor in 142 of 233 pairs that differ (61%, sign test p = 0.001) and the
  weaker novel acceptor in 149 of 244 (61%, p = 0.0007). Over all introns the
  direction is the opposite for the donor (35 of 107, 33%): the restriction
  to novel introns is what carries the signal.
- Reading. On the only truth that is independent of the annotation
  (Karousis) the sequence score adds about 0.02 AUC to the rules. The NMD
  biotype row is circular for the rules (the biotype follows the same
  50 nt rule) and leaves no room. The retained-intron row is the one with
  room (0.870 to 0.963), but intron number alone gives 0.796 there:
  retained-intron models are short fragments with few introns, so most of
  that gain is transcript structure, not splice-site sequence.
- Consequence for the stop rule below: AlphaGenome has to beat 0.805 on
  Karousis and 0.963 on retained introns with the intron-number control in
  the model, not the rules alone.

Still planned for the AlphaGenome part:

- Features per transcript: weakest annotated junction (predicted count
  relative to the gene maximum, zero when not predicted), weakest splice-site
  probability and usage, the same restricted to junctions absent from the
  canonical isoform, intron-retention ratio (intronic over flanking exonic
  coverage) for exon segments that contain a canonical intron.
- Diagnostics first: match rate of predicted junctions to annotated introns
  under offsets of -1, 0, +1 on each end; offset of the splice-site peak at
  canonical sites.
- Comparisons: rules alone, AlphaGenome alone, both (logistic regression,
  cross-validation grouped by gene).
- Truths: Karousis pairs (paired within gene), biotype, and the
  retained-intron class.
- Control: GTEx isoform fraction (`gtex_transcript_tpm`) as a feature, to
  tell splicing information from an expression proxy.
- Stop rule: if the `retained_intron` and `protein_coding` rows do not
  improve, do not integrate.

### To resume

```sh
export ALPHAGENOME_API_KEY=...        # not on this machine on 2026-11-11
cd data-raw/alphagenome
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --list-conditions
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --smoke
pixi run python 02_score.py --work-dir "$ISOFORMIC_ALPHAGENOME_DIR" --conditions <curies>
```

`test_v49/` in the work directory is a disposable trial of step 1 on v49.

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
  `release <- "49"`; retraining changes the coefficients).
- The benchmark table in `ARCHITECTURE.md` is from GENCODE v33; it was not
  recomputed on v50.

## 7. Backlog: what is relevant and implementable

Status on 2026-11-12: A1, A2, A4 and A7 are done; the others are open. Effort: S = under half a
day, M = one to two days, L = more.

| ID | Item | Effort | Blocked by | Kind |
|---|---|---|---|---|
| A1 | `annot_seq_regions.parquet` keeps 1 of 706 sequence regions | S | done | defect |
| A2 | Frozen source cache is keyed by file name, not by version | S | done | defect |
| A3 | Finish the pilot: smoke test, conditions, scoring, AlphaGenome part of step 3 (baseline done) | M | API key | benchmark |
| A4 | Splice-site strength from the annotation itself (no AlphaGenome) | M | done | feature |
| A5 | Intron-retention evidence for retained-intron transcripts | M-L | A3 or coverage data | feature |
| A6 | `add_alphagenome_splicing()` evidence tier | L | A3 positive, licence | feature |
| A7 | Recompute the rules benchmark on GENCODE 50 | S | done | docs |
| A8 | Retrain the coding potential model on GENCODE 50 | S-M | decision | data |
| A9 | Upload the rebuilt APPRIS files; Zenodo | S | decision | release |

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

**A3. Finish the pilot.** Steps in section 4. Decide the conditions after
`--list-conditions`: about four, chosen for having observed (not imputed)
`SPLICE_JUNCTIONS` and `RNA_SEQ` tracks. Budget: 500 genes x conditions requests of
1 Mb. Then add the AlphaGenome features to `03_evaluate.R` (its baseline is
in section 4) and record the result here and in `ARCHITECTURE.md`.

**A4. Done on 2026-11-12.** `R/splice-sites.R`: `extract_splice_sites()`,
`has_splice_sites()`, `splice_site_strength()`; table
`annot_data_splice_sites.parquet`, also written by
`extract_transcript_sequences(genome_fasta = )`. Not added to
`splicing_events()` / `poison_exons()`: it stays a separate table. Result
on GENCODE 50 (one run, scratch script, 698,612 unique introns, matrices
trained on 260,664 introns of canonical transcripts):

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
introns: donor 0.56, acceptor 0.59, the smaller of the two 0.62. AUC of
protein coding against NMD transcripts by weakest site: 0.55. Conclusion:
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

**A5. Intron retention.** The rules call 463 of 1,781 retained-intron
transcripts `altered_productive` (the projected ORF ends inside a last
retained intron, so no downstream junction and no NMD by the 50-nt rule).
Evidence that the intron is retained in mature, exported RNA would settle
them; neither the rules nor the annotation have it. Candidates: AlphaGenome
intronic over exonic coverage (A3), or measured intron coverage from the
user's data, which the package does not read. Keep as an open question until
A3 shows whether the predicted ratio separates anything.

**A6. Evidence tier.** Section 5. Only if A3 improves the `retained_intron`
or `protein_coding` rows over the rules plus A4.

**A7. Benchmark table. Done.** `ARCHITECTURE.md` "Benchmark of the sequence
rules" now holds the full GENCODE v50 run (646,577 transcripts with
sequence) next to the v33 values: NMD biotype 0.998, retained intron 0.598,
protein coding 0.909, canonical 0.973, lncRNA 1.000.
`protein_coding_LoF` (148 transcripts) was first left out of
`unproductive_biotypes()` and was added on 2026-11-12 (decision of the
maintainer): the calls describe the reference genome, in which a variant
stops the translation. The benchmark uses `use_biotype = FALSE`, so its
numbers do not change. Result kept in `benchmark_v50.rds` in the work
directory.

**A8. Coding model.** `data-raw/build_coding_model.R` has `release <- "49"`;
`just coding-model` retrains. Coefficients and the 0.364 cutoff may move;
tests and article numbers that depend on them must be rechecked. No evidence
that v49 training is a problem.

**A9. Release.** `appris_e116v50.parquet` and
`appris_domains_e116v50.parquet` are in the local cache only. Until they are
uploaded to the `frozen-data` release, `fetch_frozen_data()` rebuilds them
from the APPRIS server (procedure in `data-raw/README.md`). The files are 2.5
times larger than before (8.1 and 7.3 MB).
