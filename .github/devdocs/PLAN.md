# Plan: next steps for isoform-level analysis

Started 2026-10-23. Status values: `todo`, `doing`, `done`, `blocked`.
Update the table and tick the boxes as work lands; add a line to the log at
the bottom for every decision or surprise.

## Overview

| ID | Workstream | Deliverable | Status |
|----|------------|-------------|--------|
| A | Differential transcript usage | `test_usage()` with `diffsplice`, `satuRn` and `DEXSeq` backends on the uncertainty-scaled counts, benchmarked on the read-split data | done |
| B | Terminus groups | `read_terminus_groups()` and custom groups in `add_uncertainty()` / `de_confidence()`, compared with the correlation groups | done |
| C | Long-read replicates (Oarfish) | `read_oarfish()` and Parquet replicates in `read_inferential_replicates()`, validated on real Oarfish output | done |
| D | IsoformSwitchAnalyzeR v2 comparison | Switch calls of `isoform_switches()` against ISAR v2 on GSE129540, written up | done |
| E | Wrap-up | Docs, NEWS, article, ARCHITECTURE, full verification | done |
| F | Usage tests for more than two groups | `test_usage()` with several levels: each level against the reference or all pairs, one fit per backend where it allows | done |
| H | Omnibus usage test | `contrasts = "omnibus"`: DEXSeq likelihood-ratio test over all groups; Simes over the contrasts against the reference for diffSplice and satuRn | done |
| I | Oarfish single-cell | Read the single-cell output of Oarfish (real output from a simulated tagged BAM) and pseudobulk it into an IsoformicExperiment | done |
| J | Test everything | Coverage of the new code, the unrun script branches run for real, examples, edge cases | done |
| K | TreeTerminus | `read_treeterminus()` and `cut_transcript_trees()`: adaptive groups from transcript trees, cut where inferential variance falls below a cutoff | done |
| L | Stability article | 30 vs 100 Gibbs draws and the molecular features that make an isoform hard to quantify, as a new article | done |
| G | Terminus from R | `run_terminus()` (group per sample, consensus collapse) through condathis, with input checks, validated against the shell run | done |

Order: A, B, C, D, E. Each workstream ends with its own verification before
the next starts.

## Working rules

- Every new exported function: roxygen topic, tests, entry in `_pkgdown.yml`,
  NEWS bullet.
- Dependencies stay in `Suggests` and are called as `pkg::fun()`; tests skip
  when a backend is missing.
- Validate on real output of the tool being read (Terminus, Oarfish), not on
  a format guessed from documentation.
- Heavy files live in `/data/gamma/luciorq/workspaces/temp/isoformic-knockdown/`.
- Nothing is committed or uploaded without being asked.

## A. Differential transcript usage

Why: benchmarks (2024-2025) favour DEXSeq, limma/edgeR `diffSplice` and
satuRn; Baldoni et al. (2025) show that dividing out the read-ambiguity
overdispersion keeps the error rate of usage tests.

- [x] A1 Install limma, satuRn, DEXSeq; prototype each on the GSE129540 data
- [x] A2 `test_usage(self, group, levels, method, scale = TRUE)` returning
      transcript and gene level results (`gene_id`, `transcript_id`,
      `usage_log2FC`, `pvalue`, `qvalue`, `gene_pvalue`, `gene_qvalue`,
      `method`, `scaled`); `scale = NULL` scales for diffsplice when draws exist
- [x] A3 Tests: synthetic known switch, missing backend, group levels
- [x] A4 Benchmark on the read-split null contrasts and real contrast:
      false positives, replication across halves, per backend, with and
      without the overdispersion scaling
- [x] A5 Document results in the uncertainty article and `data-raw/`

Acceptance: all three backends run on the bundled data; the benchmark table
is bundled; results reported even when they are negative.

## B. Terminus groups

Why: Terminus builds data-driven groups of transcripts that the reads cannot
separate; `collapse_groups` does the same with a correlation rule.

- [x] B1 Install Terminus; run on the GSE129540 Salmon output; record the
      real output layout
- [x] B2 `read_terminus_groups(path)` mapping transcripts to groups
- [x] B3 `add_uncertainty(groups = <table>)` (or equivalent) so
      `de_confidence(collapse_groups = TRUE)` can use external groups
- [x] B4 Tests with a fixture written from real Terminus output
- [x] B5 Compare Terminus groups with the correlation groups (overlap, size,
      effect on the confidence classes)

Acceptance: groups read from real Terminus output; comparison table written.

## C. Long-read replicates (Oarfish)

Why: Oarfish leads long-read quantification benchmarks and can write
inferential replicates; Parquet fits the package.

- [x] C1 Install Oarfish, minimap2, samtools; simulate a small long-read
      dataset from the bundled transcript FASTA; run Oarfish with replicates
      to obtain the real file layout
- [x] C2 `read_oarfish(path)` for quantification, ambiguity counts and run
      metadata
- [x] C3 Parquet replicates accepted by `read_inferential_replicates()`
      (matrix per replicate, same return value as for Salmon)
- [x] C4 Tests with a tiny fixture generated from real Oarfish output
- [x] C5 Document in the preprocessing article

Acceptance: fixture comes from real Oarfish output; `add_uncertainty()` runs
on it unchanged.

## D. IsoformSwitchAnalyzeR v2 comparison

- [x] D1 Install IsoformSwitchAnalyzeR (v2 line); import the GSE129540 Salmon
      quantifications with the GENCODE v49 annotation
- [x] D2 Run its switch test; list the significant switches
- [x] D3 Compare with `isoform_switches()` and the confidence classes
      (agreement, switches only one tool calls, confidence of each)
- [x] D4 Write up in `data-raw/` and the uncertainty article, with the
      versions used

Acceptance: a table of agreement and disagreement; no claim about accuracy
without ground truth.

## E. Wrap-up

- [x] E1 `air format`, `jarl check .`, `lintr::lint_package()`, ASCII scan
- [x] E2 `devtools::document()`, full tests, `pkgdown::check_pkgdown()`,
      all articles, R CMD check
- [x] E3 NEWS, README, ARCHITECTURE, `data-raw/README.md`
- [x] E4 Stage everything; final report with the open items

## F. Usage tests for more than two groups

- [x] F1 `group_setup()` generalising `two_group_setup()` (two or more levels,
      samples outside `levels` dropped)
- [x] F2 `contrasts` argument: `NULL` (each level against the first),
      `"pairwise"`, or a list of `c(numerator, denominator)`; results stacked
      with a `contrast` column, Simes and BH per contrast
- [x] F3 Backends: diffSplice with one `~0 + group` fit and a contrast
      vector per comparison; satuRn with one fit and a contrast matrix;
      DEXSeq on the samples of each pair
- [x] F4 Tests: three-group synthetic switch, pairwise, explicit contrasts,
      two-group behaviour unchanged apart from the `contrast` column
- [x] F5 Docs, NEWS, article sentence

Acceptance: a three-group fixture finds the switches in the right contrasts
with all three backends.

## G. Terminus from R

- [x] G1 `create_terminus_env()` and `run_terminus(paths, output_dir,
      consensus_thresh, threads, keep_logs, env_name)` through condathis
- [x] G2 Input checks: replicates and equivalence-class weights in every
      Salmon directory (a class line has 2k + 2 fields with weights, k + 2
      without; `cmd_info.json` does not record the flag)
- [x] G3 Returns `read_terminus_groups()` of the consensus; removes the
      multi-gigabyte `delta.log` unless `keep_logs = TRUE`
- [x] G4 Tests with mocked `condathis::run()` (arguments, order, checks)
- [x] G5 Real run on the six GSE129540 samples; the consensus must equal the
      shell run (md5 of `clusters.txt`)
- [x] G6 Docs, NEWS, article chunk

Acceptance: identical consensus to the shell run.

## H. Omnibus usage test

- [x] H1 `contrasts = "omnibus"` in `test_usage()`; one row per transcript,
      `contrast = "omnibus"`, `usage_log2FC = NA`
- [x] H2 DEXSeq: full against reduced model on all groups; diffSplice and
      satuRn: Simes combination of the k - 1 contrasts against the reference
      (valid under positive dependence; any difference among groups implies
      a difference from the reference)
- [x] H3 Tests: three-group fixture, both switch sets found; two groups give
      the same p-values as the single contrast for diffSplice
- [x] H4 Docs

## I. Oarfish single-cell

- [x] I1 Real output: simulated reads with cell barcodes, minimap2 to the
      transcriptome, `CB:Z` tags, `oarfish --single-cell --alignments`
- [x] I2 Reader for the real layout (sparse counts, barcodes, transcripts)
- [x] I3 Pseudobulk by a cell-to-sample table into an IsoformicExperiment
- [x] I4 Tests from a tiny real fixture; docs

## J. Test everything

- [x] J1 Coverage of R/usage.R, R/terminus.R, R/oarfish.R, R/uncertainty*.R;
      tests for uncovered branches
- [x] J2 Run the computation branch of `validate_uncertainty.R` for real and
      compare with the bundled tables
- [x] J3 Make `quantify_knockdown.sh` testable at small scale (number of read
      pairs, run filter) and run it end to end
- [x] J4 Full verification: tests, lint, articles, R CMD check

## K. TreeTerminus

- [x] K1 `read_treeterminus(path, salmon_dir | transcripts)`: Newick trees of
      `group`/`consensus` (leaves are 0-based indices in `quant.sf` order,
      checked against the real output), parent-pointer table
- [x] K2 `cut_transcript_trees(trees, x, cutoff = 0.5, unresolved)`: node InfRV
      from the summed draws (bottom-up per tree), lowest ancestor below the
      cutoff for each transcript, maximal nodes kept so that groups are disjoint
- [x] K3 Groups source "treeterminus" in `add_uncertainty()` and
      `uncertainty_methods()`
- [x] K4 Tests: parser on real output lines, cut on the uncertainty fixture
- [x] K5 Validation on the six prostate samples against Terminus and the
      correlation groups; docs. The conda recipe is the user's.

## L. What makes an isoform hard to quantify

- [x] L1 100 Gibbs draws for three GSE129540 samples; TreeTerminus per sample at
      30 and 100 draws
- [x] L2 Features (data-raw/knockdown/isoform_features.R part A): length, exons,
      isoforms in gene, GC, MANE/basic/TSL, NAGNAG and near alternative
      acceptors/donors, unique exonic bases, closest sibling, containment
- [x] L3 Outcomes (part B): InfRV, InfRV disagreement and cutoff flips between
      independent 30-draw sets, tree change 30 vs 100, unique fragments
- [x] L4 Summaries and models (part C), bundled under inst/extdata/knockdown/stability
- [x] L5 Article isoformic-quantification-stability.qmd, pkgdown, NEWS; checks

## Log

- 2026-10-23: plan created. State before starting: R CMD check OK, all tests
  pass, uncertainty features (add_uncertainty, de_confidence,
  inferential_overdispersion, scale_counts) and the read-split validation are
  in the staged tree.
- 2026-10-23: A1-A3 done. All three backends run (diffSplice 5 s, satuRn 53 s,
  DEXSeq ~10 s per 800 genes on 40k transcripts). satuRn 1.18 has no
  `parameterization` argument; DEXSeq fold changes need
  `estimateExonFoldChanges()` then `DEXSeqResults()`. A4 benchmark running
  (`data-raw/knockdown/benchmark_usage.R`, log in the work directory).
- 2026-10-23: C done. Real Oarfish 0.10.3 layout: `<prefix>.quant` (tname, len,
  num_reads), `.infreps.pq` (columns `bootstrap.0`..., same row order),
  `.ambig_info.tsv` (unique_reads, ambig_reads, total_reads, no names),
  `.meta_info.json`. Oarfish maps reads itself (`--reads --annotated`), so
  minimap2 was not needed. Reads are simulated, not a real library.
- 2026-10-23: B in progress. Terminus 0.1.0 needs `--dumpEqWeights` and
  replicates in the Salmon output; the earlier runs lacked the weights, so one
  sample was re-quantified for the format check.
- 2026-10-23: A done. Read-split benchmark (technical replicates): unscaled
  diffSplice 41/59 and DEXSeq 195/244 false transcript calls at q<0.05 in the
  two null contrasts, satuRn 0; scaled 7/2, 12/4 and 0/2. Real-contrast calls
  of half 1 go the same way in half 2: unscaled diffSplice 47% (chance),
  DEXSeq 65%; scaled 100% and 94%. Scaling therefore defaults to on whenever
  replicates exist. Caveat: scaled runs test ~3x fewer transcripts (expression
  filter on divided counts).
- 2026-10-24: D done. Installed IsoformSwitchAnalyzeR is 2.10.0 (Bioconductor
  release here), not a newer 2.x. It needs the `chr_patch_hapl_scaff` GTF because
  the Salmon index includes patch transcripts. Like-for-like (q<0.05,
  |dIF|>=0.1): ISAR 40 calls of 79,463 tested; test_usage DEXSeq raw 372 of
  19,171 (8 shared with ISAR); scaled 139 of 8,779 (1 shared). Comparison is
  descriptive: different tested sets, no ground truth.
- 2026-10-24: B code done (read_terminus_groups, add_uncertainty(groups=)).
  Real Terminus output: `terminus group -d <salmon dir> -o <prefix>` then
  `terminus collapse -d <salmon dirs> -o <prefix>` writes `<prefix>/<run>/
  clusters.txt` (`NewTr<id>,<tx>,<tx>`). Pending: six-sample Terminus run and
  the comparison with the correlation groups (B5).
- 2026-10-24: B done. Six-sample Terminus 0.1.0 consensus: 14,029 groups, 52,911
  transcripts (a single sample gave 17,008 groups); identical `clusters.txt` in
  every sample folder. Against the correlation groups: 67% of correlation-grouped
  transcripts are also Terminus-grouped, 65% of correlation pairs share a
  Terminus group; median InfRV 15.5 (Terminus), 8.6 (correlation), 1.4 (neither);
  among |log2FC|>1 changes high stays 203, ~900 move moderate to low. Two bugs
  found on the way: a local variable shadowing the new `groups` argument, and
  `meta$uncertainty_group` partially matching the new metadata key (renamed to
  `uncertainty_ambiguity_source`; `$` on lists partial-matches).
- 2026-10-24: E done. 159 test blocks pass (1 skip: writexl), jarl and lintr
  clean, pkgdown check clean, all articles build, R CMD check Status OK (needs
  `_R_CHECK_FORCE_SUGGESTS_=false` because writexl is not installed). The usage
  backends (edgeR, limma, satuRn, DEXSeq, BiocParallel) are in Suggests.
  Everything is staged; nothing is committed.
- 2026-10-27: F and G added at the user's request. diffSpliceDGE takes a
  `contrast` vector, so one multi-group fit serves all comparisons.
  salmon_quant() already writes `--dumpEqWeights`.
- 2026-10-27: F done. `contrasts` = NULL (each level vs first), "pairwise" or
  pairs; `contrast` column added to the output. Three-group fixture (switches
  in b for g1-g8, in c for g9-g16) passes for all three backends. Two-group
  results on the prostate data are identical to the earlier `~group` fit.
  `group_setup()` replaces the body of `two_group_setup()`.
- 2026-10-27: G1-G4 and G6 done (run_terminus, input checks, mocked tests,
  docs). G5 running: run_terminus(workers = 6) on the six quants_terminus
  directories into terminus_r/. While linting, jarl 0.7.1 reported 18 findings
  in files that had passed before (undesirable_operator `<<-` in test mocks and
  closures, cyclomatic_complexity in resolve_ensembl_file, nzchar, library_call,
  an unused object in an example); the binary and jarl.toml are unchanged and
  the cause is unknown. All were fixed (environments instead of `<<-`,
  `ensembl_file_pattern()` helper) and the affected tests pass.
- 2026-10-27: G done. run_terminus(workers = 6) on the six GSE129540 samples took
  21 minutes and wrote a consensus identical to the command-line run (all 12
  clusters.txt share md5 038a7bf2...); logs removed, 427 MB left. 163 test
  blocks pass, lint clean, articles build, R CMD check OK. Staged, not committed.
- 2026-10-29: H, I, J added ("go for it"). edgeR's diffSplice.DGEGLM and
  diffSpliceDGE test one coefficient or contrast only, hence the Simes omnibus
  for diffSplice and satuRn.
- 2026-10-29: H done. Three-group omnibus finds both switch sets with every
  backend; with two groups the omnibus p-values equal the single contrast
  (diffSplice, DEXSeq).
- 2026-10-29: I done. Real `oarfish --single-cell` layout: `<prefix>.count.mtx`
  (MatrixMarket coordinate real, cells as rows), `.barcodes.txt`,
  `.features.txt`, `.meta_info.json`; no bootstraps. 60 cells, 9,000 reads, all
  assigned. The boosted RBPJ-203 rises 4.6-fold in the treated pseudobulk (the
  other two boosted isoforms have too few reads to show it).
- 2026-10-29: J2 done. validate_uncertainty.R computation branch rerun from
  scratch: the four bundled tables are identical; per-transcript results are
  identical except swish q-values (up to 0.0097), which come from random
  permutations. set.seed() added before swish in validate_uncertainty.R,
  build_knockdown_example.R and the preprocess article.
- 2026-10-29: J3 done. quantify_knockdown.sh gained N_PAIRS and RUNS. Running it
  end to end (100k pairs, one run per series) exposed a real bug: under
  `set -euo pipefail` the `pigz -dc | head` pipeline always fails and would
  have aborted the first download; fixed (the record check decides). Mapping
  93.0% and 86.0%, as in the full runs. Both Oarfish fixture scripts also run
  end to end: read counts reproduce to 1e-13 and single-cell counts exactly;
  Oarfish bootstraps are unseeded and barcode order varies (documented).
- 2026-10-30: J done. Coverage 93% overall; tests added for the uncovered
  branches (Oarfish path forms in the shared readers, inconsistent and
  incomplete input, Matrix Market errors, parallel Terminus failure, env
  creation); R/terminus.R now 100% and R/oarfish.R 99.7% from their own tests.
  176 test blocks pass, lint clean, all ten articles rebuilt, R CMD check OK.
  Stack imbalance warnings during pkgdown article builds: deterministic (three
  per build), traced to the Bioconductor package SparseArray 1.10.10 (Bioc
  3.22): `asNamespace("SparseArray")` alone prints a stack imbalance warning;
  none of its dependencies do. The articles reach it through downlit, which
  loads the namespaces of attached packages (`library(fishpond)`,
  `library(SummarizedExperiment)`) to link function names. Not isoformic code
  (no compiled code) and outputs are unaffected. Builds are reproducible apart
  from the per-session temporary paths printed in four articles.
  Correction: an earlier entry claimed six clean builds; those builds had
  silently failed ("No built file found") because a test script left a
  `vignettes/_quarto.yaml` behind (on.exit() does nothing at the top level of an
  Rscript script). The file was removed and the builds redone.
- 2026-10-31: K started at the user's request ("improvement against terminus
  alone for some hard cases"); the user maintains a conda recipe for
  TreeTerminus. Built here from source (commit 6a81966) with GSL from
  conda-forge.
- 2026-10-31: K done. On the six prostate samples (116,229 assessed isoforms):
  grouped 8,076 / 20,703 / 31,606 (correlation / Terminus / TreeTerminus cut),
  uncertain left alone 78,892 / 67,138 / 58,225, features after collapse
  112,258 / 103,692 / 90,326, median group InfRV 1.48 / 2.87 / 0.19. Reading and
  cutting the trees: 218 s, 4.4 GB. Leaf positions follow quant.sf, verified (80%
  of trees within one gene against 13% with the eq_classes order). Bug class
  seen three times in the comparison scripts: a tibble() column shadowing a
  variable used later in the same call; package code scanned, only identity
  cases there.
- 2026-10-31: L done. 3 samples, 87,430 isoforms with >= 5 reads (85,236 in the
  models). log InfRV: R2 0.63, led by unique_fraction (40% of explained
  deviance), reads (29%), MANE (11%); NAGNAG and near splice sites act through
  unique fragments (<0.1% once those are in). Disagreement between 30-draw sets:
  R2 0.08, mostly sampling noise. Tree change: type (62%) and isoforms in gene
  (26%); lncRNA trees change about half the time, also for uncertain isoforms.
  Article isoformic-quantification-stability.qmd with tables in
  inst/extdata/knockdown/stability/.
