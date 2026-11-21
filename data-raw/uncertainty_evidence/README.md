# Uncertainty evidence: does quantification uncertainty matter?

These scripts run a pre-registered, genome-wide evaluation of whether
transcript quantification uncertainty (Salmon Gibbs draws, InfRV, the
read-to-transcript ambiguity overdispersion of Baldoni et al. 2024)
changes differential isoform results, and whether the uncertainty flags of
isoformic predict which calls fail to replicate. Negative results are part
of the evaluation and are reported as such.

Four work packages:

- **WP-A**, placenta null: control-versus-control splits of 15 placentas
  (3 v 3, 5 v 5, 7 v 7), where every call is false, plus a gene-thinning
  power arm.
- **WP-B**, placenta audit: the real preeclampsia versus control comparison:
  how many calls are fragile, held-out replication over 30 half-splits used
  in both directions, and the FLT1, ENG and LEP story.
- **WP-C**, an independent mouse null (10 control males, all 126 5 v 5
  splits) and random technical read halves of two knockdown series.
- **WP-D**, the overdispersion prior estimated genome-wide against the prior
  of the bundled 39-gene knockdown panel.

## Pre-registration

`01_preregister.R` froze the package code (a `git archive` snapshot of one
commit, with its sha256) and wrote `prereg.yaml` (seeds, methods,
thresholds, universe, strata, decision rules R_A and R_B, story genes, stop
rules) before any test ran. Every later script loads the snapshot, never the
live `R/`, and stamps each table with the git sha, the snapshot sha256, the
seed and the script name. The one exception is the `productivity` stage of
`30_audit_prevalence.R` (amendment A18): after `productive_fraction_sd` and
`delta_z` changed definition, it was rerun with `ISOFORMIC_UE_PKG=.`, the
repository code, and its two tables are stamped `<sha>+working-tree`.

`prereg.yaml` is frozen. The committed copy differs from the original only
in its absolute paths, replaced by `${ISOFORMIC_UE_OUT}`-style placeholders
(amendment A10 gives the sha256 of the original). Every change to the plan
is in `prereg_amendments.yaml`; the ones made or recorded after the results
were seen say so (`post_hoc: true`). The mouse null
(WP-C) is not in the plan: it applies the R_A rule of the placenta null to a
dataset the plan does not name (A19).

## Data

All reads are public and come from NCBI SRA:

| Dataset | Accessions | Used for |
|---|---|---|
| Preeclampsia placenta, 8 preeclampsia and 17 controls | GSE148241, PRJNA622166 (SRR11498039 to SRR11498046, SRR11498063 to SRR11498079) | WP-A, WP-B, WP-D |
| Mouse adBNST, chronic stress, 10 control and 9 stressed males | PRJNA984605 | WP-C, WP-D |
| ESRP2 knockdown | GSE129540 (SRR8870780 to SRR8870785) | WP-C, WP-D |
| ESRP1/2 knockdown, H358 | GSE75489 (SRR2966412 to SRR2966417) | WP-C, WP-D |

`00_inventory.R` also lists GSE139962 and the prostate long-read runs
(GSE259334) when they are present; the evaluation does not test them.

The quantifications (Salmon, decoy-aware GENCODE 50 index, 30 Gibbs draws;
GENCODE vM39 for the mouse) come from the scripts of `data-raw/placenta/`,
`data-raw/case_studies/` and `data-raw/knockdown/`. The scripts expect them
in one work directory, `ISOFORMIC_V50_DIR`:

```
human/ref/            GENCODE 50 GTF, transcripts.fa.gz, salmon_index_v50/
human/envs/salmon/    Salmon used by 43_rhalf_quant.sh
human/placenta/quants/<run>/
human/quants_sub/<GSE>_<run>/    12 M read pairs (quantify_knockdown.sh)
human/quants_half/<run>_h1|h2/   consecutive halves (split_halves.sh)
human/g100/<run>/                100 draws (gibbs_terminus.sh)
cases/mouse_bnst/quants/<run>/, cases/mouse_bnst/ref/
```

The annotation directory of isoformic for GENCODE 50 (from
`prepare_annotation()`; default in the isoformic cache, or
`ISOFORMIC_UE_ANNOT`) is read for the productivity classes and the story.

## Environment variables

| Variable | Meaning |
|---|---|
| `ISOFORMIC_UE_OUT` | output root; required by every script |
| `ISOFORMIC_V50_DIR` | work directory above; required by the scripts that read quantifications |
| `ISOFORMIC_KNOCKDOWN_DIR` | work directory of the GSE139962 runs (`00_inventory.R`) and default parent of `KD_SUB_DIR` (`42_split_reads.sh`) |
| `ISOFORMIC_UE_THREADS` | DuckDB and data.table threads (default 4) |
| `ISOFORMIC_UE_RLIB` | optional extra R library (for example with edgeR and limma) put first on `.libPaths()` by the overdispersion-prior scripts (`50_` to `59_`); default `$ISOFORMIC_V50_DIR/rlib`, used only when that directory exists |
| `ISOFORMIC_UE_PKG` | package code to load (default: the frozen snapshot of `prereg.yaml`; `.` loads the repository, amendment A18) |
| `ISOFORMIC_UE_ANNOT` | isoformic GENCODE 50 annotation directory (default: the one in the isoformic cache) |

A script stops with a message when a variable it needs is unset. Each
script header lists its own options (workers, stages, bootstrap sizes).

## Running

Run every script from the package root, for example
`pixi run --frozen Rscript data-raw/uncertainty_evidence/20_null_splits.R`.
Each header gives the measured (or estimated) run time and peak memory. Run
the heavy steps one at a time inside a memory cap a little above that peak,
for example `systemd-run --user --scope -p MemoryMax=8G <command>`, and set
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for the scripts with forked
workers.

Run order:

1. `00_inventory.R`, `01_preregister.R`
2. S0: `10_placenta_prepare.R` (the longest step, about 4 h)
3. WP-A: `20_null_splits.R`, `21_null_dte.R` (variants main, member,
   draws16, unfolded), `22_null_dtu.R`, `23_thin_dte.R`,
   `29_null_summarise.R`
4. WP-B: `30_audit_prevalence.R` (stages), `31_halfsplit_replication.R`
   (stages), `32_story.R`, `39_audit_summarise.R`, `33_post_hoc_checks.R`
   (post hoc: A9, and the same-size controls of A21)
5. WP-C: `40_mouse_prepare.R`, `41_mouse_null.R` (fit, then analyse),
   `42_split_reads.sh`, `43_rhalf_quant.sh`, `44_rhalf_replicate.R`,
   `49_replication_summarise.R`
6. WP-D: `50_od_prior.R`, then `51` to `54`, then `59_od_summarise.R`
7. `34_cluster_bootstrap.R` (post hoc, amendment A11; reads the outputs of
   `31_halfsplit_replication.R` and `44_rhalf_replicate.R`) and
   `92_package_flag_check.R` (post hoc, amendment A24: recomputes the
   sign agreement of every discovery half with the repository code, by
   default `ISOFORMIC_UE_PKG=.`, and counts its disagreements with the
   evaluated flag), then `90_article_tables.R`, then `91_bundle_article_data.R` (copies the
   article tables into `vignettes/articles/data/uncertainty_evidence/`)

`ue_common.R`, `null_common.R`, `38_audit_common.R`, `lib_null_methods.R`,
`30b_prevalence_summary.R` and `31b_replication_summary.R` are sourced by
the others, not run.

## Results

Everything is written under `ISOFORMIC_UE_OUT` as Parquet: `store/` and
`S0/` (S0), `null_placenta/` (WP-A), `placenta_audit/` (WP-B),
`null_mouse/`, `rhalf/` and `wpc_summary/` (WP-C), `od_prior/` (WP-D).
`90_article_tables.R` copies every number the article states, with its 95%
CI and its source table, into `article_tables/`; those small tables are the
ones the article reads. The article
`vignettes/articles/isoformic-uncertainty-evidence.qmd` reads a copy of
them, written by `91_bundle_article_data.R`, in `vignettes/articles/data/uncertainty_evidence/`, with
`evidence_numbers.parquet` pruned to the keys the article uses (see the
header of that article).
