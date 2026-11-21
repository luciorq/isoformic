# Data processing: reproduce and update

Everything under `inst/extdata/` (except the coding model inputs downloaded
by its script), the data of the articles under `vignettes/articles/data/`
(except the complete Reactome GMT file) and every frozen evidence
dataset is built by code in this directory or in `R/frozen-data.R`. All
examples use GENCODE release 50 (human) and M39 (mouse); the shell scripts
take the release from `GENCODE_RELEASE`. Run the scripts from the
package root with `Rscript`. Downloads go to the isoformic cache
(`get_isoformic_cache()`, usually `~/.cache/R/isoformic/`), so a second run
is offline.

| Script | Produces | Sources |
|---|---|---|
| `placenta/`, `prepare_extdata.R` | `inst/extdata/example_experiment/` (with transcript sequences), `example.gff3.gz`, `example_transcripts.fa.gz`, `mini_quants/`, `mini_draws.parquet`, `mini_annotation/`, `reactome_subset.gmt`; `vignettes/articles/data/preeclampsia/pe_*.parquet` | Salmon output of PRJNA622166 (below), GENCODE v50 GFF3 and transcript FASTA, `vignettes/articles/data/preeclampsia/c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt` (MSigDB 2023.1) |
| `build_coding_model.R` | `inst/extdata/coding_model/hexamer_logodds.parquet`, `coding_model_coefficients.parquet` | GENCODE v50 GFF3, `pc_transcripts.fa.gz` and `transcripts.fa.gz` |
| `build_frozen_data.R` | Parquet files of the frozen datasets in `get_isoformic_cache("frozen")` | see `frozen_datasets()` |
| `build_knockdown_example.R`, `knockdown/` | `inst/extdata/knockdown/`, `knockdown_annotation/` | Salmon output of GSE129540 and GSE75489 (below) |
| `oarfish/` | `inst/extdata/oarfish/`, `oarfish_sc/` | simulated reads (below) |
| `case_studies/` | `vignettes/articles/data/<study>/` | NCBI SRA reads (`fetch_sra.sh`, `quantify_*.sh`; usage in the header of `build_case_studies.R`) |
| `uncertainty_evidence/` | `vignettes/articles/data/uncertainty_evidence/` (the `evidence_numbers`, `flt1_isoforms`, `replication_flags` and `story_genes` Parquet tables of `isoformic-uncertainty-evidence.qmd`, copied by `91_bundle_article_data.R` from the tables that `90_article_tables.R` writes; the larger outputs stay outside the repository) | Salmon output of PRJNA622166, PRJNA984605, GSE129540 and GSE75489 from NCBI SRA; pre-registration in `prereg.yaml` and `prereg_amendments.yaml` (run order in `data-raw/uncertainty_evidence/README.md`) |
| `alphagenome/` | nothing in the package | AlphaGenome pilot (`.github/devdocs/ALPHAGENOME-BENCHMARK.md`) |
| `benchmark_memory.R`, `benchmark_rules.R` | nothing in the package | memory benchmark of the draws and genome-wide annotation code paths (`just benchmark-memory`, full mode) and benchmark of the productivity sequence rules (see Benchmarks below; headers give the variables) |

`.gitignore` keeps the scripts (`*.R`, `*.sh`, `*.py`), the sample sheets
(`*.tsv`), this file, the pixi files of `alphagenome/` and the README and
pre-registration YAML files of `uncertainty_evidence/` in git, and ignores
everything else under `data-raw/` (downloads, work files, environments).

## Bundled example data

```sh
data-raw/knockdown/quantify_knockdown.sh /path/to/workdir   # builds the index
data-raw/placenta/quantify_placenta.sh /path/to/workdir
ISOFORMIC_PLACENTA_DIR=/path/to/workdir Rscript data-raw/prepare_extdata.R
```

The example data come from the preeclampsia placenta dataset (PRJNA622166):
8 placentas of early-onset severe preeclampsia and 17 normal placentas
(`placenta/placenta_samples.tsv`). `quantify_placenta.sh` downloads every run
from NCBI SRA and quantifies it with Salmon 2.8.0 and 30 Gibbs draws against
the decoy-aware GENCODE index built by `quantify_knockdown.sh` in the same
work directory (pass another index with `INDEX`). Three preeclampsia and
three control runs are also quantified, with 20 Gibbs draws, against an
index of the transcripts of the 15 genes of `placenta/mini_genes.tsv` only:
these are the mini quantifications.

`prepare_extdata.R` then runs `swish` (fishpond) on the transcripts and on
the genes (draws summed over the transcripts of each gene) of the complete
quantifications and writes `pe_*.parquet` (TPM, results of every tested
transcript and gene, sample table) to `vignettes/articles/data/preeclampsia/`:
only the articles read them, and at about 12 MB they would put the package
over the 5 MB size limit of CRAN. The example experiment keeps the
showcase genes and the 60 genes with the lowest q-value, with the TPM and the
Salmon counts as assays and the transcript sequences in its annotation;
`example.gff3.gz` and `example_transcripts.fa.gz` hold three of its genes so
that the annotation and sequence examples run offline. The mini
quantifications are copied without their draws, which go to
`mini_draws.parquet` (read it with `as_inferential_replicates()`), and with
the work directory paths replaced by `/work`.

To move the examples to a newer GENCODE release, set `GENCODE_RELEASE` for
the shell scripts and `release` in `prepare_extdata.R`, rerun them and update
the numbers that `README.md`, the tests and the articles print for
`example_isoformic()`.

## Coding potential model

```sh
Rscript data-raw/build_coding_model.R   # or: just coding-model
```

Trains the Fickett + hexamer logistic model used by `coding_potential()`
and `transcript_orfs()` (`coding_probability`). Training data are a
balanced random sample (seed 1, 15,000 per class) of GENCODE human
`protein_coding` transcripts with a complete CDS (positive class; the
hexamer table is built from their CDS) and `lncRNA` transcripts (negative
class). One third is held out and the script prints the held-out
sensitivity and specificity at the 0.364 cutoff.

Current model: GENCODE release 50, 9,938 coding and 10,000 non-coding
training transcripts, held-out sensitivity 0.941, specificity 0.919 (4,969
coding and 5,000 non-coding transcripts, AUC 0.978). The release 49 model
scores 0.928 and 0.936 on the same transcripts (AUC 0.977; some of them were
in its training sample): the two are equivalent, the new one trades
specificity for sensitivity at the cutoff.

To update: set `release` in the script to the new GENCODE release, rerun,
record the numbers above in this file and in `NEWS.md`. For another
organism, point `gff`, `fasta` and `pc_fasta` to that organism's files
(the GENCODE mouse files work unchanged with `organism = "mouse"`).

## Frozen evidence datasets

The registry in `R/frozen-data.R` (`frozen_registry()`) pins, for every
dataset, the source URLs, a version label, the output file name and the MD5
of the frozen Parquet file. `fetch_frozen_data()` downloads the file from
the GitHub release tagged `frozen-data` of this repository and verifies the
MD5; when the release asset is missing it rebuilds the file from the sources
with the `build_frozen_*()` function of the dataset.

| Dataset | Source version | Update when |
|---|---|---|
| `appris`, `appris_domains` | APPRIS `2026_08.v50` release, `e116v50` (Ensembl 116 / GENCODE 50; the `e<N>` part is the Ensembl release, the `v<N>` part the APPRIS annotation version) | a new APPRIS release matches the GENCODE release you annotate with (`https://apprisws.bioinfo.cnio.es/pub/releases/`) |
| `riboseq_orfs` | GENCODE Ribo-seq ORFs phase 1 (`data/Ribo-seq_ORFs.bed`, 26-column bigGenePred) and phase 2 comprehensive BED12 | GENCODE publishes a new phase (`https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/`) |
| `nmd_exons` | Karousis et al. 2021 Genome Biology, Supplementary Table 1 (`13059_2021_2439_MOESM2_ESM.xlsx`) | a newer validated NMD target catalog is adopted |
| `polyasite` | PolyASite 2.0 atlas, GRCh38.96 clusters BED | PolyASite releases a new atlas |
| `gtex_transcript_tpm` | GTEx v8 RSEM transcript TPM (`GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz`, 3.6 GB) and sample attributes | GTEx releases new transcript-level data |

### Rebuild all datasets

```sh
Rscript data-raw/build_frozen_data.R            # all datasets (just frozen-data)
Rscript data-raw/build_frozen_data.R appris     # one dataset (just frozen-data appris)
```

The script downloads the sources into `get_isoformic_cache("sources", <name>, <version>)`
(the GTEx file takes a few minutes to download and about six minutes to
reduce to tissue medians), writes the Parquet files into
`get_isoformic_cache("frozen")` and prints the MD5 of every file.
`build_frozen_nmd_exons()` needs the `readxl` package.

### Publish a new version

1. Edit `frozen_registry()` in `R/frozen-data.R`: update `version` and the
   `sources` URLs of the dataset (for APPRIS, change the release path and the
   `e<ensembl>v<gencode>` directory; for GENCODE Ribo-seq ORFs, the phase
   directories).
2. Run `Rscript data-raw/build_frozen_data.R <name>` and check the row counts
   and a few rows of the output (`read_frozen_data(<name>)`).
3. Copy the printed MD5 into the `md5` field of the dataset in
   `frozen_registry()`.
4. Upload the Parquet file as an asset of the `frozen-data` GitHub release
   (create the release once; `--clobber` replaces an existing asset):

   ```sh
   gh release create frozen-data --title "Frozen evidence datasets" --notes "See data-raw/README.md"
   gh release upload frozen-data ~/.cache/R/isoformic/frozen/<file>.parquet --clobber
   ```

   Keep the previous file name when the version label changes, or add the
   new file and leave the old one so that older package versions keep
   working.
5. Deposit the snapshot on Zenodo (below) so that the files stay available
   under a DOI independently of GitHub.
6. Run `just test` and `just check`, add a `NEWS.md` entry naming the new
   source version, and bump the package version.

### Zenodo deposit of the processed files

```sh
Rscript data-raw/zenodo_deposit.R            # or: just zenodo-deposit
```

`bundle_frozen_data()` copies every frozen Parquet file and the coding
potential model into one directory (default
`get_isoformic_cache("zenodo_deposit")`) with an `MD5SUMS` file and a README
listing dataset versions and sources. Upload that directory as a new version
of the Zenodo deposit, keeping the file names, publish it and set the record
id in `frozen_zenodo_record()` in `R/frozen-data.R`. From then on
`fetch_frozen_data()` tries, in order: any base URL in
`options(isoformic.frozen_mirrors = )`, the `frozen-data` GitHub release and
the Zenodo record (`https://zenodo.org/records/<id>/files/<file>?download=1`),
and only rebuilds from sources when none has the file. Cite the versioned DOI
in `NEWS.md` and `inst/CITATION`. To reproduce an analysis with a specific
deposit, pin it with
`options(isoformic.frozen_mirrors = "https://zenodo.org/records/<id>/files/")`;
the MD5 in `frozen_datasets()` must match the deposited file, so a deposit
and the package version that lists its checksums go together.

### Column contracts

The analysis functions depend on these columns; keep them when updating a
builder:

- `appris`: `transcript_id` (no version), `appris_annotation`,
  `appris_principal`, `appris_score`, the method score columns.
- `appris_domains`: `transcript_id`, `seqid`, `start_pos`, `end_pos`,
  `strand`, `pfam_name`.
- `riboseq_orfs`: `orf_id`, `seqid`, `strand`, `start_pos`, `end_pos`,
  `n_blocks`, `block_starts`, `block_ends` (1-based, `;`-separated),
  `orf_type`, `transcript_ids` (`;`-separated, no version, may be `NA`).
- `nmd_exons`: `seqid`, `start_pos`, `end_pos`, `strand`,
  `nmd_transcript_ids`.
- `polyasite`: `seqid`, `cleavage_pos`, `strand`, `site_id`, `site_class`,
  `tpm`.
- `gtex_transcript_tpm`: `transcript_id` (no version), `gene_id`, `tissue`,
  `median_tpm`, `n_samples`.

Identifiers are stored without version (`strip_id_version()`, Ensembl and
RefSeq patterns only) and chromosome names with the `chr` prefix; the join
functions add or remove the prefix to match the annotation.

## Benchmarks

Neither script changes the package; run them on a shared machine under a
memory cap, and record their output outside the repository.

- `benchmark_memory.R` (`just benchmark-memory`, `pixi run benchmark-memory`)
  runs each step in its own R process and writes `benchmark_memory.parquet`
  (peak RSS above the start of the step, time, and the peak as a multiple
  of one copy of the draws) plus one Parquet file per step. The draws steps
  read Salmon output with Gibbs draws (`ISOFORMIC_BENCHMARK_QUANTS`), the
  annotation steps use GENCODE (`ISOFORMIC_FULL=true`). To compare a change,
  run it once on the old code (`ISOFORMIC_BENCHMARK_PKG` pointing at a copy
  of the package before the change) and once on the new one, with the same
  samples and draws (`ISOFORMIC_BENCHMARK_SAMPLES`, `_DRAWS`) and each with
  its own copy of the annotation directory (`ISOFORMIC_BENCHMARK_ANNOT`).
  `ISOFORMIC_BENCHMARK_STEPS` reruns a subset of the steps. A fresh copy of
  an annotation directory has no valid ORF table, so the first
  `productive_fraction` step also computes the ORFs of the whole annotation
  (about 4 minutes and 6 GB on GENCODE 50): rerun that step to measure it
  with the ORFs stored. Steps that peak below about 0.3 GB vary from run to
  run by as much as they differ between versions.
- `benchmark_rules.R` writes an RDS file with the agreement of the sequence
  rules of `transcript_productivity()` with the GENCODE biotypes over a full
  annotation; the table "Benchmark of the sequence rules" in
  `.github/devdocs/ARCHITECTURE.md` comes from it. Rerun it when the rules,
  the ORF code or the coding potential model change.

## ESRP knockdown examples

`inst/extdata/knockdown/<series>/` and `inst/extdata/knockdown_annotation/`
hold two real datasets in which an isoform switch defines the phenotype, used
by the uncertainty article: GSE129540 (prostate cancer cells, `ESRP2`
siRNA; Munkley et al. 2019) and GSE75489 (H358 cells, `ESRP1/2` shRNA; Yang
et al. 2016), three control and three knockdown samples each. They are
rebuilt in two steps; the Salmon output is not part of the package.

```sh
data-raw/knockdown/quantify_knockdown.sh /path/to/workdir   # reads, index, Salmon
ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/build_knockdown_example.R
```

The first step downloads each run from NCBI SRA and keeps its first 12
million read pairs, builds a decoy-aware Salmon 2.8.0 index of GENCODE v50
and quantifies every run with 30 Gibbs draws and the equivalence classes
(about 60 GB of disk space, 20 GB of it the index and 1.1 GB the
quantifications). The second runs `swish` on the
complete quantifications, so the bundled q-values are those of the full
analysis, and keeps the transcripts of 39 genes (19 known ESRP targets plus
the 12 genes with the most significant transcripts of each series): 2,001
annotated transcripts of 52 gene IDs, because `HLA-B`, `HLA-C`, `IFI27` and
`SCRIB` also have copies on alternative haplotype and patch sequences. Each series
folder holds `abundance`, `counts`, `det`, `col_data` (with the Salmon run
metadata and the `library_size` of the complete quantification),
`equivalence` and `draws` (long tables, see `inferential_replicates_table()`) and a `README.txt`. The two
series folders and the annotation take about 1.1 MB.

The Salmon directories of the work directory are the larger dataset: they are
not shipped, and can be deposited on Zenodo (`data-raw/zenodo_deposit.R`
shows the procedure used for the frozen datasets) if the quantifications
should be citable.

## Read-split validation, usage tests and tool comparisons

These scripts check the uncertainty tools on GSE129540 and write small result
tables to `inst/extdata/knockdown/read_split/`. They need the work directory of
`quantify_knockdown.sh` (`ISOFORMIC_KNOCKDOWN_DIR`); heavy files stay there.

| Script | Needs | Writes |
|---|---|---|
| `knockdown/split_halves.sh WORKDIR` | `sub/` reads | `quants_half/` (12 Salmon runs on half of the reads) |
| `knockdown/validate_uncertainty.R` | `quants_half/` | `null_calls`, `null_stable`, `replication_flag`, `replication_rv` |
| `knockdown/benchmark_usage.R` | `quants_half/` | `usage_null`, `usage_replication` |
| `knockdown/compare_isar.R` | `quants_sub/`, the GENCODE v50 `chr_patch_hapl_scaff` GTF, IsoformSwitchAnalyzeR | `isar_comparison` |
| `knockdown/compare_terminus.R` | Terminus run on all six samples (below) | `terminus_comparison` |
| `knockdown/compare_treeterminus.R` | Terminus run and TreeTerminus Mean trees (below) | `grouping_comparison` |

The Terminus runs need Salmon output with `--dumpEqWeights`, which the main
quantification does not write. `knockdown/gibbs_terminus.sh WORKDIR` makes
them: the six GSE129540 runs again with the weights (`quants_terminus/`),
three of them with 100 Gibbs draws (`g100/`), `terminus group` and
`collapse` on the six runs (`terminus_all/out/`), TreeTerminus per sample at
30 and 100 draws (`tt_single/`) and on the six runs (`tt_mean/`, about 18 GB
of memory). TreeTerminus has no Bioconda package yet; build it from source
with Rust and GSL and pass its path in `TREETERMINUS`. `terminus group`
writes a 3.8 GB `delta.log` per sample, deleted by the script.

`knockdown/isoform_features.R` writes `inst/extdata/knockdown/stability/`
(the tables of the quantification stability article) from the GENCODE v50
annotation and three GSE129540 samples quantified with 100 and 30 Gibbs
draws, with TreeTerminus run per sample (`gibbs_terminus.sh`); its header lists the
folders of the work directory it needs.

## Oarfish examples

`inst/extdata/oarfish_sc/` holds real Oarfish 0.10.3 single-cell output
(`out.count.mtx`, `out.barcodes.txt`, `out.features.txt`,
`out.meta_info.json`) for 60 simulated cells in three samples, and the table of
cells (`cells.tsv`). Rebuild with
`data-raw/oarfish/run_oarfish_sc.sh /path/to/workdir`
(`simulate_single_cell.R` makes the reads; minimap2 aligns them and the
barcode is moved to the `CB:Z` tag).

`inst/extdata/oarfish/` holds real Oarfish 0.10.3 output (`.quant`,
`.ambig_info.tsv`, `.meta_info.json`, `.infreps.pq` with 30 bootstrap
replicates) for three samples of simulated nanopore reads. Rebuild with
`data-raw/oarfish/run_oarfish.sh /path/to/workdir` from the package root
(`simulate_long_reads.R` makes the reads). The reads of both examples are
simulated from the 141 GENCODE v50 transcripts of
`inst/extdata/example_transcripts.fa.gz`; they only provide a realistic file
layout.

Rerunning both scripts reproduces the read counts (to 1e-13) and the
single-cell counts exactly, but not byte for byte: Oarfish 0.10.3 draws its
bootstrap replicates without a fixed seed (it has no seed option), and writes
the barcodes in the order its threads finish. The readers do not depend on
that order.
