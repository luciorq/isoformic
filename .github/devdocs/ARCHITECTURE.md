# isoformic architecture notes

Last updated: 2026-11-18 (first written 2026-10-15 around
`IsoformicExperiment`).

## Layers

| Layer | Files | Notes |
|-------|-------|-------|
| DuckDB / Parquet IO | `R/utils-duckdb.R` | `duckdb_connect()` (self-closing), `duckdb_run()`, `query_registered()`, exported `read_parquet()`, `write_parquet()`, `read_tabular()`, and `read_parquet_matching()` (rows whose key is in a set of values: the values are registered as a table and joined with a SEMI JOIN on `IS NOT DISTINCT FROM`, so `NA` matches `NA` as with `%in%`; file order is restored in R from `file_row_number`, not with `ORDER BY`, which would sort the sequence strings). The only parquet engine; `arrow` is not used. Every connection gets `duckdb_limits_config()`: the `isoformic.duckdb_threads` and `isoformic.duckdb_memory_limit` options, or defaults from Slurm and the CPU affinity mask (capped by cgroup v2 limits), which DuckDB does not detect itself. |
| Annotation directory | `R/annotation.R` | `prepare_isoformic_annotation()` parses GFF3/GTF with DuckDB into `annot_data_{genes,transcripts,exons}.parquet`, `annot_metadata.parquet`, `annot_seq_regions.parquet` (`build_annotation_sql()`): one read of the file into a narrow `features` table (exon and CDS rows keep only transcript, exon ID and number; other rows keep the attribute string), then the attribute extractions on the gene and transcript rows only. The raw regex values are extracted in an inner query and the placeholders removed in the outer one: inside a `CASE`, DuckDB evaluates both arguments of each nested `nullif()`, so cleaning there runs every regex 32 times. `read_annotation_table(filter = list(<column> = values))` reads only the matching rows through `read_parquet_matching()`; `exon_table_from_dir()` reads genes, then the transcripts of the selected genes, then their exons. The parser is structural (genes = parents of transcripts, transcripts = parents of exons, plus type lists) so GENCODE, Ensembl, NCBI RefSeq and T2T Liftoff/CAT files all work; `source_gene_id`/`source_transcript_id` keep lifted identifiers. `resolve_annot_path()` accepts a directory or a file. |
| Readers | `R/read_salmon.R`, `R/make_tx_to_gene.R`, `R/download_reference.R`, `R/salmon_utils.R` | `read_salmon()` pivots `quant.sf` files in DuckDB; `make_tx_to_gene()` parses FASTA headers in DuckDB; `download_reference()` defaults to GENCODE 50 (mouse M39), returns a cached file without contacting the server when the request determines its name, and lists the server otherwise (`list_reference_dir()`, with `cached_fallback()` to one matching cached file when the listing fails); `salmon_index()` / `salmon_quant()` run Salmon through `condathis`. |
| Cache and options | `R/get_isoformic_cache.R`, `R/isoformic-options.R` | `get_isoformic_cache()` is the user cache of downloads and frozen data; `?"isoformic-options"` lists every option and environment variable (see `CLAUDE.md`). |
| Class | `R/IsoformicExperiment-class.R` | S7 class, custom constructor, validator (works on stored attributes, never on computed properties), `print`, `dim`, `dimnames`, `[`. Two caches: package-level `annotation_cache` keyed by file signature (path, size, mtime) and a per-object `.cache` environment keyed by hashes of the inputs. |
| Accessors | `R/IsoformicExperiment-methods.R` | SummarizedExperiment-style generics. File name sorts after `-class.R` on purpose (collation). |
| DE results | `R/dea.R` | `standardize_de_table()` harmonises tool-specific column names; cutoffs live in `dea$cutoffs`. Table helpers `join_DEG_DET()`, `is_deg_sig()`. |
| Diversity | `R/IsoformicExperiment-diversity.R` | `rowsum()` / `order()` + `!duplicated()` tricks instead of per-group closures; `exon_diversity()` aggregates in DuckDB. |
| Persistence | `R/IsoformicExperiment-io.R` | `write_isoformic()` / `read_isoformic()`, `as_isoformic()`, `as_summarized_experiment()`, `example_isoformic()`. `store_annotation()` copies the annotation into the directory with `copy_annotation = TRUE` or when it lies inside the directory being replaced, and drops the copied `annot_data_orfs.parquet`. The draws in `metadata$inf_reps` go to `inf_reps/` (`write_inf_reps()`/`read_inf_reps()`: one wide Parquet file per draw, `draw_0001.parquet`, ..., the moments in `inf_reps/moments/`, and `inf_reps_n_draws`, `inf_reps_samp_type`, `inf_reps_moments_n` in `metadata.parquet`), one draw at a time; a directory without `inf_reps/` reads without draws. |
| Plots | `R/plot_log2fc.R`, `R/profile.R`, `R/plot_genomic_context.R`, `R/enrichment.R`, `R/uncertainty-plots.R`, `R/tx_type_palette.R` | Pure ggplot2; shared `isoformic_theme()` and `complete_palette()`. |
| Transcript features | `R/transcript-features.R` | `transcript_features()` (ORF/UTR/NMD rule/TSS-TES groups/canonical comparison), `splicing_events()` (compares each transcript with the canonical isoform of its gene), `add_annotation_columns()` evidence hook, `add_fantom_cat_tables()`. `compute_splicing_events()` joins on integer keys (`tid`, the first row of the transcript; `ctid`, the row of the gene's canonical transcript; strand as 1, -1 or 0), computes exon rank, count and transcript span in R, and runs one materialized overlap join: novel and skipped exons are anti joins on it. The character methods of `transcript_features()` and `transcript_orfs()` read only the genes asked for (`annotation_of_genes()`). |
| Productivity | `R/productivity.R` | `transcript_productivity()`, `productive_fraction()`, `run_productivity_enrichment()`, `poison_exons()`, `productivity_palette()`. `productive_fraction()` keeps the genes of `gene_name` and `min_gene_count` before summing the assay and the draws. `poison_exons()` uses the loaded exon table, or reads the needed columns of the object's transcripts (`exons_of_transcripts()`); `poison_exon_candidates()` counts with `tabulate()`. |
| Annotation tools | `R/annotation-tools.R` | `merge_annotations()`, `add_sqanti_classification()`. `merge_annotations()` never loads the tables into R: `merge_annotation_table()` runs one DuckDB `COPY (SELECT ... UNION ALL BY NAME ...)` per table and relies on insertion order (no `ORDER BY`). The sequences table is copied with one DuckDB thread and `ROW_GROUP_SIZE 20000`: an order-preserving parallel `COPY` of long strings buffers out-of-order batches (6 GB on two GENCODE 50 annotations, 0.5 GB with one thread). |

## Conventions

- Every dependency call is `pkg::fun()` (see `CLAUDE.md`).
- S7 generics are declared with explicit signatures (`S7::new_generic(name, "self", function(self, ...) S7::S7_dispatch())`) so R CMD check's codoc matches the documentation; methods are not given their own roxygen block.
- Computed properties (`annot_data_*`, `annot_metadata`, `gene_assays`) read raw values with `attr(self, "<prop>", exact = TRUE)` inside getters and the validator to avoid recursion and unwanted computation.
- Bundled data is built on GENCODE 50 (mouse M39) by the scripts of
  `data-raw/` (excluded from the build by `.Rbuildignore`):
  `data-raw/placenta/quantify_placenta.sh` quantifies the preeclampsia
  placenta runs (PRJNA622166, `placenta_samples.tsv`) with Salmon 2.8.0, and
  `data-raw/prepare_extdata.R` turns that output and the cached GENCODE
  references into `example_experiment/`, the `pe_*.parquet` tables (in
  `vignettes/articles/data/preeclampsia/`, with the complete Reactome GMT:
  only the articles read them, and they would put the tarball over the
  5 MB CRAN limit), `example.gff3.gz`, `example_transcripts.fa.gz`, `mini_quants/` (genes of
  `mini_genes.tsv`), `mini_draws.parquet`, `mini_annotation/` and
  `reactome_subset.gmt`. `data-raw/knockdown/` (Salmon) with
  `data-raw/build_knockdown_example.R` builds the ESRP knockdown examples,
  `data-raw/oarfish/` simulates reads from `example_transcripts.fa.gz` and
  runs Oarfish 0.10.3 for `oarfish/` and `oarfish_sc/`,
  `data-raw/case_studies/` builds the case-study bundles, and
  `data-raw/alphagenome/` holds the AlphaGenome benchmark
  (`ALPHAGENOME-BENCHMARK.md`). The coding model and the frozen datasets
  have their own scripts (`data-raw/README.md`).
- Sample groups: every function with `group` and `levels` gets its samples
  from `group_samples()` (`R/IsoformicExperiment-diversity.R`; the pure core
  is `group_factor()`), and subsets its matrices to them before any
  computation. The result must equal the one of `x[, samples]`: samples
  outside `levels` and samples without a group value are dropped, a
  requested level without samples is an error, factor levels without
  samples are not groups. Tests in `tests/testthat/test-group-levels.R`.
  Exception: `add_uncertainty()` and `de_confidence()` restrict the contrast
  only; their uncertainty tables describe all the samples of the object.
- `[` subsets by ID everything stored by transcript or sample (assays,
  `col_data`, stored gene assays, `metadata$inf_reps`, the rows of the
  uncertainty tables). It does not recompute the uncertainty tables when
  samples are dropped; it informs with class `isoformic_stale_uncertainty`.
  The draws are a view: `subset_inf_reps()` keeps the parent's matrices and
  sets the `selection` attribute (`transcripts`, `samples`) of
  `metadata$inf_reps`, which a nested subset replaces. Every consumer reads
  them through `inf_reps_reader()`, `inf_reps_draw_index()` (used by
  `draws_view()`) or `inf_reps_moments()`; `write_inf_reps()` writes the
  selected rows and columns, and `inf_reps_materialize()` gives the plain
  copy. The parent's draws stay alive while a subset exists; to compare two
  groups of a larger object, pass `levels =` instead of subsetting. A subset with other transcripts gets a new
  `.cache` (only `tx_gene_annotation` is carried over), so it neither
  replaces the entries of the parent nor keeps them alive.

## Performance reference (72k transcripts x 25 samples, GENCODE v33)

Measured on 2026-10-15 with `Rscript` on a laptop, on the complete placenta
quantification of GENCODE v33 used before the bundled data moved to GENCODE
50; not rerun since.

| Operation | Time |
|-----------|------|
| `prepare_isoformic_annotation()` full human GFF3 (gz) | ~5 s |
| Constructor (validation against 249k transcripts) | 0.7 s |
| `gene_assays()` (derive, then cached) | 0.16 s / 0.02 s |
| `isoform_fraction()` | 0.1 s |
| `tx_diversity()` | 2.3 s |
| `dominant_isoform(group = )` | 1.6 s |
| `exon_diversity()` (1.5M exon rows, DuckDB) | 1.9 s |
| `isoform_switches()` | 0.9 s |
| `write_isoformic()` / `read_isoformic()` | 0.5 s / 0.4 s |
| `make_tx_to_gene()` GENCODE FASTA (gz) | 1.7 s |

## Reference assemblies

- GRCh38/GRCm39 (GENCODE) and T2T-CHM13v2.0 (NCBI RefSeq through
  `download_reference(reference = "t2t")`, or consortium Liftoff/CAT GFF3)
  are supported today; the object is assembly-agnostic as long as assay row
  names match the annotation's transcript IDs.
- Trio / haplotype-resolved assemblies (HPRC, Verkko/hifiasm trio binning):
  `merge_annotations()` combines per-haplotype annotation files into one
  directory (`haplotype` column, `_<haplotype>` identifier suffix) and
  `haplotype_summary()` compares haplotypes through `source_gene_id` /
  `source_transcript_id`. `download_hprc_annotation()` fetches the HPRC
  year-1 CAT files.
- Pangenome graphs (Minigraph-Cactus/vg, `vg rna` + `rpvg`): haplotype
  specific transcript quantifications can be loaded through
  `IsoformicExperiment(assays = ...)` once collapsed to the reference
  transcript IDs; graph-native transcript paths are not supported.

## Sequences, ORFs and productivity evidence

- `R/sequences.R`: `annot_data_sequences.parquet` is an optional table of the
  annotation directory. `stream_fasta()` reads any FASTA in 64 MB byte chunks
  (`readChar`, not `readLines`: a 3 GB genome takes minutes, not an hour) and
  only keeps records the caller asks for; exon sequences are cut with
  `substring()` per chromosome and assembled per transcript in DuckDB
  (`string_agg ... ORDER BY`, reverse complement with `translate`/`reverse`).
- `R/splice-sites.R`: `annot_data_splice_sites.parquet` is an optional table
  of the unique introns (`seqid`, `strand`, `intron_start`, `intron_end`)
  with the donor window (3 exonic + 6 intronic bases) and the acceptor
  window (20 intronic + 3 exonic), in the direction of transcription. It is
  written by `sequences_from_genome()` while a chromosome is in memory, so
  `extract_transcript_sequences(genome_fasta = )` gives it for free and
  `extract_splice_sites()` adds it alone. `splice_site_strength()` trains
  log-odds matrices on the introns of canonical transcripts (all introns
  below 100), background = base composition of the training windows, half a
  pseudocount per base, and scores every intron; windows with a base other
  than A, C, G, T are `NA`. The matrices are trained once per directory and
  session (`splice_site_model()`, key `splice_pwm::<dir>` of
  `annotation_cache`, stale when the splice sites, exons or transcripts
  change) on the whole annotation, so a gene gets the same scores in a
  restricted call. With `gene_name` or an object, only those transcripts,
  their exons and their sites are read; on an object the `by = "transcript"`
  table is cached and later gene calls take their rows from it. The
  per-transcript summary (`splice_strength_by_transcript()`) uses
  `tabulate()` and one radix ordering, and `mean()` only for transcripts with
  two or more scored introns. GENCODE 50: 698,612 introns, 97.6% GT-AG,
  260,664 training introns, median donor 8.3 and acceptor 7.7 bits,
  extraction 94 s, file 8.4 MB; the genome-wide table by transcript 39 s
  and 1.3 GB, one gene 9.5 s the first time and 0.3 s after.
  `merge_annotations()` does not carry the table.
- `R/orfs.R`: `transcript_orfs()` works on positions, never on translated
  strings: stop codon positions come from one `gregexpr` per sequence (stop
  codons cannot overlap), the first in-frame stop is a filtered vector
  lookup, the longest ORF per frame uses `findInterval`. Coordinates move
  between transcript and genome with the DuckDB containment joins in
  `transcript-features.R` (`genomic_to_transcript_pos()` is the row-indexed
  variant to use when the same transcript appears several times, e.g. the
  canonical of every isoform; `genomic_to_transcript_coord()` assumes unique
  IDs and silently drops repeats).
- Enrichment: the three functions (`run_enrichment()` in
  `R/enrichment.R`, `run_productivity_enrichment()` in `R/productivity.R`,
  `run_apa_enrichment()` in `R/polyadenylation.R`) share the helpers of
  `R/enrichment.R`:
  `geneset_collections()` (one collection or a named list of them),
  `enrich_collections()` (one fgsea run per collection, `collection` column
  only for several), `gene_rank_enrichment()` (gene tables ranked by one
  column) and `certainty_z()` (estimate over its standard deviation over the
  replicates; no deviation is the most certain value, not a missing one).
  A new enrichment function must use them, return `experiment` and `rank_by`
  (the columns of `empty_fgsea()` when nothing is tested) and take a table
  method next to the object method (`run_productivity_enrichment()` takes a
  `productive_fraction()` table). `rank_by` choices: `log2FC`, `signed_p`,
  `lfc_z` (`run_enrichment()`); `delta`, `delta_z`
  (`run_productivity_enrichment()`); `delta`, `utr3_length`
  (`run_apa_enrichment()`). `rank_by = "lfc_z"` has the contrast of
  the stored fold-change table (`dea$confidence`); `run_enrichment()` warns
  (`isoformic_contrast_mismatch`) when it is the reverse of the stored
  results. The fgsea calls (collections x transcript categories) are one
  task list run by `run_tasks()` (`R/parallel.R`): in the session, in the
  old order, with one worker; on `mirai` daemons of the compute profile
  `"isoformic"` with `workers` above 1 (`resolve_workers()`, `"auto"` via
  `parallelly::availableCores()`), each task with its own seed and only
  its ranking and the gene sets cut to it (`sets_in_ranking()`). Task
  functions go through `worker_function()` (base environment, no srcref)
  so daemons never load isoformic; warnings come back as values and are
  raised again in the session. The daemons are kept for the session
  (`ensure_workers()`): started on first use with one per task up to
  `workers`, grown with `mirai::launch_local()`, restarted when a call asks
  for fewer than are running, and replaced when they are no longer connected
  (`mirai::status()$connections`). Every start waits for the connections
  and runs under `withr::with_preserve_seed()`. `wait_tasks()` polls the
  tasks every 10 ms and the connections once a second, and cancels the
  remaining tasks (`mirai::stop_mirai()`) on an interrupt, an error or no
  connected daemon for `the_workers$grace` (30 s). Daemons the user started
  on the profile are used as they are and never stopped by isoformic, except
  by `isoformic_workers(0)`; `.onUnload` stops isoformic's own. The
  evaluation is in `eval-parallelism/NOTES.md`.
- The ORF rows are kept in `annot_data_orfs.parquet` of the annotation
  directory (`stored_transcript_orfs()`): an ORF depends on the transcript
  and on the canonical isoform of its gene among the rows given, so a stored
  row is reused when its `canonical_transcript_id` is the canonical of the
  call. The genes that differ are computed; of those, only the rows computed
  with the canonical isoform of the whole annotation are stored, so that an
  object without that isoform does not replace the rows the others read. The
  NMD class is derived again on every call (`orf_nmd_class()`), so
  `nmd_rule`, `start_proximal` and `long_utr3` do not invalidate the table;
  `min_orf_length`, the package version, the coding potential model in
  use, the size and time of the transcripts, exons and sequences files and
  `orf_table_version()` do (so a directory copied without its modification
  times, as `fs::dir_copy()` does, computes the table again). The names
  (`transcript_name`, `gene_id`, `gene_name`) are always taken from the
  transcripts of the call, because the genes table can change without the
  others.
  Increase `orf_table_version()` whenever `compute_transcript_orfs()` or one
  of its helpers changes. Nothing is written into
  the files of the installed package (the bundled example annotation), into
  a read-only directory or with `options(isoformic.orf_cache = FALSE)`.
- Coding potential: Fickett tables as in CPAT, hexamer log-odds and the
  logistic coefficients live in `inst/extdata/coding_model/` and are
  retrained by `data-raw/build_coding_model.R` (`just coding-model`). The
  shipped model is trained on GENCODE 50; release and held-out accuracy are
  in `data-raw/README.md`.
- `R/frozen-data.R`: the registry pins versions, URLs and MD5s; builders are
  exported functions so the processing is testable with tiny synthetic source
  files (`tests/testthat/test-frozen.R`). `read_tabular_lines()` is the
  DuckDB reader for header-less tab files (BED, GTF, GFF) and counts leading
  `#` lines itself, since the DuckDB sniffer does not honour `comment` when
  columns are given. DuckDB `least()`/`greatest()` ignore NULLs, so
  coverage sums over LEFT JOINs must guard the NULL side explicitly.
- `R/productivity-evidence.R` writes evidence columns into the transcript
  table through `add_annotation_columns()` and clears the object cache;
  `transcript_productivity()` reads them by name (`appris_principal`,
  `riboseq_main_orf`, `nmd_validated_exon`, `productivity_evidence`).
- `R/utils-stats.R`: `two_group_proportion_test()` is the closed-form
  quasi-binomial GLM shared by `test_productive_fraction()` and `test_apa()`
  (pooled proportions, Pearson dispersion floored at 1, t on
  `n_samples - 2` df, 0.5 continuity correction).
- `R/polyadenylation.R`: 3' end groups come from `transcript_features()`
  (`tes_group`); ranks follow the direction of transcription; the atlas join
  is a DuckDB window query taking the nearest site within `window` bp.
  `test_apa()` and `run_apa_enrichment()` live in the same file.

## Quantification uncertainty

Definitions in the literature, how the code compares and the backlog
(U1 to U10; U1 is a deviation of `inferential_rv()` from fishpond):
`UNCERTAINTY-LITERATURE.md`.

- `R/uncertainty-io.R`: readers for the files Salmon writes. `bootstraps.gz`
  is raw little-endian doubles, `num_bootstraps` replicates of the
  `names.tsv.gz` transcripts (same layout in Salmon 1.x and 2.x);
  `eq_classes.txt.gz` is tab-separated (class size, 0-based transcript
  indices, optional weights, count) and parsed with `strsplit()`/`unlist()`
  without a per-line loop; `meta_info.json` is read with a regex, there is no
  JSON dependency.
- `R/uncertainty.R`: draws live in `metadata$inf_reps` (a list of transcript x
  sample count matrices, the layout `productive_fraction()` already used),
  which `write_isoformic()` writes one wide file per draw; the derived tables live in `dea$uncertainty`,
  `dea$ambiguity` and `dea$confidence`, which `standardize_dea()` keeps and
  `write_isoformic()` persists.
- Inferential relative variance follows `fishpond::computeInfRV()`
  (`max(var - mean, 0) / (mean + 5) + 0.01` per sample, averaged; one
  function, `inferential_rv()`, also used by `cut_transcript_trees()`). With 30 draws
  the sampling noise of a variance puts the value of a well-determined
  transcript near 0.1, hence the default cutoff of 0.5.
- Sibling correlation is computed per sample from standardised draws for all
  pairs inside a gene (pairs are generated vectorised, 50,000 at a time) and
  averaged; ambiguity groups are the connected components of the pairs at or
  below the cutoff (label propagation).
  The same pass gives `collapse_gain` per pair (InfRV of the sum, from the
  two variances and the correlation, minus the mean InfRV of the two; the
  Terminus score). `cor_cutoff = "adaptive"` is resolved by
  `resolve_cor_cutoff()` (minus the 97.5th percentile of the largest mean
  correlation per isoform); the number used and the rule (`fixed` or
  `adaptive`) are stored in the metadata.
- `R/identifiability.R`: `equivalence_identifiability()`, structural
  identifiability from class membership (unit vector in the row space of the
  class by transcript matrix). `identifiable_from_classes()` works per
  sample: `reduce_class_system()` peels classes with one unresolved member
  and merges transcripts with identical classes (hashed columns);
  `solve_class_blocks()` splits the rest into connected blocks, solved by
  `dense_identifiable()` (`svd()`, up to `max_dense` columns) or
  `sparse_identifiable()` (CGLS projection of two random vectors on the null
  space, `Matrix`). Off by default in `add_uncertainty()` because of its
  cost: `identifiability = TRUE` (or the table) goes through
  `join_identifiability()` into the read support table, so
  `stored_equivalence()` keeps it when the tables are recomputed, and
  `compute_de_confidence()` reports `identifiable` and adds its flag to
  `flags`, never combining it with `unstable` (groups of `collapse_groups`
  are not assessed). Weights of
  `--dumpEqWeights` are ignored. Evaluation and the reasoning about ranges
  of optima: `UNCERTAINTY-LITERATURE.md` section 5.4, U5.
- `low_information` is `unique_fraction == 0 & uncertain`, present only when
  equivalence classes were read.
- `de_confidence()` repeats a Welch t-test on `log2(count scaled to a common
  library size + 5)` per draw and keeps the median and 2.5/97.5% quantiles of
  the log2 fold-change, the share of draws agreeing on its sign and the share
  in which it is called. Timing: 60,000 transcripts x 12 samples x 30 draws,
  15 s for tables plus contrast, 10 s for the collapsed version.
- `unstable` is `sign_agreement < min_agreement` (`NA` below `min_count`)
  and `unstable_de` a significant call with an unstable direction; `flags`
  lists the descriptive flags. The former `high`/`moderate`/`low` class was
  removed on 2026-11-20 (evidence: `UNCERTAINTY-LITERATURE.md` 8.3 and 8.5).
  Tables stored with the class go through `upgrade_confidence_table()`:
  `read_isoformic()` upgrades them with one message
  (`isoformic_legacy_confidence`), `stored_confidence()` silently; read the
  stored table through `stored_confidence()`, never `self@dea$confidence`.
- `add_uncertainty()` keeps draws for transcripts with a mean count of at
  least `min_count` (1); `draw_productive_fraction()` matches each draw to the
  transcripts by row name, one draw at a time, so the others count as zero.
- `productive_fraction()` puts the count draws on the scale of `assay`
  (`draw_assay_scale()`: `assay / counts` per cell; zero cells from the
  `length` assay when TPM-like, else an additive log fit whose sample terms
  come from at most 5000 evenly spaced transcripts of the object, so a
  `gene_name` subset gives the same values). `productive_fraction_sd` is the
  per-sample SD over the draws (running moments per gene x sample cell, only
  cells with a point fraction), averaged over samples.
- The replicates are stored as a list with one transcript x sample matrix
  per draw and never copied into a dense array: `draws_view()`
  (`R/uncertainty.R`) indexes the stored list, and its consumers read one
  draw (`draw_slice()`) or one sample (`draws_sample()`) at a time. The same
  accessors take an array, which `draws_array()` still builds for the tests.
  `read_inferential_replicates()` (`R/uncertainty-io.R`) holds one sample's
  replicate file in memory at a time: it allocates the matrices of the kept
  draws once (shared dimnames) and fills one column per sample.
  `replicate_moments()` sums with `.rowSums()` on the first columns in
  place (identical to `rowSums()` of the subset) and `gc_pacer()` collects
  after each sample; samples whose `names.tsv` has the order of the first
  skip the `match()`. Its peak is about 2.2 times the draws it returns, 1.4
  of which it keeps.
  `x[i, j]` keeps a view of the draws, not a copy (see Conventions).
  `cut_transcript_trees()` keeps the mean and variance of the tree leaves
  and of the node sums only.

- Real example data: `inst/extdata/knockdown/<series>/` (GSE129540 and
  GSE75489, ESRP knockdowns) with `knockdown_annotation/`, built by
  `data-raw/build_knockdown_example.R` from Salmon output produced by
  `data-raw/knockdown/quantify_knockdown.sh`. The bundled objects hold a gene
  subset, so `col_data$library_size` (total count of the complete
  quantification) overrides the library size taken from the counts assay or
  the draws in `compute_de_confidence()`. Equivalence classes cannot be
  subset, so the per-transcript table of `equivalence_support()` is bundled
  and passed to `add_uncertainty(equivalence = <table>)`.
- Case study data: `vignettes/articles/data/<study>/` holds one
  self-contained `write_isoformic()` directory per case study article
  (`mouse_bnst`, `prostate_ont`), outside the package (`vignettes/articles`
  is in `.Rbuildignore`). `data-raw/case_studies/quantify_*.sh` download the
  reads from NCBI SRA (S3 mirror, because ENA is far slower from the US) and
  run Salmon or Oarfish.
  `build_case_studies.R` runs swish and DESeq2 on the complete quantifications.
  It keeps the full DE tables, the annotation of the tested genes and the
  assays and exons of a gene panel (8.9 MB for `prostate_ont`, 9.7 MB for
  `mouse_bnst`).
  The mouse study adjusts both tests for two batches found in the
  sample PCA (`col_data$batch`).
  The work directory of the scripts is large (60 GB per study).
  Keep it on a data disk.
  Oarfish takes about 5 hours per nanopore sample on 8 threads.

## Benchmark of the sequence rules (GENCODE v50, 646,577 transcripts with sequence)

`benchmark_productivity()` with `use_biotype = FALSE`, `use_evidence = FALSE`
on the full annotation (primary assembly genome; the 24,093 scaffold
transcripts have no sequence), by `data-raw/benchmark_rules.R`. Rerun
2026-11-17 with the coding potential model of GENCODE 50 (commit 4a6e983);
the first v50 run (2026-11-11, GENCODE 49 model) differs only in the
`retained_intron` row (0.598), see `ALPHAGENOME-BENCHMARK.md` A7. The last
two columns are the earlier run on GENCODE v33 (227,979 transcripts).

| metric | n | value | n (v33) | value (v33) |
|---|---|---|---|---|
| `nonsense_mediated_decay` called unproductive | 91,818 | 0.998 | 15,925 | 0.989 |
| `retained_intron` called unproductive | 34,250 | 0.588 | 28,463 | 0.611 |
| `protein_coding` called productive or altered | 278,484 | 0.909 | 84,010 | 0.679 |
| canonical (MANE/Ensembl canonical) called productive | 20,134 | 0.973 | 20,019 | 0.968 |
| `lncRNA` called noncoding | 191,362 | 1.000 | 46,683 | 1.000 |

GENCODE grew mostly in complete `protein_coding` (88,305 productive, 164,955
altered productive) and `nonsense_mediated_decay` transcripts, which is why
the `protein_coding` row rose. Its shortfall is still almost entirely
`uncertain` (24,400): 18,353 transcripts whose 3' end is annotated incomplete
(`mRNA_end_NF`, no stop reachable) and 6,047 with an incomplete CDS; of the
824 called unproductive, 551 carry a premature stop by the 50-nt rule, 262 a
stop close to the start codon and 11 have no in-frame stop.
Retained-intron transcripts: 20,144 unproductive, 9,424 `altered_productive`,
3,726 uncertain, 867 noncoding (genes without a coding isoform) and 89
productive (with the GENCODE 49 model: 20,491, 9,254 and 3,549 for the first
three). Those whose projected ORF ends inside a last retained intron are
`altered_productive` (no downstream junction, so no NMD by the rule), which
is what the rule says, not what the biotype says; `use_biotype = TRUE` (the
default) keeps the annotation's call for them.
`protein_coding_LoF` (148 transcripts, 86 with sequence) is in
`unproductive_biotypes()` since 2026-11-12 (decision of the maintainer): the
reference genome has a variant that stops its translation, and the calls
describe the reference. With `use_biotype = FALSE` the rules decide (39
productive, 17 altered productive, 28 unproductive, 2 uncertain); in other
haplotypes the transcript can be translated.

The first v50 run took 9.3 minutes (16 cores, machine in use). The rerun
read the ORFs from the table cached in the annotation directory on
2026-11-13, and its classification took 49 s, so the two times are not
comparable.

Timings on v33 (one core): genome FASTA extraction 1.6 min (gzip, 3.1 GB), transcript
FASTA matching 50 s, `transcript_orfs()` 2.6 min, classification < 5 s.
Genome-derived sequences were identical to the GENCODE transcript FASTA for
all 227,751 shared transcripts.

## Not implemented (from the isoform-annotation review)

- Detained-intron candidates (need intron-level coverage, not available
  from the inputs the package reads). `ALPHAGENOME-BENCHMARK.md` has a pilot
  benchmark of AlphaGenome predictions for this (2026-11-12, rerun
  2026-11-17 with five seeds and a bootstrap over genes; 500 genes): no
  detectable gain in pooled AUC on validated NMD pairs (-0.001, 95%
  interval -0.016 to 0.013), a small gain on the retained-intron biotype
  (+0.007, 0.004 to 0.010), predicted intron retention AUC 0.554 for that
  biotype; not integrated, on judgement rather than by the stop rule.
  Splice-site strength is implemented with first-order matrices
  (`R/splice-sites.R`); a model with dependencies between positions
  (MaxEntScan) is not.
- Ranges of optima (only the structural identifiability of U5 is
  implemented): `UNCERTAINTY-LITERATURE.md` sections 5.2 and 5.4. Coverage
  anomalies and junction coverage compatibility: section 5.5 there.
- Proteomics evidence (PeptideAtlas / OpenProt peptide-to-isoform mapping):
  join through `add_annotation_columns()` once a table exists; no fetcher.
- Transcript-level NMD-inhibition fold changes (UPF1/SMG6/SMG7 knockdown):
  only the validated exon catalog is frozen; per-transcript responses need
  reprocessing of the raw data.
- Of the frozen files only `appris_e116v50.parquet` and
  `appris_domains_e116v50.parquet` are in the `frozen-data` GitHub release
  (2026-11-12, `ALPHAGENOME-BENCHMARK.md` A9), and nothing is deposited on
  Zenodo (`frozen_zenodo_record()` is `NA`); `fetch_frozen_data()` rebuilds
  the others from the sources (`data-raw/README.md`,
  `data-raw/zenodo_deposit.R`).

## Open items

- `draws_view()` reads the draws from memory; reading draw k from the
  `inf_reps/` files of an experiment directory on demand is open in the
  backlog of `PLAN.md`.
- The per-draw contrast is a screening statistic, not a calibrated test; a
  permutation-based version (as `swish`) would allow a calibrated score.
- `plot_genomic_context()` has no ideogram; sequence lengths are available in
  `annot_seq_regions.parquet` if one is wanted later.

## Usage tests, Terminus groups and Oarfish

- `R/usage.R`: `test_usage()` filters transcripts (`filter_usage()`), divides
  the counts by `overdispersion_for()` when draws exist and calls one of three
  backends (`usage_diffsplice()`, `usage_saturn()`, `usage_dexseq()`), then
  adds Simes gene p-values and Benjamini-Hochberg adjustments so that the
  backends return the same columns. Backends are in Suggests and checked with
  `check_usage_backend()`. Scaling defaults to on whenever draws are stored:
  the read-split benchmark (`data-raw/knockdown/benchmark_usage.R`) shows that
  unscaled `diffSplice` and `DEXSeq` call hundreds of false transcripts on
  technical replicates.
- `R/terminus.R`: `read_terminus_groups()` parses `clusters.txt`;
  `add_uncertainty(groups = )` stores the groups in `ambiguity_group`, which
  `ambiguous` and `collapse_groups` read, so the metadata key
  `uncertainty_ambiguity_source` records where they came from. Do not name
  metadata keys with the prefix of another key: `$` on a list partially
  matches.
- `R/oarfish.R`: Oarfish output is recognised by `quant_kind()` (a prefix with
  a `.quant` file); `resolve_quant_dirs()` returns a `kind` column and refuses
  mixed quantifiers. `read_inferential_replicates()`, `equivalence_support()`
  and `salmon_run_qc()` dispatch on it. Oarfish has no equivalence classes, so
  `unique_fraction` comes from `.ambig_info.tsv` and `n_classes` is `NA`.
- `test_usage()` takes several groups: `resolve_contrasts()` turns `NULL`,
  `"pairwise"` or pairs into named `c(numerator, denominator)` pairs; diffSplice
  and satuRn fit `~ 0 + group` once and test each contrast vector
  (`contrast_vector()`), DEXSeq is fitted per pair (its test is full against
  reduced model). Adjustments are per contrast (`adjust_usage()`).
- `run_terminus()` wraps `terminus group` and `terminus collapse` with
  `condathis` (env `terminus-env`, Terminus 0.1.0). `check_terminus_input()`
  detects equivalence-class weights from the field count of the first class
  line (2k + 2 with weights); Salmon 2 does not record the flag in
  `cmd_info.json`. Sample folders are named after the directory basename, so
  basenames must be unique. Parallel grouping uses `parallel::mclapply()`.
- Oarfish single-cell: `read_oarfish_sc_counts()` reads the coordinate Matrix
  Market file with DuckDB (`read_matrix_market()`, cells as rows, as Oarfish
  writes them; the transposed layout is accepted) and `read_oarfish_sc()`
  pseudobulks by a barcode to sample table with `tapply()`, without a sparse
  matrix dependency. Lengths come from the annotation (`tx_length`).
- Omnibus usage test: DEXSeq native likelihood-ratio test; Simes over the
  contrasts against the reference for diffSplice and satuRn
  (`omnibus_usage()`), because edgeR tests one coefficient or contrast.

## TreeTerminus trees

- `R/treeterminus.R`: `parse_newick_forest()` turns one Newick tree per line
  into parent pointers; internal nodes get ids in pre-order, so a child always
  has a larger id than its parent. `read_treeterminus()` maps the 0-based leaf
  positions to the transcript names of `quant.sf` (not the order of
  `eq_classes.txt.gz`, which also lists the decoys).
- `cut_transcript_trees()` reads one draw at a time (`draw_slice()`): it
  sums the leaves into their parents with `rowsum()`, then the nodes of every
  depth into their parents, from the deepest up, and `slice_moments()` keeps
  the mean and variance over the draws; InfRV of a node uses the same formula
  as `inferential_rv()`. Two top-down passes find, for every node, the lowest
  well-determined node above it and the topmost chosen node above it; the
  second makes the groups disjoint. Memory is (nodes + leaves) x samples,
  whatever the number of draws. Parent pointers that do not form trees are
  an error.
- The groups enter `add_uncertainty(groups = )` with the attributes `source`
  (`"treeterminus"`) and `cutoff`, stored as `uncertainty_ambiguity_source` and
  `uncertainty_tree_cutoff` for `uncertainty_methods()`.
