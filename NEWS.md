## Isoformic 0.2.0 (Development Version)

Release Date: Unreleased

Development Changelog: [dev](https://github.com/luciorq/isoformic/compare/v0.1.2...HEAD)

This release rebuilds the package around the `IsoformicExperiment` class.
Table-based helpers remain available as data-frame methods, but every
analysis and plot is now first a method on the object.

### Added

* `IsoformicExperiment` is an S7 class modelled after `SummarizedExperiment`
  with three consistent levels: transcript assays, gene assays (derived by
  summing isoforms, cached, overridable with `gene_assays<-`) and exon
  structures from the annotation. Properties are validated on every change.
* Annotation is stored on disk as a directory of Parquet files written by
  DuckDB (`prepare_isoformic_annotation()`, with `n_exons`, `tx_length`,
  metadata and sequence region tables). The parser is vocabulary-agnostic:
  besides GENCODE GFF3/GTF it reads Ensembl and Ensembl Rapid Release GFF3,
  NCBI RefSeq GFF3 and the T2T-CHM13 Liftoff/CAT annotations, keeping the
  `source_gene`/`source_transcript` attributes of lifted annotations as
  `source_gene_id`/`source_transcript_id` columns.
* `download_reference()` covers Ensembl (any release, species and
  division, including Ensembl Plants/Metazoa/Fungi/Protists), any NCBI
  assembly by accession, T2T-CHM13v2.0 (`reference = "t2t"`) and MANE
  releases; `download_reference_url()` caches any file by URL and
  `download_hprc_annotation()` fetches the HPRC per-haplotype CAT
  annotations.
* Annotation tables carry transcript `tags` and an `is_canonical` flag
  (MANE Select, Ensembl canonical, RefSeq Select), used by
  `dominant_isoform()` to report `dominant_is_canonical`.
* Haplotype-resolved annotations: `prepare_isoformic_annotation(haplotype = )`
  labels one haplotype, `merge_annotations()` combines per-haplotype files
  (HPRC CAT, trio-binned assemblies) into one annotation directory with
  unique identifiers, and `haplotype_summary()` reports per-gene haplotype
  fractions and dominant-isoform switches between haplotypes.
* FANTOM CAT / FANTOM6 support: `download_reference(reference = "fantom")`
  fetches a FANTOM CAT level (GTF), `download_fantom_cat_tables()` the
  companion tables, and the parser keeps `geneClass` (as gene and fallback
  transcript type), `geneCategory`, `DHS_type`, `coding_status`,
  `cumulative_support` and `TIEScore`; `add_fantom_cat_tables()` attaches
  CAGE cluster promoters and the gene information table.
* Per-transcript structure: the parser aggregates `CDS` features
  (`cds_start`, `cds_end`, `cds_length`) and keeps `transcript_support_level`;
  `transcript_features()` computes ORF and UTR lengths, the 50-nt NMD rule
  (`stop_to_last_junction`, `nmd_predicted`), transcription start/end groups
  (alternative promoters and 3' ends), the canonical isoform of each gene
  and splicing events relative to it (`splicing_events()`: skipped and
  novel exons, retained introns, alternative splice sites, first and last
  exons).
* Productivity as a graded, computed property: `transcript_productivity()`
  (`productive`, `altered_productive`, `unproductive`, `noncoding`,
  `uncertain`, with the reason for each call), `productive_fraction()` per
  sample or group with the between-group change, `run_productivity_enrichment()`
  (gene sets ranked by change in productive fraction), `poison_exons()`
  and `highlight_poison_exons` in `plot_genomic_context()`.
* `add_annotation_columns()` attaches external per-transcript or per-gene
  evidence (Ribo-seq, NMD inhibition, GTEx dominant isoforms, domain
  retention, ...) to the annotation; a `productivity_evidence` column
  overrides the computed productivity.
* Transcript sequences: `extract_transcript_sequences()` stores the spliced
  sequence of every transcript in the annotation directory
  (`annot_data_sequences.parquet`), from a genome FASTA (streamed in chunks,
  gzip or plain, chromosome names matched with or without `chr`) or a
  transcript FASTA; `has_sequences()` and `transcript_sequences()` read it.
* Sequence-aware ORFs: `transcript_orfs()` projects the canonical start
  codon of each gene onto every isoform and translates it (the ORFanage
  approach), so novel transcripts without CDS get an ORF; annotated CDS are
  re-read from the sequence (`annotation_cds_agrees`), with stop codon,
  genomic coordinates, uORFs, frame relative to the canonical ORF and an
  NMD class (`ptc`, `ptc_start_proximal`, `non_stop`, `long_3utr`, `none`).
  `coding_potential()` scores the longest ORF with the Fickett statistic,
  an in-frame hexamer score and a frozen logistic model trained on GENCODE
  49 (`build_coding_model()`, `data-raw/build_coding_model.R`; held-out
  sensitivity 0.96, specificity 0.93).
* `transcript_productivity()` uses the sequence-derived ORFs when present
  (`use_sequence`), reports `source` (`sequence`, `annotation`, `evidence`),
  distinguishes alternative starts, in-frame changes, frameshifts and stops
  in novel sequence among `altered_productive` isoforms, calls non-stop
  decay and start-proximal stops, and applies evidence columns in a fixed
  priority (`productivity_evidence`, validated NMD exon, Ribo-seq main ORF,
  APPRIS principal).
* Frozen evidence datasets (`frozen_datasets()`, `fetch_frozen_data()`,
  `read_frozen_data()`, `build_frozen_*()`): APPRIS principal isoforms and
  method scores plus SPADE Pfam domain coordinates (e116v50, Ensembl 116 / GENCODE 50), GENCODE
  Ribo-seq ORFs phase 1 and 2, NMD-inducing exons validated by NMD
  inhibition (Karousis et al. 2021), the PolyASite 2.0 atlas and GTEx v8
  median transcript TPM per tissue. Files are pinned by version and MD5,
  downloaded from the first mirror that has them (`frozen_mirrors()`: user
  option, the `frozen-data` GitHub release, the Zenodo snapshot of the
  processed files) or rebuilt from the sources; `bundle_frozen_data()`
  assembles the Zenodo deposit and `data-raw/README.md` documents the
  update and deposit procedure.
* New article "What makes an isoform hard to quantify?": three samples
  quantified with 30 and 100 Gibbs draws, and how transcript type, reads, length,
  exons, isoforms in the gene, GC content, MANE and support tags, NAGNAG and
  nearby alternative splice sites, unique exonic bases, the closest sibling and
  unique fragments relate to the inferential uncertainty, to how reliably 30
  draws measure it, and to the stability of TreeTerminus trees. The tables come
  from `data-raw/knockdown/isoform_features.R` and ship in
  `inst/extdata/knockdown/stability/`.
* `read_treeterminus()` reads the transcript trees of TreeTerminus (Newick;
  leaves are positions in `quant.sf`, named from a Salmon directory) and
  `cut_transcript_trees()` cuts them where the inferential relative variance of
  the summed replicates falls below the cutoff: an isoform the reads determine
  stays itself, one they do not is replaced by its lowest well-determined
  ancestor, and nested choices are merged so that the groups never overlap.
  The groups go to `add_uncertainty(groups = )`, so the confidence classes and
  `collapse_groups` work at the resolution the data support, isoform by
  isoform. Tested on real TreeTerminus output.
* `run_terminus()` runs Terminus from R (`terminus group` on every Salmon
  directory, then the consensus `terminus collapse`) through `condathis`,
  checks first that every directory has replicates and equivalence-class
  weights, deletes the multi-gigabyte `delta.log` of each sample as soon as it
  is grouped, and returns the consensus groups. On the six GSE129540 samples it
  reproduces the consensus of the command-line run.
* `read_terminus_groups()` reads the transcript groups of Terminus
  (`terminus collapse`, `clusters.txt`), and `add_uncertainty(groups = )`
  accepts them (or any table of `transcript_id` and `group_id`) in place of the
  groups built from the correlation of the replicates; the confidence classes,
  `collapse_groups` and `uncertainty_methods()` follow the new groups. The
  parser is tested on real output of Terminus 0.1.0 on Salmon 2.8.0.
* `test_usage()` tests differential transcript usage with `edgeR`
  `diffSplice`, `satuRn` or `DEXSeq` (all in Suggests) on counts divided by the
  overdispersion of the replicates, and reports transcript and gene level
  p-values in one layout. With more than two groups it tests every level
  against the first, every pair (`contrasts = "pairwise"`) or a list of
  `c(numerator, denominator)` pairs, and stacks the results with a `contrast`
  column; `diffSplice` and `satuRn` fit all groups once, `DEXSeq` fits each
  pair. `contrasts = "omnibus"` tests for a difference among any of the groups:
  the DEXSeq likelihood-ratio test on all samples, or the Simes combination of
  the contrasts against the reference for `diffSplice` and `satuRn`.
* `read_oarfish_sc()` reads the single-cell output of Oarfish
  (`--single-cell`: `.count.mtx`, `.barcodes.txt`, `.features.txt`) and sums
  the cells of each sample, given a table of barcodes and samples, into an
  IsoformicExperiment (pseudobulk; TPM from the annotated transcript lengths).
  `read_oarfish_sc_counts()` returns the cell-level counts as a long table.
  The fixture is real Oarfish 0.10.3 output on simulated barcoded reads.
* `read_oarfish()` reads long-read quantifications from Oarfish (counts,
  lengths, TPM computed from them). `read_inferential_replicates()`,
  `equivalence_support()` (from the unique and ambiguous read counts),
  `salmon_run_qc()` and `add_uncertainty()` accept Oarfish output prefixes,
  including the Parquet bootstrap replicates (`.infreps.pq`). The bundled
  fixture is real Oarfish 0.10.3 output on simulated nanopore reads. Mixing
  Oarfish and Salmon in one call is an error (`isoformic_mixed_quantifiers`).
* `inferential_overdispersion()` and `scale_counts()` put the quantification
  uncertainty into the counts as `edgeR::catchSalmon()` does (Baldoni et al.
  2024), from the replicates stored in the object; the estimates match
  `catchSalmon()` on a full 517,038-transcript quantification. A read-split
  experiment on GSE129540 (data-raw/knockdown/validate_uncertainty.R) backs
  the confidence classes: of the isoforms that are large and stable over the
  draws in null contrasts, almost none are classed high, and the isoforms
  called in one half of the reads reproduce in the other half less often as the
  class or the inferential variance worsens. The small result tables are
  bundled and shown in the uncertainty article.
* Splice-site strength without an external model. `extract_splice_sites()`
  (and `extract_transcript_sequences()` when given a genome FASTA) store the
  donor and acceptor windows of every intron in the annotation directory
  (`annot_data_splice_sites.parquet`); `splice_site_strength()` scores them
  with position weight matrices trained on the introns of the canonical
  transcripts of the same annotation and reports `weakest_donor` and
  `weakest_acceptor` per transcript, or the scores per intron. On GENCODE 50
  (698,612 introns) the extraction takes about 1.5 minutes and the scoring
  25 s. The scores describe the sites; on their own they separate NMD from
  protein coding transcripts only slightly (AUC 0.55), so they are evidence
  to read next to `transcript_productivity()`, not a classifier.
* More evidence in the uncertainty tables, none of which changes a default
  result. `isoform_ambiguity()` has `collapse_gain`, the score Terminus
  groups on (inferential relative variance of the sum of two isoforms minus
  the mean of their own; negative when the sum is better determined).
  `tx_uncertainty()` has `low_information`, `TRUE` for an uncertain isoform
  without a single fragment of its own (the selection of `mmcollapse`), when
  equivalence classes were read. `add_uncertainty(cor_cutoff = "adaptive")`
  derives the correlation cutoff from the data as `mmcollapse` does; the
  default stays -0.5, and on the six GSE129540 samples the adaptive value is
  -0.25, which calls 36,910 isoforms ambiguous in place of 8,076.
* `equivalence_identifiability()` reports which transcripts the equivalence
  classes of Salmon can tell apart: a transcript is identifiable when its
  abundance is fixed once the total of every class is fixed. It is a
  property of the observed classes, not of the posterior draws, and close to
  independent of them. On the six GSE129540 samples 13,790 of 116,389
  expressed isoforms are non-identifiable in every sample, and their
  estimates differ more between two halves of the reads (median absolute
  log2 ratio 2.4 against 0.7 at 20 to 100 reads). Takes one to three minutes
  per human sample and needs `Matrix` (new in `Suggests`) for the largest
  block of transcripts.
* `add_uncertainty(identifiability = TRUE)` runs it and adds `identifiable`
  and `identifiable_fraction` to `tx_uncertainty()`; the table of
  `equivalence_identifiability()` can be passed in its place.
  `de_confidence()` then counts a non-identifiable isoform as a fourth
  condition against the call (reason `not identifiable from the equivalence
  classes`) and returns an `identifiable` column. Off by default, so
  existing results do not change.
* `transcript_productivity()` treats the GENCODE biotype
  `protein_coding_LoF` (148 transcripts in release 50) as unproductive when
  `use_biotype = TRUE`: a variant of the reference genome stops its
  translation.
* Fixes from a code review: `write_isoformic(overwrite = TRUE)` writes to a
  sibling directory and swaps it in, so overwriting the directory an object
  was read from no longer deletes it (and refuses directories that are not
  experiments, `isoformic_not_experiment_dir`); `de_table()` and
  `plot_log2fc()` accept results that already have `gene_name`, `gene_id` or
  `transcript_name` columns; `test_productive_fraction()` and `test_apa()`
  contrast two of several groups with `levels`; `join_DEG_DET()` uses the same
  cutoff rule as the stored cutoffs (`>=` and `<=`) and gains `use_fdr`; a
  frozen dataset that fails its checksum after download is discarded and the
  next mirror is tried; `plot_log2fc()` no longer drops its title when a
  requested feature has no data; the assay hash behind the caches is computed
  once when the assays are set; `tools` is declared in Imports.
* Inferential replicates convert to and from a long table
  (`inferential_replicates_table()`, `as_inferential_replicates()`), and the
  bundled draws are Parquet tables instead of `.rds` files.
* Real ESRP knockdown examples for the uncertainty article: GSE129540
  (prostate cancer, `ESRP2` siRNA) and GSE75489 (H358 cells, `ESRP1/2`
  shRNA), re-quantified with Salmon 2.8.0 and 30 Gibbs draws, reduced to 41
  genes (under 1 MB). `data-raw/knockdown/quantify_knockdown.sh` and
  `data-raw/build_knockdown_example.R` rebuild them. `add_uncertainty()` scales
  the counts with the `library_size` column of `col_data()` when present (for
  objects holding a subset of the transcripts) and accepts the table of
  `equivalence_support()` in `equivalence`.
* Quantification uncertainty from Salmon: `read_inferential_replicates()`,
  `equivalence_support()` and `salmon_run_qc()` read the Gibbs/bootstrap
  replicates, the equivalence classes and the run metadata;
  `add_uncertainty()` (or `read_salmon(inferential = TRUE)`) stores, per
  transcript, the inferential relative variance, the correlation with the
  isoforms of the same gene (ambiguity groups that cannot be told apart) and
  the share of uniquely mapping fragments. `de_confidence()` repeats a
  contrast on every draw and reports the interval and sign agreement of the
  fold-change with a `high`/`moderate`/`low` confidence class, optionally
  testing ambiguous isoforms as one feature. `uncertainty_report()` and
  `uncertainty_methods()` summarise the results and the settings for a
  manuscript; `plot_de_confidence()`, `plot_isoform_correlation()` and
  `plot_inf_rv()` draw them. `isoform_switches(confident_only = TRUE)`,
  `dominant_isoform()` and `productive_fraction()` use them. The tables are
  written by `write_isoformic()`; the draws stay in memory.
* `salmon_quant()` draws 30 Gibbs samples by default (it was 100;
  `num_gibbs = 0` skips them).
* `salmon_index()` and `salmon_quant()` target Salmon 2 (the Rust rewrite;
  the conda environment installs 2.8.0): options that Salmon 2 ignores
  (`--no-version-check`, `--type puff`, `--disableChainingHeuristic`,
  `--useVBOpt`) are no longer passed, single-end reads are passed with
  `--unmatedReads`, and the documentation explains that indices built with
  Salmon 1.x must be rebuilt.
* The documentation site is split into an introduction, a preprocessing
  article and one article per analysis (diversity, productivity,
  polyadenylation, enrichment, plots), plus articles on other annotations and
  on the table-based functions.
* `?isoformic-options` documents every option and environment variable the
  package reads or is affected by: the `isoformic.frozen_mirrors` option, the
  base R download options, the cache location variables and the settings of
  DuckDB and `condathis`.
* Evidence joins: `add_appris()`, `add_riboseq_orfs()` (listed transcripts
  and exon-block compatibility, `riboseq_main_orf`, `riboseq_uorf`),
  `add_nmd_exons()`, `domain_retention()` (domains retained, partial or
  lost relative to the canonical isoform, projected onto novel isoforms)
  and `gtex_reference()` / `gtex_tissues()` (healthy-tissue productive
  fraction and dominant isoform per gene).
* Productivity statistics: `productivity_fractions()` (expression share per
  class and unproductive mechanism: NMD, retained intron, no ORF, non-stop),
  `test_productive_fraction()` (quasi-binomial two-group test on counts),
  `productive_fraction(min_expr = , inf_reps = , reference = )` (expression
  floor, uncertainty from Salmon inferential replicates, difference to the
  GTEx baseline), `dominant_isoform(min_expr = )` and
  `benchmark_productivity()` (recall of NMD and retained-intron biotypes by
  the rules alone, productive rate of coding and canonical transcripts).
* Alternative polyadenylation: `polya_sites()` ranks the 3' ends of each
  gene, classifies `tandem_utr`, `alternative_last_exon` and `intronic`
  events relative to the canonical isoform, finds poly(A) signal hexamers
  in the sequence and validates ends against the PolyASite atlas or the
  GENCODE poly(A) GFF (`download_reference(file_type = "polya")`,
  `read_polya_gff()`); `apa_usage()` quantifies distal usage and weighted
  3' UTR length per sample or group and `test_apa()` tests the shift.
* Long-read and assembled transcript catalogs: StringTie/Bambu/FLAIR style
  GTF attributes (`reference_id`, `ref_gene_id`, `ref_gene_name`) are kept,
  and `add_sqanti_classification()` writes SQANTI3 structural categories
  into `transcript_type`; `tx_type_palette()` includes those categories and
  the `unknown`/`novel` labels. Tables are read on first use and cached keyed by file
  signature, so objects sharing an annotation share one copy.
* Experiments persist as Parquet directories: `write_isoformic()` and
  `read_isoformic()` (optionally self-contained with `copy_annotation = TRUE`).
* Readers: `read_salmon()` reads `quant.sf` files straight into an
  `IsoformicExperiment` with DuckDB; `read_parquet()`, `write_parquet()` and
  `read_tabular()` are exported; `make_tx_to_gene()` parses FASTA headers with
  DuckDB and also accepts GFF3/GTF files and annotation directories.
* `SummarizedExperiment`-style accessors: `assays()`, `assay()`, `assay<-`,
  `assay_names()`, `gene_assays()`, `gene_assay()`, `gene_data()`,
  `gene_names()`, `row_names()`, `col_names()`, `col_data<-`, `exon_table()`,
  `annot_metadata()`, `clear_cache()`, plus `dim()`, `dimnames()` and `[`.
* Differential expression results are standardised on assignment
  (`de_tx<-`, `de_gene<-` accept DESeq2, edgeR, limma and fishpond column
  names) and flagged with stored cutoffs (`de_cutoffs()`, `de_cutoffs<-`,
  `de_table()`).
* Transcriptomic diversity summaries: `isoform_fraction()`,
  `tx_diversity()` (expressed isoforms, dominant isoform, effective number of
  isoforms), `dominant_isoform()` (per group, with dominance switches),
  `isoform_switches()` (transcript versus gene level calls) and
  `exon_diversity()` (constitutive versus alternative exons).
* `as_summarized_experiment()` converts back to `SummarizedExperiment` at the
  transcript or gene level; `as_isoformic()` also accepts `tximport` lists.
* `example_isoformic()` loads a bundled experiment (preeclampsia placenta
  subset, GENCODE v33) used by all examples and tests.
* `profile_data()` and `plot_tx_expr()` methods for the object;
  `plot_enrichment()` and `read_gmt()`.
* `summarize_to_gene()` is exported and sums isoforms (abundance-weighted
  mean for `length` assays, as `tximport` does).
* Case study articles built from complete public datasets: chronic stress
  in the male mouse adBNST (PRJNA984605, Salmon and swish with batch
  adjustment) and a CLK inhibitor in prostate cancer cells sequenced with
  Oxford Nanopore long reads (GSE259334, Oarfish and swish on the bootstrap
  replicates) (#12). The processed bundles live in `vignettes/articles/data/`
  (not installed) and are rebuilt from the SRA reads by the scripts in
  `data-raw/case_studies/`.
* The articles cite the isoformic paper (#24, #25) and `inst/CITATION`
  carries its DOI.
* `pixi.toml` and `pixi.lock` define an isolated development environment
  for Linux, macOS (Apple silicon) and Windows: `pixi install`, then
  `pixi run test`. Linux and macOS take all R packages from conda-forge and
  bioconda. On Windows, `pixi run r-packages` installs the packages that
  conda lacks there (`duckdb`, `condathis` and the Bioconductor packages) as
  binary packages from CRAN, r-universe and Bioconductor 3.22. On Linux and
  macOS the environment also has LaTeX (a pinned TinyTeX release, packaged
  by the recipe in `.github/pixi/tinytex`), so `pixi run check` builds the
  PDF manual.

### Changed

* `add_uncertainty()` reports when the replicates give fewer than 60
  degrees of freedom over the study (samples x (replicates - 1); Baldoni,
  Chen and Smyth 2024), the number is stored and stated by
  `uncertainty_methods()`, and `uncertainty_report()` has a new `replicates`
  table (replicate type, draws, degrees of freedom, median overdispersion).
* The documentation calls the correlation rule of `add_uncertainty()`
  "ambiguity between isoforms" instead of "identifiability", credits
  `mmcollapse` (Turro et al. 2014) and states what the replicates do not
  measure: non-identifiable abundances (Zheng et al. 2022) and estimates
  that are wrong because the annotation is incomplete (Ma and Kingsford
  2019).
* `plot_genomic_context()` is a pure `ggplot2` implementation reading the
  exon table of the annotation (one row per transcript, introns with strand
  arrows, expression labels and highlighted differentially expressed
  isoforms); it no longer needs `plotgardener`, `txdbmaker`, `TxDb`, `OrgDb`
  or `BSgenome` packages. `plot_tx_context()` is now an alias of its data
  frame method and accepts the legacy column names.
* `plot_genomic_context()` labels each transcript on the y axis, aligned
  with its track (the labels floated above the tracks and overlapped
  transcripts starting at the left edge), and prints expression labels
  without padded zeros. As a plain `ggplot2` object it renders the same on
  every graphics device, including Cairo and the RStudio plot panel (#21).
* `plot_log2fc()` facets several genes, orders the gene bar first and uses a
  common theme; `plot_log2FC()` is kept as a wrapper.
* `plot_tx_expr()` and `run_enrichment()` are S7 generics whose first
  argument is the object (or a table); `run_enrichment()` no longer takes the
  unused `lfc_cutoff` argument and expands gene sets with a single join.
* `prepare_profile_data()` returns the same long table as `profile_data()`
  (`feature_name`, `feature_type`, `gene_name`, `group`, `mean_expr`,
  `sd_expr`, `is_de`); gene expression is the sum of the isoforms (it was the
  mean).
* `prepare_annotation()` and `prepare_exon_annotation()` read the Parquet
  annotation directory and return standard column names (`transcript_id`,
  `start_pos`, `end_pos`, ...).
* Bundled example tables are now Parquet (`pe_abundance.parquet`,
  `pe_det.parquet`, `pe_deg.parquet`, `pe_col_data.parquet`).
* `download_reference()` discovers files by listing the release directory
  on the server instead of relying on hard-coded file names, so any current
  or future GENCODE release (and genome patch, e.g. `GRCh38.p15`) works.
  New `region` argument (`"all"`, `"chr"`, `"primary"`) selects the
  annotation flavour; versions are normalised (`49`, `"v49"`, `"38"` for
  mouse). Errors are specific: no internet, release not found (listing the
  latest releases), file not found in the release (listing the candidates),
  and failed or empty transfers (partial files are removed). Defaults:
  GENCODE release 50 (M39 for mouse), `file_type = "gff"`, `output_path = ":cache:"`.
* Development tooling: `jarl.toml` (all jarl rules), a `jarl-check` GitHub
  workflow, lint and test-coverage workflows, and justfile recipes aligned
  with `condathis`.

* `group` and `levels` choose the samples in the same way in every function
  that has them (`productive_fraction()`, `productivity_fractions()`,
  `test_productive_fraction()`, `run_productivity_enrichment()`,
  `apa_usage()`, `test_apa()`, `profile_data()`, `plot_tx_expr()`,
  `dominant_isoform()`, `haplotype_summary()`, `test_usage()`,
  `de_confidence()`, `add_uncertainty()`, `prepare_profile_data()`): samples
  of groups outside `levels` are ignored, as if the object had been subset
  to them, so a contrast between two of several groups needs no sub-object.
  The summary functions aborted in this case. Every returned value,
  including `as_matrix = TRUE` and the inferential standard deviations,
  comes from the samples kept. A level in `levels` without samples, or
  `levels` without `group`, is an error (`isoformic_group_levels`); samples
  without a group value are ignored with a message; factor levels without
  samples are not groups. `dominant_isoform()` and `haplotype_summary()`
  gain `levels`.
* `run_productivity_enrichment()` gains `rank_by = c("delta", "delta_z")`
  to rank genes by the change in productive fraction divided by its
  standard deviation over the inferential replicates, and passes
  `include_altered`, `min_expr` and `inf_reps` on to
  `productive_fraction()`. The result has a `rank_by` column. `delta_z` is 0,
  not `NA`, for a gene whose productive fraction changes in no replicate; a
  change that is the same in every replicate is ranked with the largest
  `delta_z`.
* `productive_fraction()` reads the inferential replicates one at a time
  and keeps only running sums, restricted to the transcripts of coding genes
  and the samples of the `levels`, instead of holding the fractions of all
  the replicates: about 18 times faster with a fraction of the memory on
  295,000 transcripts, 12 samples and 100 replicates, with the same results.
  `inf_reps = FALSE` skips the stored replicates. With the default
  `rank_by = "delta"`, `run_productivity_enrichment()` does not read them.
* `run_productivity_enrichment()` also takes the table of
  `productive_fraction()`, so that one table serves several gene set
  collections.
* `add_uncertainty()`, `de_confidence()`, `inferential_overdispersion()`,
  `scale_counts()`, `uncertainty_report()` and `test_usage()` read the
  inferential replicates one at a time from the stored list instead of
  copying them into a dense array: on 534,000 transcripts, 12 samples and
  100 replicates the memory above the object went from 19 GB to under 4 GB
  for `add_uncertainty(group = )` and `de_confidence()`, with the same
  results.
* `run_enrichment()` gains `rank_by`: besides the fold-change
  (`"log2FC"`, the default) the transcripts can be ranked by
  `-log10(pvalue)` with the sign of the fold-change (`"signed_p"`) or by the
  median fold-change over the inferential replicates divided by its standard
  deviation (`"lfc_z"`, from the confidence table of
  `add_uncertainty(x, group = )`), so that a fold-change the quantification
  cannot determine ranks low. The result has a `rank_by` column.
* `run_apa_enrichment()` runs a gene set enrichment analysis on the genes
  ranked by their change in distal 3' end usage or in weighted 3' UTR
  length, from an object or from the table of `apa_usage()`.
* `run_enrichment()`, `run_productivity_enrichment()` and
  `run_apa_enrichment()` take several collections of gene sets at once (a
  named list of gene set lists, for example
  `list(hallmark = read_gmt(h), reactome = read_gmt(r))`): each collection is
  tested on its own and the result has a `collection` column, by which
  `plot_enrichment()` splits its panels. Every gene set now needs a name: a
  list with a missing or empty name is an error
  (`isoformic_invalid_genesets`) instead of a set dropped without notice.
  An enrichment without result is a table with the `fgsea` columns and no
  rows.
* `productive_fraction()` with inferential replicates stops when a level of
  `group` is named `sd`, whose column `productive_fraction_sd` used to be
  overwritten by the standard deviation over the replicates.
* `transcript_orfs()` keeps the ORFs it computes in
  `annot_data_orfs.parquet` in the annotation directory, and later calls,
  sessions and objects with the same annotation read them back: for the
  534,000 transcripts of GENCODE 49 the ORFs take 5 seconds instead of more
  than 5 minutes, and `transcript_productivity()` about 35 seconds instead
  of 7 minutes, once the table is written. The results are the same.
  `options(isoformic.orf_cache = FALSE)` turns the table off.
* `cut_transcript_trees()` reads the inferential replicates one at a time
  and keeps only the mean and variance of the tree nodes, instead of a dense
  array of all the replicates.
* `x[i, j]` subsets everything stored by transcript or by sample: the
  inferential replicates in `metadata(x)$inf_reps` and stored gene assays
  follow the assays, and the uncertainty tables are restricted to the kept
  transcripts. Dropping samples leaves the uncertainty tables as computed
  and says so (`add_uncertainty(x)` recomputes them from the subset
  replicates). Subscripts are checked (`isoformic_invalid_subset`), and must
  keep at least one transcript and one sample.
* `test_usage()` estimates the inferential overdispersion on the samples of
  the test, not on all the samples of the object.
* The imported packages have minimum versions, those the package is
  developed and tested with: `cli` 3.6.6, `DBI` 1.3.0, `dplyr` 1.2.1,
  `duckdb` 1.5.4, `fs` 2.1.0, `ggplot2` 4.0.3, `ggrepel` 0.9.8, `glue` 1.8.1,
  `rlang` 1.3.0, `S7` 0.2.2, `scales` 1.4.0, `tibble` 3.3.1 and `withr`
  3.0.3. The tests also pass with `duckdb` 2.0.0.
* `jsonlite` and `readr` are in `Suggests`. The test that compares the
  inferential overdispersion with `edgeR::catchSalmon()` needs them.

### Fixed

* `prepare_isoformic_annotation()` kept only the first `##sequence-region`
  of a GENCODE GFF3 file (one row, `chr1`), because GENCODE writes one pragma
  before the features of each sequence and only the header was read. The
  whole file is now scanned (`annot_seq_regions.parquet` has 706 rows for
  GENCODE 50 `chr_patch_hapl_scaff`), and a sequence without a pragma keeps
  the largest feature end. Annotation directories made before this fix hold
  the single row: rebuild them with `overwrite = TRUE` if you use the table.
* The inferential relative variance (`mean_inf_rv`, `max_inf_rv`,
  `group_inf_rv` of `cut_transcript_trees()`) now is the formula of
  `fishpond::computeInfRV()`, `max(var - mean, 0) / (mean + 5) + 0.01`, as
  the documentation said. It was `max(0.01, var - mean) / (mean + 5)`: about
  0.01 lower, and far below 0.01 for transcripts whose draws vary no more
  than counts. A transcript within 0.01 below `inf_rv_cutoff` is now
  `uncertain`. The tables of the uncertainty and stability articles were
  rebuilt.
* The source files of the frozen datasets are cached per dataset version
  (`get_isoformic_cache("sources", <name>, <version>)`). With one directory
  per dataset, a new version in the registry reused the cached file of the
  old one, which has the same name, and rebuilt the old data.
* `x[, j]` no longer fails on an object with stored gene assays, and a
  factor subscript is read as IDs, not as its integer codes.
* `dominant_isoform()` and `haplotype_summary()` failed with a single group.
* `de_confidence(collapse_groups = TRUE)` labelled a transcript with the ID
  of its ambiguity group when it was the only member left in the object.
* `isoform_fraction()` returned cached fractions after `gene_assays<-`.
* Stored gene assays whose columns are in another order than the transcript
  assays gave wrong fractions in `isoform_fraction()`, `dominant_isoform()`
  and `profile_data()`; the matrices are now matched by sample ID.
* `run_enrichment()` and `run_productivity_enrichment()` run `fgsea` with
  `BiocParallel::SerialParam()` unless `BPPARAM` or `nproc` is given; the
  default socket cluster failed (`serialize()` errors) or hung on Windows
  (#9).
* `stream_fasta()` drops the carriage returns of Windows (CRLF) FASTA files.
* DuckDB connections set their extension, secret and temporary directories
  under the session temporary directory, which silences the storage warnings
  of DuckDB 2.0.
* `build_frozen_nmd_exons()` no longer replaces the Ensembl `gene_id` with
  the StringTie `GeneID` of the source table.
* The default tests, examples and site build use only the bundled data;
  everything that downloads references or runs external tools runs only with
  `ISOFORMIC_FULL=true`.

### Removed

* `ContextData`, `create_context_data()` and `prepare_annotation_db()`.
* Dependencies `arrow`, `vroom`, `readr`, `purrr`, `stringr`, `tidyr` and
  `ggpubr`; suggested packages `plotgardener`, `txdbmaker`,
  `GenomicFeatures`, `IRanges`, `AnnotationHub`, `org.Hs.eg.db`,
  `BSgenome.Hsapiens.NCBI.GRCh38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `viridis`, `carrier` and `mirai`.

## Isoformic 0.1.2

Release Date: 2025-10-06

Development Changelog: [0.1.2](https://github.com/luciorq/isoformic/compare/v0.1.1...v0.1.2)

### Added

* New DuckDB based `parse_annotation()`.
* New `plot_genomic_context()` function to visualize the genomic context of a
  gene with its transcripts.
* New `tx_type_palette()` function to provide a default color palette for
  transcript types.

### Fixed

* In `make_tx_to_gene()` output, wrong column name `entrez_id` replaced for `tx_length`.

## Isoformic 0.1.1

Release Date: 2024-06-18

Development Changelog: [0.1.1](https://github.com/luciorq/isoformic/compare/v0.1.0...v0.1.1)

### Added

* `download_reference()`: now supports `organism = c("human", "mouse")`,
  with `organism = "human"` being the default.
* `download_reference()`: argument `file_type = "gtf"` is the default.
* `prepare_annotation()`: Parse both GTF and GFF file formats into
  required annotation data.

### Fixed

* `prepare_profile_data()`: accepts `matrix` and `data.frame` as input for the `txi_transcript` argument.

## Isoformic 0.1.0

Release Date: 2024-06-11

Development Changelog: [0.1.0](https://github.com/luciorq/isoformic/compare/v0.0.1...v0.1.0)

### Added

* Release of the initial workflow.

## Isoformic 0.0.1

Release Date: 2024-06-08

Development Changelog: [0.0.1](https://github.com/luciorq/isoformic/releases/tag/v0.0.1)

### Added

* Original workflow style code added.
