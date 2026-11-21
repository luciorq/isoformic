## Isoformic 0.1.3 (Development Version)

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
  50 (`build_coding_model()`, `data-raw/build_coding_model.R`; 9,938 coding
  and 10,000 non-coding transcripts, held-out sensitivity 0.941 and
  specificity 0.919 at the 0.364 cutoff, AUC 0.978).
* `transcript_productivity()` uses the sequence-derived ORFs when present
  (`use_sequence`), reports `source` (`sequence`, `annotation`, `evidence`),
  distinguishes alternative starts, in-frame changes, frameshifts and stops
  in novel sequence among `altered_productive` isoforms, calls non-stop
  decay and start-proximal stops, and applies evidence columns in a fixed
  priority (`productivity_evidence`, validated NMD exon, Ribo-seq main ORF,
  APPRIS principal).
* Frozen evidence datasets (`frozen_datasets()`, `fetch_frozen_data()`,
  `read_frozen_data()`, `build_frozen_*()`): APPRIS principal isoforms and
  method scores plus SPADE Pfam domain coordinates (e116v50, Ensembl 116 /
  GENCODE 50), GENCODE Ribo-seq ORFs phase 1 and 2, NMD-inducing exons
  validated by NMD inhibition (Karousis et al. 2021), the PolyASite 2.0
  atlas and GTEx v8 median transcript TPM per tissue. Files are pinned by
  version and MD5, downloaded from the first mirror that has them
  (`frozen_mirrors()`: user option, the `frozen-data` GitHub release, the
  Zenodo snapshot of the processed files) or rebuilt from the sources;
  `bundle_frozen_data()` assembles the Zenodo deposit and
  `data-raw/README.md` documents the update and deposit procedure.
* New article "What makes an isoform hard to quantify?": three samples
  quantified with 30 and 100 Gibbs draws, and how transcript type, reads,
  length, exons, isoforms in the gene, GC content, MANE and support tags,
  NAGNAG and nearby alternative splice sites, unique exonic bases, the
  closest sibling and unique fragments relate to the inferential
  uncertainty, to how reliably 30 draws measure it, and to the stability of
  TreeTerminus trees. The tables come from
  `data-raw/knockdown/isoform_features.R` and ship in
  `inst/extdata/knockdown/stability/`.
* New article "Is quantification uncertainty worth checking? Genome-wide
  evidence": an evaluation pre-registered for a placenta cohort
  (control-versus-control null splits and held-out halves of the samples;
  `data-raw/uncertainty_evidence/`, `prereg.yaml` with its amendments),
  with two checks added outside the plan: mouse control-versus-control
  comparisons (amendment A19) and technical read halves of two knockdown
  series. Negative results are included. With 7 versus 7 samples no excess of false positives
  from ignoring the inferential uncertainty was detected (rate ratio 1.00;
  the conservative raw jackknife interval, 0.56 to 1.76, does not rule one
  out). The excess with 3 or 5 samples per group held only at p < 0.01, not
  at the FDR level, and the second pre-registered result is fragile once
  genes are treated as clusters. An unstable direction (`sign_agreement`
  below 0.95; a secondary flag of the analysis scripts, with only its 0.95
  threshold in the plan, made the headline after the results, one cohort,
  edgeR calls) predicted non-replication in held-out
  halves of the placenta cohort (odds ratio 0.36, gene-clustered interval
  0.31 to 0.41; lower raw replication in all 60 directed pairs of halves)
  but found only about 10% of the calls that did not replicate. The
  uncertainty defaults now follow this article (see Changed).
* `read_treeterminus()` reads the transcript trees of TreeTerminus (Newick;
  leaves are positions in `quant.sf`, named from a Salmon directory) and
  `cut_transcript_trees()` cuts them where the inferential relative variance of
  the summed replicates falls below the cutoff: an isoform the reads determine
  stays itself, one they do not is replaced by its lowest well-determined
  ancestor, and nested choices are merged so that the groups never overlap.
  The groups go to `add_uncertainty(groups = )`, so `ambiguous` and
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
  groups built from the correlation of the replicates; `ambiguous`,
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
  Its `workers` (forked `BiocParallel` workers for `satuRn` and `DEXSeq`,
  default 1) accepts `"auto"` and is checked like the `workers` of
  `run_enrichment()`; it does not read the `isoformic.workers` option.
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
  2024). `add_uncertainty()` keeps only the replicates of the transcripts
  with a mean count of at least `min_count`, but
  `read_inferential_replicates()` records the mean and variance of the
  replicates of every transcript (`moments` attribute), and the prior and
  the estimates come from them: on six GSE129540 samples with 68,245 of
  517,038 transcripts kept (`min_count = 5`), every estimate matches
  `catchSalmon()` (largest relative difference 4e-12), and the prior is 1.61,
  where the kept transcripts alone would give 5.23. A read-split
  experiment on GSE129540 (data-raw/knockdown/validate_uncertainty.R) checks
  the uncertainty columns against technical replication: the isoforms called
  in one half of the reads reproduce in the other half less often as the
  inferential variance rises. Between read halves every flag predicts
  non-replication, so this supports the information in the draws, not one
  flag in particular. The small result tables are bundled and shown in the
  uncertainty article.
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
  log2 ratio 2.4 against 0.7 at 20 to 100 reads). Takes half a minute to
  three minutes per human sample and needs `Matrix` (new in `Suggests`) for
  the largest block of transcripts.
* `add_uncertainty(identifiability = TRUE)` runs it and adds `identifiable`
  and `identifiable_fraction` to `tx_uncertainty()`; the table of
  `equivalence_identifiability()` can be passed in its place.
  `de_confidence()` then reports `identifiable` (and the flag `not
  identifiable from the equivalence classes`). Off by default, so existing
  results do not change. `uncertainty_methods()` describes it as a check on
  class membership, related to but not the same as the non-identifiability
  of Zheng et al. 2022, which is defined through ranges of optimal estimates.
* `transcript_productivity()` treats the GENCODE biotype
  `protein_coding_LoF` (148 transcripts in release 50) as unproductive when
  `use_biotype = TRUE`: a variant of the reference genome stops its
  translation.
* Inferential replicates convert to and from a long table
  (`inferential_replicates_table()`, `as_inferential_replicates()`), and the
  bundled draws are Parquet tables instead of `.rds` files. Their moments
  convert in the same way (`inferential_moments_table()`, the `moments`
  argument of `as_inferential_replicates()`).
* Real ESRP knockdown examples for the uncertainty article: GSE129540
  (prostate cancer, `ESRP2` siRNA) and GSE75489 (H358 cells, `ESRP1/2`
  shRNA), re-quantified with Salmon 2.8.0 and 30 Gibbs draws, reduced to 39
  genes (about 1.1 MB). `data-raw/knockdown/quantify_knockdown.sh` and
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
  contrast on every draw and reports the interval of the fold-change and the
  share of draws that agree on its direction (`sign_agreement`), with
  `unstable` flagging a direction that fewer than 95% of the draws agree on,
  optionally testing ambiguous isoforms as one feature. `uncertainty_report()` and
  `uncertainty_methods()` summarise the results and the settings for a
  manuscript; `plot_de_confidence()`, `plot_isoform_correlation()` and
  `plot_inf_rv()` draw them. `isoform_switches()` counts the significant
  isoforms with an unstable direction (`n_det_unstable`) and
  `drop_unstable = TRUE` leaves them out; `dominant_isoform()` and
  `productive_fraction()` use the tables too. The tables are
  written by `write_isoformic()`, and so are the draws (see below).
* The documentation site is split into an introduction, a preprocessing
  article and one article per analysis (diversity, productivity,
  polyadenylation, enrichment, plots), plus articles on other annotations and
  on the table-based functions.
* `?isoformic-options` documents every option and environment variable the
  package reads or is affected by: the `isoformic.duckdb_threads`,
  `isoformic.duckdb_memory_limit`, `isoformic.frozen_mirrors`,
  `isoformic.orf_cache` and `isoformic.workers` options, the base R download
  options, the cache location variables and the settings of DuckDB, `mirai`,
  `parallelly` and `condathis`.
* `options(isoformic.duckdb_threads = , isoformic.duckdb_memory_limit = )`
  bound the threads and the memory of the DuckDB connections that read and
  write Parquet and run the annotation queries. When they are unset, the
  connections follow a Slurm job (`SLURM_CPUS_PER_TASK`,
  `SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU`) and the CPU affinity mask
  (`taskset`), which DuckDB does not detect on its own, instead of using
  every core and 80% of the RAM of the node; a cgroup CPU quota or memory
  limit below the Slurm allocation still applies. Invalid values (including
  a memory limit below 1 MiB) raise an error that names the option.
* `run_enrichment()`, `run_productivity_enrichment()` and
  `run_apa_enrichment()` gain `workers` (default: the `isoformic.workers`
  option, or 1): above 1 their `fgsea` calls, one per gene set collection
  and transcript category, run in parallel on `mirai` daemons (`mirai` in
  `Suggests`). The daemons receive only the ranking of their call and the
  gene sets cut down to it, so their memory does not grow with the session,
  unlike forked workers (`nproc`, `BiocParallel::MulticoreParam()`), which
  tools that sum resident memory count once per worker. `workers = "auto"`
  uses `parallelly::availableCores()` (`parallelly` in `Suggests`; without
  it, one worker) and so follows a Slurm allocation, the CPU affinity mask
  and a cgroup quota. With one worker the results are those
  of earlier versions; above one every call gets its own seed, drawn from
  the session's stream, so the result does not depend on the number of
  workers but its sampled p-values differ slightly from a serial run.
  Warnings of `fgsea` raised on a worker are raised again in the session,
  naming the collection and the category. The daemons are started on first
  use, with the session's library paths and one thread each, and kept for
  the session, so a loop of small calls does not start them every time;
  the new `isoformic_workers()` starts a set number of them ahead or stops
  them (`isoformic_workers(0)`), and they also stop when the session ends.
  Results do not depend on whether the daemons were already running.
  Daemons that died (for example killed for their memory) are replaced by
  the next call, and a call whose daemons all die fails instead of waiting
  forever. Daemons started with `mirai::daemons(n, .compute = "isoformic")`
  are used as they are. `mirai` 2.3.0 or later is needed.
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
* `run_apa_enrichment()` runs a gene set enrichment analysis on the genes
  ranked by their change in distal 3' end usage (`rank_by = "delta"`) or in
  weighted 3' UTR length (`"utr3_length"`), from an object or from the table
  of `apa_usage()`.
* Long-read and assembled transcript catalogs: StringTie/Bambu/FLAIR style
  GTF attributes (`reference_id`, `ref_gene_id`, `ref_gene_name`) are kept,
  and `add_sqanti_classification()` writes SQANTI3 structural categories
  into `transcript_type`; `tx_type_palette()` includes those categories and
  the `unknown`/`novel` labels. Tables are read on first use and cached
  keyed by file signature, so objects sharing an annotation share one copy.
* Experiments persist as Parquet directories: `write_isoformic()` and
  `read_isoformic()` (optionally self-contained with `copy_annotation = TRUE`).
  The inferential replicates in `metadata$inf_reps` are written to
  `inf_reps/`, one wide Parquet file per draw, with their type and moments,
  and restored as they were (other properties, such as names of the draws,
  are dropped with a warning). One draw is written and read at a time, so the
  long table of `inferential_replicates_table()` (several times the memory
  of the draws) is no longer needed to keep them. Directories written
  before have no `inf_reps/` and are read without replicates. An object
  read back now has its draws, so `productive_fraction()`, `test_usage()`
  and `uncertainty_report()` use them as they do on the original object;
  to save without them, use `add_uncertainty(keep_draws = FALSE)` or set
  `x@metadata$inf_reps` to `NULL` before writing.
* Readers: `read_salmon()` reads `quant.sf` files straight into an
  `IsoformicExperiment` with DuckDB; `read_parquet()`, `write_parquet()` and
  `read_tabular()` are exported; `make_tx_to_gene()` parses FASTA headers with
  DuckDB and also accepts GFF3/GTF files and annotation directories.
* `SummarizedExperiment`-style accessors: `tx_assays()`, `tx_assay()`,
  `tx_assay<-`, `assay_names()`, `gene_assays()`, `gene_assay()`,
  `gene_data()`, `gene_names()`, `row_names()`, `col_names()`, `col_data<-`,
  `exon_table()`, `annot_metadata()`, `clear_cache()`, plus `dim()`,
  `dimnames()` and `[`.
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
  subset, GENCODE v50) used by all examples and tests.
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

### Changed

* Breaking (development version only; no release had these functions): the
  uncertainty defaults now follow the genome-wide evaluation in the article
  "Is quantification uncertainty worth checking? Genome-wide evidence" (see
  Added). `de_confidence()` no
  longer assigns a `high`/`moderate`/`low` confidence class. In a placenta
  cohort the class marked 83% of the calls (a filter on it would have
  removed 81% of the calls that replicated), separated replicating from
  non-replicating calls only modestly once count and test statistic were
  accounted for (odds ratio 0.68), and in control-versus-control
  comparisons keeping only its `high` class removed no more false calls than
  dropping as many of the least expressed transcripts. The columns
  `confidence` and `confident_de` are removed;
  `unstable` (`sign_agreement` below `min_agreement`, 0.95) and
  `unstable_de` (a significant call with an unstable direction) are added,
  and `reason` is renamed `flags`. Calls with an unstable direction
  replicated less often in held-out halves of that cohort (odds ratio 0.36;
  a secondary flag chosen as the headline after the results, one cohort),
  but the flag found only about 10% of the calls that did not replicate, so
  a stable direction is not a sign that a call is right.
  `isoform_switches(confident_only = )` is defunct, also as the third
  positional argument; `drop_unstable = TRUE` (by name) leaves out the
  unstable calls, members of a collapsed group included, and
  `n_det_unstable` replaces `n_det_confident`. Without `group`,
  `de_confidence(min_agreement = , inf_rv_cutoff = )` applies the new
  thresholds to the stored table, and its `n_det_sig` and `unstable_de`
  follow the current `de_tx()` results; `min_agreement` and
  `inf_rv_cutoff` must be a single number in range there and in
  `add_uncertainty()`, or they raise an error. A feature below `min_count`
  has a missing `mean_inf_rv` in `de_confidence()`, as in
  `tx_uncertainty()`, and is never flagged as "uncertain quantification".
  `uncertainty_report()$switches` needs only the transcript-level table:
  it is no longer `NULL` without a gene-level table (its `n_det_not_deg`
  is then `NA`) or without the stored fold-change stability (its "without
  unstable directions" row is then `NA`).
  `uncertainty_report()$de` counts the
  significant features per flag (a collapsed group once), with the rule of
  the ambiguity groups that were used, `uncertainty_methods()` describes the sign agreement
  instead of the classes, `plot_de_confidence()` colours isoforms by the
  direction over the replicates, and `confidence_palette()` is removed.
  `inf_rv_cutoff` (0.5) and `cor_cutoff` (-0.5) keep their values but no
  longer judge a call: an InfRV of at least 0.5 marked 87% of the expressed
  isoforms of that cohort. `read_isoformic()` upgrades fold-change tables
  written by earlier development versions.
* Breaking (development version only; no release had these functions): the
  transcript-level accessors are renamed to `tx_assay()`, `tx_assays()`
  and `tx_assay<-` (matching `gene_assay()`), so that they no longer mask
  `SummarizedExperiment::assay()` and related generics when both packages
  are attached. The old names `assay()`, `assays()` and `assay<-` are
  removed.
* Printing an `IsoformicExperiment` shows its directory relative to the
  working directory when it lies under it (the stored path stays absolute).
* `run_enrichment()` applies `min_size` and `max_size` to the number of
  genes of a set that have a ranked transcript of the category, as
  documented, rather than to its number of transcripts. Large sets that
  were dropped because their genes have many transcripts (about half of
  the Reactome sets for protein-coding transcripts on GENCODE 50, and
  hallmark sets such as E2F targets and G2M checkpoint) are now tested, so
  results change for large sets.
* `polya_signal_distance` of `polya_sites()` is the distance of the
  occurrence of the hexamer closest to the 3' end, including one that
  overlaps an earlier occurrence (as in `AATAAAATAAA`); such signals were
  missed before, so the distance can get smaller (about 1.6% of GENCODE 50
  transcripts with a signal).
* `polya_sites()` reads only the last `signal_window` nucleotides of each
  sequence and is cached in the object, so `apa_usage()`, `test_apa()` and
  `run_apa_enrichment()` reuse it.
* `equivalence_support()` sums the equivalence classes in DuckDB and
  `equivalence_identifiability()` parses them as numbers, one sample at a
  time: both use a fraction of the memory and time they did.
* `add_uncertainty(group = )` computes the median, interval and standard
  deviation of the fold-changes of all features at once, and
  `as_inferential_replicates()` builds the draws of a table in the order of
  `inferential_replicates_table()` without indexing the whole table.
* `salmon_quant()` draws 30 Gibbs samples by default (it was 100;
  `num_gibbs = 0` skips them).
* `salmon_index()` and `salmon_quant()` target Salmon 2 (the Rust rewrite;
  the conda environment installs 2.8.0): options that Salmon 2 ignores
  (`--no-version-check`, `--type puff`, `--disableChainingHeuristic`,
  `--useVBOpt`) are no longer passed, single-end reads are passed with
  `--unmatedReads`, and the documentation explains that indices built with
  Salmon 1.x must be rebuilt.
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
* The bundled example data are rebuilt as Parquet: `example_isoformic()`
  and the files under `system.file("extdata", package = "isoformic")`
  replace the CSV tables of 0.1.2 (see Removed).
* `download_reference()` discovers files by listing the release directory
  on the server instead of relying on hard-coded file names, so any current
  or future GENCODE release (and genome patch, e.g. `GRCh38.p15`) works.
  New `region` argument (`"all"`, `"chr"`, `"primary"`) selects the
  annotation flavour; versions are normalised (`49`, `"v49"`, `"38"` for
  mouse). Errors are specific: no internet, release not found (listing the
  latest releases), file not found in the release (listing the candidates),
  and failed or empty transfers (partial files are removed). Defaults:
  GENCODE release 50 (M39 for mouse), `file_type = "gff"`,
  `output_path = ":cache:"`.
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
  the replicates: about 14 times faster with about half the peak memory on
  295,000 transcripts, 12 samples and 100 replicates, with the same results.
  `inf_reps = FALSE` skips the stored replicates. With the default
  `rank_by = "delta"`, `run_productivity_enrichment()` does not read them.
* `run_productivity_enrichment()` also takes the table of
  `productive_fraction()`, so that one table serves several gene set
  collections.
* `run_productivity_enrichment(rank_by = "statistic")` ranks genes by the
  test of `test_productive_fraction()` on the `"counts"` assay, between the
  first and the last level, as the signed normal quantile of its p-value
  (the degrees of freedom differ between genes). It accounts for the reads
  of each gene and the variation between its samples, which `delta_z`, a
  measure of the quantification uncertainty only, does not. It also takes
  the table of `test_productive_fraction()`.
* `productive_fraction()`, `test_productive_fraction()` and
  `run_productivity_enrichment()` gain `min_gene_count`: genes whose mean
  count (the sum of their isoforms in the `"counts"` assay) is below it in
  any group are dropped before the statistics and the BH adjustment. Off by
  default (0).
* `add_uncertainty()`, `de_confidence()`, `inferential_overdispersion()`,
  `scale_counts()`, `uncertainty_report()` and `test_usage()` read the
  inferential replicates one at a time from the stored list instead of
  copying them into a dense array: `add_uncertainty(group = )` and
  `de_confidence()` held up to five arrays of the stored replicates at once
  and now hold one replicate at a time, with the same results.
* `run_enrichment()` gains `rank_by`: besides the fold-change
  (`"log2FC"`, the default) the transcripts can be ranked by
  `-log10(pvalue)` with the sign of the fold-change (`"signed_p"`) or by the
  median fold-change over the inferential replicates divided by its standard
  deviation (`"lfc_z"`, from the fold-change table of
  `add_uncertainty(x, group = )`), so that a fold-change the quantification
  cannot determine ranks low. The result has a `rank_by` column.
* `run_enrichment()`, `run_productivity_enrichment()` and
  `run_apa_enrichment()` take several collections of gene sets at once (a
  named list of gene set lists, for example
  `list(hallmark = read_gmt(h), reactome = read_gmt(r))`): each collection is
  tested on its own and the result has a `collection` column, by which
  `plot_enrichment()` splits its panels. Every gene set and every collection
  now needs a name: a list with a missing or empty name is an error
  (`isoformic_invalid_genesets`) instead of a set dropped without notice.
  An enrichment without result is a table with the `fgsea` columns and no
  rows.
* `transcript_orfs()` keeps the ORFs it computes in
  `annot_data_orfs.parquet` in the annotation directory, and later calls,
  sessions and objects with the same annotation read them back: for the
  534,000 transcripts of GENCODE 49 the ORFs take 15 to 45 seconds instead
  of about 11 minutes, and `transcript_productivity()` about one minute
  instead of 10, once the table is written. The results are the same.
  `options(isoformic.orf_cache = FALSE)` turns the table off.
* Computing ORFs and coding potential (`transcript_orfs()`,
  `coding_potential()`, and the productivity, polyadenylation and
  genomic context functions that use them) is about twice as fast, with the
  same results: the Fickett and hexamer scores are 3 to 7 times faster, and
  the coordinate conversions share one DuckDB connection.
* `cut_transcript_trees()` reads the inferential replicates one at a time
  and keeps only the mean and variance of the tree nodes, instead of a dense
  array of all the replicates.
* `x[i, j]` subsets everything stored by transcript or by sample: the
  inferential replicates in `x@metadata$inf_reps` and stored gene assays
  follow the assays, and the uncertainty tables are restricted to the kept
  transcripts. Dropping samples leaves the uncertainty tables as computed
  and says so (`add_uncertainty(x)` recomputes them from the subset
  replicates). Subscripts are checked (`isoformic_invalid_subset`), and must
  keep at least one transcript and one sample.
* `x[i, j]` no longer copies the inferential replicates: the subset shares
  the draws of `x` and records its transcripts and samples in the
  `selection` attribute of `metadata$inf_reps`. Every function reads the
  draws through that selection, so the results are the same as before, and
  `write_isoformic()`, `inferential_replicates_table()` and
  `inferential_moments_table()` write or return only the selected
  transcripts and samples. A subset of a contrast of 6 samples no longer
  takes a copy of their draws (about 1.3 GB with 295,000 transcripts and
  100 draws); in exchange, the draws of `x` stay in memory while the subset
  exists, and `saveRDS()` of a subset (or a copy sent to a parallel
  worker) carries them. Code that indexes `metadata$inf_reps` of a subset
  directly sees the draws of `x`, and so does a list operation such as
  `lapply()` or `[`, which drops the `selection` attribute: the functions
  that take the object still read only its kept transcripts and samples.
  `inferential_replicates_table()` and `inferential_moments_table()` take
  the object as a new `x` argument (alone, or with such a list as `reps`)
  to read the draws of its transcripts and samples in every case.
* `extract_transcript_sequences()` warns when it keeps existing sequences
  (without `overwrite = TRUE`) that come from another source than the FASTA
  it is given, for example sequences of a cached annotation directory cut
  from a genome while a transcript FASTA is given, and says how many
  transcripts of the annotation have no sequence. The stored ORF tables of
  annotation directories are computed again on first use, because the ORF
  search changed.
* `test_usage()` estimates the inferential overdispersion on the samples of
  the test, not on all the samples of the object.
* The imported packages have minimum versions, those the package is
  developed and tested with: `cli` 3.6.6, `DBI` 1.3.0, `dplyr` 1.2.1,
  `duckdb` 1.5.4, `fs` 2.1.0, `ggplot2` 4.0.3, `ggrepel` 0.9.8, `glue` 1.8.1,
  `rlang` 1.3.0, `S7` 1.0.0, `scales` 1.4.0, `tibble` 3.3.1 and `withr`
  3.0.3. The tests also pass with `duckdb` 2.0.0.
* `S7` 1.0.0 is required. It stores a `NULL` property differently, which
  made `IsoformicExperiment` objects with no `data_path` fail validation
  under the earlier code. The package now uses the `S7_on_load()`,
  `S7_on_unload()` and `S7_on_build()` hooks of that release.
* `jsonlite` and `readr` are in `Suggests`. The test that compares the
  inferential overdispersion with `edgeR::catchSalmon()` needs them.
* Lower memory use and faster subsetting on large objects, with the same
  results. `read_inferential_replicates()` (and `add_uncertainty(paths = )`)
  holds one sample at a time instead of every sample twice, sums the
  replicates of a sample without copying them and frees each sample before
  the next one is read: on 12 placenta samples and GENCODE 50 (25 draws of
  197,000 transcripts) its peak falls from about 3 to about 2.2 times the
  size of the replicates it returns, about 1.4 of which it keeps (the
  replicates and the moments of every transcript). `x[i, j]` matches the
  transcripts once instead of twice per replicate (about 3 times faster
  with 295,000 transcripts and 20 replicates), validates the object once,
  and returns the replicates without a copy when the subset keeps
  everything. A subset with other transcripts gets its own cache of
  derived tables, so it no longer replaces those of the object it came
  from. `gene_assays()`, `as_isoformic()` (IDs without `|`) and
  `write_isoformic()` no longer make spare copies of the assays.
* Calls that need a few transcripts or genes of an annotation directory read
  only those, with the same results. `transcript_sequences(transcript_id = )`
  (and so `transcript_orfs()` when it computes ORFs that are not stored yet)
  no longer loads every sequence: on GENCODE 50, 1,000 transcripts take 0.5
  seconds and 740 MB instead of 8.5 seconds and 2.8 GB.
  `prepare_exon_annotation()`, `exon_table()` (when the exons of the object
  are not loaded yet) and `transcript_features()` and `transcript_orfs()` on
  a directory with `gene_name` read only the transcripts and exons of those
  genes: for six GENCODE 50 genes, 0.4 to 2.3 seconds and under 450 MB
  instead of 7 to 42 seconds and 1 to 3.6 GB; `transcript_orfs()` without a
  stored ORF table takes 1.2 seconds and about 700 MB instead of 4 minutes
  and 6 GB.
* `splice_site_strength()` trains its position weight matrices once per
  annotation directory and session, and with `gene_name` reads and scores
  only the introns of those genes; on an `IsoformicExperiment` the table by
  transcript is cached in the object. The scores are the same. On GENCODE 50,
  one gene takes 9.5 s and 590 MB the first time and the next gene 0.3 s,
  instead of 115 s and 2.3 GB each, and the genome-wide table by transcript
  takes 39 s and 1.3 GB instead of 115 s and 2.2 GB.
* `poison_exons()` reads only the exons of the transcripts it looks at (when
  the exon table is not loaded yet) and lists the transcripts of the
  candidate exons only, with the same results: 177,000 transcripts of
  GENCODE 50 take 4 s and 560 MB instead of 30 s and 1.2 GB.
* `transcript_features()` and `splicing_events()` compare each transcript
  with the canonical isoform of its gene using integer keys and a single
  overlap join, with the same results. On all of GENCODE 50 the splicing
  events step peaks at 1.4 GB instead of 2.5 GB, and the whole of
  `transcript_features()` at 2.1 GB instead of 3.3 GB;
  `options(isoformic.duckdb_threads = 4)` lowers that step to 1 GB.
* `prepare_isoformic_annotation()` extracts the gene and transcript
  attributes only from the gene and transcript features, and from the exons
  and CDS only their transcript, exon ID and number, with the same tables.
  The comprehensive GENCODE 50 annotation with patches and haplotypes
  converts in 29 seconds with a peak of 2.1 GB, instead of 90 seconds and
  5.0 GB (GFF3), or in 41 seconds with a peak of 3.2 GB instead of 105
  seconds and 5.1 GB (GTF), with four DuckDB threads.
* `merge_annotations()` copies each table from the Parquet files of the
  haplotypes to the merged file in DuckDB instead of reading them into R,
  with the same tables, and copies the transcript sequences with one DuckDB
  thread. Merging two GENCODE 50 annotations with their transcript
  sequences takes 10 seconds and 1.8 GB instead of 54 seconds and 4.7 GB.
* `add_uncertainty()` keeps the cached derived tables of the object (such as
  those of `transcript_features()` and `transcript_productivity()`), which
  do not depend on the uncertainty tables, instead of discarding them, so
  they are not recomputed after it (about 50 seconds for the features of
  295,000 transcripts). It no longer empties the cache of the object it was
  given either.
* `as_summarized_experiment()` leaves the inferential replicates out of the
  metadata of the result.
* All the bundled examples (`example_isoformic()`, `example.gff3.gz`,
  `example_transcripts.fa.gz`, `mini_quants/`, `oarfish/`, `oarfish_sc/`),
  the ESRP knockdown examples, the placenta tables of the articles and the
  case-study bundles are rebuilt on GENCODE 50 (mouse M39), with Salmon
  2.8.0 and Oarfish 0.10.3 (the case studies with Salmon 2.11.0 and Oarfish
  0.11.0), from scripts now tracked in `data-raw/`
  (`placenta/`, `knockdown/`, `case_studies/`, `oarfish/`, with
  `prepare_extdata.R` and `build_knockdown_example.R`).
* The mini example (`mini_quants/`, `mini_draws.parquet`,
  `mini_annotation/`) now holds three preeclampsia (SRR11498039 to
  SRR11498041) and three control placentas (SRR11498077 to SRR11498079),
  quantified against the transcripts of 15 genes with 20 Gibbs draws,
  instead of SRR11498039 to SRR11498041 and SRR11498079 to SRR11498081, so
  the examples compare preeclampsia with control.
* The identifiers of the bundled examples changed with the release: code
  that hardcodes GENCODE v33 transcript or gene IDs of the example
  experiment, `example.gff3.gz` or `example_transcripts.fa.gz` (or v49 IDs of
  `mini_annotation/`) must be updated to their GENCODE 50 versions.

### Fixed

* `add_uncertainty(paths = )` on an object without a `counts` assay reads
  the replicates of the object's transcripts only. It read every transcript
  of the quantification, so the tables recomputed from the stored
  replicates and the replicates written by `write_isoformic()` covered fewer
  transcripts than the tables built from `paths`.
* `download_reference()` uses a file already in `output_path` without
  contacting the server when the request determines its name (GENCODE
  annotations, transcript and poly(A) files, MANE with an explicit
  `version`, FANTOM CAT), so a cached reference works offline and no longer
  fails on a transient error listing the release directory. For the other
  requests (GENCODE genomes, Ensembl annotations of an explicit release,
  NCBI and T2T assemblies), when the listing fails and exactly one cached
  file matches, that file is used with a message: a GENCODE genome must be
  of the release's assembly (GRCh37 up to release 19, GRCm38 up to M25), and
  an Ensembl `region = "all"` annotation must be the `chr_patch_hapl_scaff`
  file. Requests for the current Ensembl or MANE release still need the
  listing.
* `example_isoformic()` copies the bundled annotation directory to a new
  temporary directory on every call and points the object there, so the
  functions that store tables in the annotation
  (`extract_transcript_sequences()`, `transcript_orfs()`, `add_appris()`,
  `add_riboseq_orfs()`, `add_nmd_exons()`) no longer write into, or fail on,
  the installed package.
  `write_isoformic()` copies an annotation directory that lies under
  `tempdir()` into an experiment directory written outside `tempdir()` (with
  a message), so an experiment saved from `example_isoformic()` can be read
  in a later session; before, it pointed at a directory removed when the
  session ended.
* `salmon_index(decoy_fasta = )` builds a decoy-aware index (usually with the
  genome FASTA of the same release, plain or compressed with gzip, xz or
  bzip2): it writes the gentrome (transcripts followed by the decoys) and
  `decoys.txt` to a temporary directory and passes them to `salmon index`.
  The argument was accepted and ignored before. A decoy compressed with xz
  or bzip2, or with CRLF line endings anywhere in the file (plain or
  gzipped), is decompressed line by line instead of being appended as raw
  bytes, and decoy names that are also transcript names are an error.
* Subsetting an object to transcripts that have no stored inferential
  replicates (those below `min_count` in `add_uncertainty()`) keeps the
  replicate moments, so `scale_counts()` and `inferential_overdispersion()`
  still work on the subset. The draws themselves are not kept (zero-row
  draws), and the functions that need them treat the subset as an object
  without inferential replicates.
* `test_usage(method = "DEXSeq")` no longer prints the DESeq2 note about the
  dispersion trend fit.
* The Salmon and Terminus environments are created without the deprecation
  warning of the `packages` argument of `condathis::create_env()`: the
  package spec is passed unnamed, which every `condathis` version accepts.
* The `...` argument of the generics whose methods do not use it is
  documented as not used.
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
* `productive_fraction()` with inferential replicates stops when a level of
  `group` is named `sd`, whose column `productive_fraction_sd` used to be
  overwritten by the standard deviation over the replicates.
* `productive_fraction_sd` of `productive_fraction()` is now the standard
  deviation of each gene's productive fraction over the inferential
  replicates, around their mean, averaged over the samples. It was the root
  mean square distance of the draws to the fraction of the assay, which also
  counted the gap between the mean of Gibbs draws and the EM estimate, and,
  with the default `assay = "abundance"`, compared fractions of TPM with
  fractions of the counts the draws hold. The draws are now put on the scale
  of `assay` first (multiplied by `assay / counts` of each transcript and
  sample; where a count is zero, from the `length` assay or a fit over the
  other transcripts and samples, see `?productivity`), so
  `delta_productive_fraction_sd`, `delta_z` and `rank_by = "delta_z"` of
  `run_productivity_enrichment()` compare TPM changes with TPM spread, and
  `min_expr` applies on one scale. The draws are taken to be counts, as
  `read_inferential_replicates()` returns them: draws on another scale given
  to an object with a `counts` assay are rescaled too, while an object
  without one uses them as they are. The spread counts only the samples
  where the assay gives the gene a fraction. Expect smaller `productive_fraction_sd`
  values: on 25 placenta samples with 30 Gibbs draws (GENCODE 50, chromosome
  1), the old value was a median 2.2 times the new one (90% quantile 5.2
  times), 1.9 times with `assay = "counts"`. The median
  `delta_productive_fraction_sd` with `assay = "abundance"` rises from 0.017
  to 0.023, so fewer genes reach `abs(delta_z) >= 1.96` (697 instead of 813 of
  2,064 genes; Spearman correlation of `delta_z` before and after 0.985); with
  `assay = "counts"` it barely changes (677 instead of 679).
* `productive_fraction(gene_name = )` and `min_gene_count` keep only the
  genes asked for before summing the assay and each inferential replicate,
  with the same results: one gene of chromosome 1 with 30 draws takes 0.14 s
  and about 60 MB instead of 0.6 s and 410 MB.
* `write_isoformic(copy_annotation = TRUE)` does not copy the stored ORF
  table, which was not valid for the copy.
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
* `extract_transcript_sequences()` and `extract_splice_sites()` drop the
  carriage returns of Windows (CRLF) FASTA files.
* DuckDB connections set their extension, secret and temporary directories
  under the session temporary directory, which silences the storage warnings
  of DuckDB 2.0.
* `build_frozen_nmd_exons()` no longer replaces the Ensembl `gene_id` with
  the StringTie `GeneID` of the source table.
* The default tests, examples and site build use only the bundled data;
  everything that downloads references or runs external tools runs only with
  `ISOFORMIC_FULL=true`.
* `test_productive_fraction()` and `test_apa()` no longer test genes with
  fewer than 5 reads in their smaller class (productive or not, distal or
  not), pooled over both groups, or without reads in one group. The
  continuity correction gave such genes, including those with a stray
  fraction of a read in one class, a log odds ratio of about the log
  fold-change of the gene: a change in expression, not in usage. They have
  no statistic or p-value and do not enter the BH adjustment; their log
  odds ratio is 0, or `NA` without reads in one group.
* DuckDB queries that spill to disk (under a memory limit set by the
  `isoformic.duckdb_memory_limit` option, Slurm or a cgroup) no longer fail
  because the spill directory could not be created.
* `run_enrichment()` on an object whose stored transcript results already
  have `transcript_type`, `transcript_name` or `gene_name` columns takes
  them from the annotation instead of failing.
* `as_inferential_replicates()` and the `moments` table it reads reject a
  table where one transcript and sample (and draw) has two rows, instead of
  silently leaving the missing cell at zero.

### Removed

* `assay()`, `assays()` and `assay<-`, renamed to `tx_assay()`,
  `tx_assays()` and `tx_assay<-` (see Changed).
* `ContextData` and `create_context_data()`.
* Imported packages `forcats`, `GenomicFeatures`, `ggpubr`, `IRanges`,
  `plotgardener`, `purrr`, `stringr`, `tidyr`, `txdbmaker` and `vroom`;
  `fgsea` and `readr` move to `Suggests`. Suggested packages `rmarkdown` and
  `viridis`.
* The `isoformic_intro` vignette, replaced by the articles of the
  documentation site.
* The example files of 0.1.2 under `extdata` (`DEGs_PE_fixed2.csv`,
  `DETs_fixed.csv`, `PE_1_DETs.csv`, `PE_1_counts.csv`,
  `hsapiens_REAC_subset.gmt` and the complete MSigDB Reactome collection
  `c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt`), to keep the package small.
  `reactome_subset.gmt` replaces the gene set subset; the complete placenta
  tables and the Reactome collection that the articles use are in the
  source repository, under `vignettes/articles/data/preeclampsia/`.

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
