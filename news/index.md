# Changelog

## Isoformic 0.2.0 (Development Version)

Release Date: Unreleased

Development Changelog:
[dev](https://github.com/luciorq/isoformic/compare/v0.1.2...HEAD)

This release rebuilds the package around the `IsoformicExperiment`
class. Table-based helpers remain available as data-frame methods, but
every analysis and plot is now first a method on the object.

### Added

- `IsoformicExperiment` is an S7 class modelled after
  `SummarizedExperiment` with three consistent levels: transcript
  assays, gene assays (derived by summing isoforms, cached, overridable
  with `gene_assays<-`) and exon structures from the annotation.
  Properties are validated on every change.
- Annotation is stored on disk as a directory of Parquet files written
  by DuckDB
  ([`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md),
  with `n_exons`, `tx_length`, metadata and sequence region tables). The
  parser is vocabulary-agnostic: besides GENCODE GFF3/GTF it reads
  Ensembl and Ensembl Rapid Release GFF3, NCBI RefSeq GFF3 and the
  T2T-CHM13 Liftoff/CAT annotations, keeping the
  `source_gene`/`source_transcript` attributes of lifted annotations as
  `source_gene_id`/`source_transcript_id` columns.
- [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)
  covers Ensembl (any release, species and division, including Ensembl
  Plants/Metazoa/Fungi/Protists), any NCBI assembly by accession,
  T2T-CHM13v2.0 (`reference = "t2t"`) and MANE releases;
  [`download_reference_url()`](https://luciorq.github.io/isoformic/reference/download_reference_url.md)
  caches any file by URL and
  [`download_hprc_annotation()`](https://luciorq.github.io/isoformic/reference/download_hprc_annotation.md)
  fetches the HPRC per-haplotype CAT annotations.
- Annotation tables carry transcript `tags` and an `is_canonical` flag
  (MANE Select, Ensembl canonical, RefSeq Select), used by
  [`dominant_isoform()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  to report `dominant_is_canonical`.
- Haplotype-resolved annotations:
  `prepare_isoformic_annotation(haplotype = )` labels one haplotype,
  [`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md)
  combines per-haplotype files (HPRC CAT, trio-binned assemblies) into
  one annotation directory with unique identifiers, and
  [`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  reports per-gene haplotype fractions and dominant-isoform switches
  between haplotypes.
- FANTOM CAT / FANTOM6 support:
  `download_reference(reference = "fantom")` fetches a FANTOM CAT level
  (GTF),
  [`download_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/download_fantom_cat_tables.md)
  the companion tables, and the parser keeps `geneClass` (as gene and
  fallback transcript type), `geneCategory`, `DHS_type`,
  `coding_status`, `cumulative_support` and `TIEScore`;
  [`add_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/add_fantom_cat_tables.md)
  attaches CAGE cluster promoters and the gene information table.
- Per-transcript structure: the parser aggregates `CDS` features
  (`cds_start`, `cds_end`, `cds_length`) and keeps
  `transcript_support_level`;
  [`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)
  computes ORF and UTR lengths, the 50-nt NMD rule
  (`stop_to_last_junction`, `nmd_predicted`), transcription start/end
  groups (alternative promoters and 3’ ends), the canonical isoform of
  each gene and splicing events relative to it
  ([`splicing_events()`](https://luciorq.github.io/isoformic/reference/transcript_features.md):
  skipped and novel exons, retained introns, alternative splice sites,
  first and last exons).
- Productivity as a graded, computed property:
  [`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  (`productive`, `altered_productive`, `unproductive`, `noncoding`,
  `uncertain`, with the reason for each call),
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  per sample or group with the between-group change,
  [`run_productivity_enrichment()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  (gene sets ranked by change in productive fraction),
  [`poison_exons()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  and `highlight_poison_exons` in
  [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md).
- [`add_annotation_columns()`](https://luciorq.github.io/isoformic/reference/add_annotation_columns.md)
  attaches external per-transcript or per-gene evidence (Ribo-seq, NMD
  inhibition, GTEx dominant isoforms, domain retention, …) to the
  annotation; a `productivity_evidence` column overrides the computed
  productivity.
- Transcript sequences:
  [`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  stores the spliced sequence of every transcript in the annotation
  directory (`annot_data_sequences.parquet`), from a genome FASTA
  (streamed in chunks, gzip or plain, chromosome names matched with or
  without `chr`) or a transcript FASTA;
  [`has_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  and
  [`transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  read it.
- Sequence-aware ORFs:
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  projects the canonical start codon of each gene onto every isoform and
  translates it (the ORFanage approach), so novel transcripts without
  CDS get an ORF; annotated CDS are re-read from the sequence
  (`annotation_cds_agrees`), with stop codon, genomic coordinates,
  uORFs, frame relative to the canonical ORF and an NMD class (`ptc`,
  `ptc_start_proximal`, `non_stop`, `long_3utr`, `none`).
  [`coding_potential()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  scores the longest ORF with the Fickett statistic, an in-frame hexamer
  score and a frozen logistic model trained on GENCODE 49
  ([`build_coding_model()`](https://luciorq.github.io/isoformic/reference/build_coding_model.md),
  `data-raw/build_coding_model.R`; held-out sensitivity 0.96,
  specificity 0.93).
- [`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  uses the sequence-derived ORFs when present (`use_sequence`), reports
  `source` (`sequence`, `annotation`, `evidence`), distinguishes
  alternative starts, in-frame changes, frameshifts and stops in novel
  sequence among `altered_productive` isoforms, calls non-stop decay and
  start-proximal stops, and applies evidence columns in a fixed priority
  (`productivity_evidence`, validated NMD exon, Ribo-seq main ORF,
  APPRIS principal).
- Frozen evidence datasets
  ([`frozen_datasets()`](https://luciorq.github.io/isoformic/reference/frozen_data.md),
  [`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md),
  [`read_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md),
  `build_frozen_*()`): APPRIS principal isoforms and method scores plus
  SPADE Pfam domain coordinates (e113v49), GENCODE Ribo-seq ORFs phase 1
  and 2, NMD-inducing exons validated by NMD inhibition (Karousis et
  al. 2021), the PolyASite 2.0 atlas and GTEx v8 median transcript TPM
  per tissue. Files are pinned by version and MD5, downloaded from the
  first mirror that has them
  ([`frozen_mirrors()`](https://luciorq.github.io/isoformic/reference/frozen_data.md):
  user option, the `frozen-data` GitHub release, the Zenodo snapshot of
  the processed files) or rebuilt from the sources;
  [`bundle_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  assembles the Zenodo deposit and `data-raw/README.md` documents the
  update and deposit procedure.
- New article “What makes an isoform hard to quantify?”: three samples
  quantified with 30 and 100 Gibbs draws, and how transcript type,
  reads, length, exons, isoforms in the gene, GC content, MANE and
  support tags, NAGNAG and nearby alternative splice sites, unique
  exonic bases, the closest sibling and unique fragments relate to the
  inferential uncertainty, to how reliably 30 draws measure it, and to
  the stability of TreeTerminus trees. The tables come from
  `data-raw/knockdown/isoform_features.R` and ship in
  `inst/extdata/knockdown/stability/`.
- [`read_treeterminus()`](https://luciorq.github.io/isoformic/reference/treeterminus.md)
  reads the transcript trees of TreeTerminus (Newick; leaves are
  positions in `quant.sf`, named from a Salmon directory) and
  [`cut_transcript_trees()`](https://luciorq.github.io/isoformic/reference/treeterminus.md)
  cuts them where the inferential relative variance of the summed
  replicates falls below the cutoff: an isoform the reads determine
  stays itself, one they do not is replaced by its lowest
  well-determined ancestor, and nested choices are merged so that the
  groups never overlap. The groups go to `add_uncertainty(groups = )`,
  so the confidence classes and `collapse_groups` work at the resolution
  the data support, isoform by isoform. Tested on real TreeTerminus
  output.
- [`run_terminus()`](https://luciorq.github.io/isoformic/reference/run_terminus.md)
  runs Terminus from R (`terminus group` on every Salmon directory, then
  the consensus `terminus collapse`) through `condathis`, checks first
  that every directory has replicates and equivalence-class weights,
  deletes the multi-gigabyte `delta.log` of each sample as soon as it is
  grouped, and returns the consensus groups. On the six GSE129540
  samples it reproduces the consensus of the command-line run.
- [`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md)
  reads the transcript groups of Terminus (`terminus collapse`,
  `clusters.txt`), and `add_uncertainty(groups = )` accepts them (or any
  table of `transcript_id` and `group_id`) in place of the groups built
  from the correlation of the replicates; the confidence classes,
  `collapse_groups` and
  [`uncertainty_methods()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  follow the new groups. The parser is tested on real output of Terminus
  0.1.0 on Salmon 2.8.0.
- [`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
  tests differential transcript usage with `edgeR` `diffSplice`,
  `satuRn` or `DEXSeq` (all in Suggests) on counts divided by the
  overdispersion of the replicates, and reports transcript and gene
  level p-values in one layout. With more than two groups it tests every
  level against the first, every pair (`contrasts = "pairwise"`) or a
  list of `c(numerator, denominator)` pairs, and stacks the results with
  a `contrast` column; `diffSplice` and `satuRn` fit all groups once,
  `DEXSeq` fits each pair. `contrasts = "omnibus"` tests for a
  difference among any of the groups: the DEXSeq likelihood-ratio test
  on all samples, or the Simes combination of the contrasts against the
  reference for `diffSplice` and `satuRn`.
- [`read_oarfish_sc()`](https://luciorq.github.io/isoformic/reference/read_oarfish_sc_counts.md)
  reads the single-cell output of Oarfish (`--single-cell`:
  `.count.mtx`, `.barcodes.txt`, `.features.txt`) and sums the cells of
  each sample, given a table of barcodes and samples, into an
  IsoformicExperiment (pseudobulk; TPM from the annotated transcript
  lengths).
  [`read_oarfish_sc_counts()`](https://luciorq.github.io/isoformic/reference/read_oarfish_sc_counts.md)
  returns the cell-level counts as a long table. The fixture is real
  Oarfish 0.10.3 output on simulated barcoded reads.
- [`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)
  reads long-read quantifications from Oarfish (counts, lengths, TPM
  computed from them).
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  (from the unique and ambiguous read counts),
  [`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  and
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  accept Oarfish output prefixes, including the Parquet bootstrap
  replicates (`.infreps.pq`). The bundled fixture is real Oarfish 0.10.3
  output on simulated nanopore reads. Mixing Oarfish and Salmon in one
  call is an error (`isoformic_mixed_quantifiers`).
- [`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  and
  [`scale_counts()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  put the quantification uncertainty into the counts as
  [`edgeR::catchSalmon()`](https://rdrr.io/pkg/edgeR/man/catchSalmon.html)
  does (Baldoni et al. 2024), from the replicates stored in the object;
  the estimates match `catchSalmon()` on a full 517,038-transcript
  quantification. A read-split experiment on GSE129540
  (data-raw/knockdown/validate_uncertainty.R) backs the confidence
  classes: of the isoforms that are large and stable over the draws in
  null contrasts, almost none are classed high, and the isoforms called
  in one half of the reads reproduce in the other half less often as the
  class or the inferential variance worsens. The small result tables are
  bundled and shown in the uncertainty article.
- Fixes from a code review: `write_isoformic(overwrite = TRUE)` writes
  to a sibling directory and swaps it in, so overwriting the directory
  an object was read from no longer deletes it (and refuses directories
  that are not experiments, `isoformic_not_experiment_dir`);
  [`de_table()`](https://luciorq.github.io/isoformic/reference/dea.md)
  and
  [`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  accept results that already have `gene_name`, `gene_id` or
  `transcript_name` columns;
  [`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  and
  [`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  contrast two of several groups with `levels`;
  [`join_DEG_DET()`](https://luciorq.github.io/isoformic/reference/join_DEG_DET.md)
  uses the same cutoff rule as the stored cutoffs (`>=` and `<=`) and
  gains `use_fdr`; a frozen dataset that fails its checksum after
  download is discarded and the next mirror is tried;
  [`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  no longer drops its title when a requested feature has no data; the
  assay hash behind the caches is computed once when the assays are set;
  `tools` is declared in Imports.
- Inferential replicates convert to and from a long table
  ([`inferential_replicates_table()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
  [`as_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)),
  and the bundled draws are Parquet tables instead of `.rds` files.
- Real ESRP knockdown examples for the uncertainty article: GSE129540
  (prostate cancer, `ESRP2` siRNA) and GSE75489 (H358 cells, `ESRP1/2`
  shRNA), re-quantified with Salmon 2.8.0 and 30 Gibbs draws, reduced to
  41 genes (under 1 MB). `data-raw/knockdown/quantify_knockdown.sh` and
  `data-raw/build_knockdown_example.R` rebuild them.
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  scales the counts with the `library_size` column of
  [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  when present (for objects holding a subset of the transcripts) and
  accepts the table of
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  in `equivalence`.
- Quantification uncertainty from Salmon:
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  and
  [`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  read the Gibbs/bootstrap replicates, the equivalence classes and the
  run metadata;
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  (or `read_salmon(inferential = TRUE)`) stores, per transcript, the
  inferential relative variance, the correlation with the isoforms of
  the same gene (ambiguity groups that cannot be told apart) and the
  share of uniquely mapping fragments.
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  repeats a contrast on every draw and reports the interval and sign
  agreement of the fold-change with a `high`/`moderate`/`low` confidence
  class, optionally testing ambiguous isoforms as one feature.
  [`uncertainty_report()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  and
  [`uncertainty_methods()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  summarise the results and the settings for a manuscript;
  [`plot_de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md),
  [`plot_isoform_correlation()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  and
  [`plot_inf_rv()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  draw them. `isoform_switches(confident_only = TRUE)`,
  [`dominant_isoform()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  and
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  use them. The tables are written by
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md);
  the draws stay in memory.
- [`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
  draws 30 Gibbs samples by default (it was 100; `num_gibbs = 0` skips
  them).
- [`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md)
  and
  [`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
  target Salmon 2 (the Rust rewrite; the conda environment installs
  2.8.0): options that Salmon 2 ignores (`--no-version-check`,
  `--type puff`, `--disableChainingHeuristic`, `--useVBOpt`) are no
  longer passed, single-end reads are passed with `--unmatedReads`, and
  the documentation explains that indices built with Salmon 1.x must be
  rebuilt.
- The documentation site is split into an introduction, a preprocessing
  article and one article per analysis (diversity, productivity,
  polyadenylation, enrichment, plots), plus articles on other
  annotations and on the table-based functions.
- `?isoformic-options` documents every option and environment variable
  the package reads or is affected by: the `isoformic.frozen_mirrors`
  option, the base R download options, the cache location variables and
  the settings of DuckDB and `condathis`.
- Evidence joins:
  [`add_appris()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md),
  [`add_riboseq_orfs()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
  (listed transcripts and exon-block compatibility, `riboseq_main_orf`,
  `riboseq_uorf`),
  [`add_nmd_exons()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md),
  [`domain_retention()`](https://luciorq.github.io/isoformic/reference/domain_retention.md)
  (domains retained, partial or lost relative to the canonical isoform,
  projected onto novel isoforms) and
  [`gtex_reference()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md)
  /
  [`gtex_tissues()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md)
  (healthy-tissue productive fraction and dominant isoform per gene).
- Productivity statistics:
  [`productivity_fractions()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  (expression share per class and unproductive mechanism: NMD, retained
  intron, no ORF, non-stop),
  [`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  (quasi-binomial two-group test on counts),
  `productive_fraction(min_expr = , inf_reps = , reference = )`
  (expression floor, uncertainty from Salmon inferential replicates,
  difference to the GTEx baseline), `dominant_isoform(min_expr = )` and
  [`benchmark_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  (recall of NMD and retained-intron biotypes by the rules alone,
  productive rate of coding and canonical transcripts).
- Alternative polyadenylation:
  [`polya_sites()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  ranks the 3’ ends of each gene, classifies `tandem_utr`,
  `alternative_last_exon` and `intronic` events relative to the
  canonical isoform, finds poly(A) signal hexamers in the sequence and
  validates ends against the PolyASite atlas or the GENCODE poly(A) GFF
  (`download_reference(file_type = "polya")`,
  [`read_polya_gff()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md));
  [`apa_usage()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  quantifies distal usage and weighted 3’ UTR length per sample or group
  and
  [`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  tests the shift.
- Long-read and assembled transcript catalogs: StringTie/Bambu/FLAIR
  style GTF attributes (`reference_id`, `ref_gene_id`, `ref_gene_name`)
  are kept, and
  [`add_sqanti_classification()`](https://luciorq.github.io/isoformic/reference/add_sqanti_classification.md)
  writes SQANTI3 structural categories into `transcript_type`;
  [`tx_type_palette()`](https://luciorq.github.io/isoformic/reference/tx_type_palette.md)
  includes those categories and the `unknown`/`novel` labels. Tables are
  read on first use and cached keyed by file signature, so objects
  sharing an annotation share one copy.
- Experiments persist as Parquet directories:
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
  and
  [`read_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
  (optionally self-contained with `copy_annotation = TRUE`).
- Readers:
  [`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
  reads `quant.sf` files straight into an `IsoformicExperiment` with
  DuckDB;
  [`read_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  [`write_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  and
  [`read_tabular()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  are exported;
  [`make_tx_to_gene()`](https://luciorq.github.io/isoformic/reference/make_tx_to_gene.md)
  parses FASTA headers with DuckDB and also accepts GFF3/GTF files and
  annotation directories.
- `SummarizedExperiment`-style accessors:
  [`assays()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`assay()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  `assay<-`,
  [`assay_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`gene_assays()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`gene_assay()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`gene_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`gene_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`row_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`col_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  `col_data<-`,
  [`exon_table()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`annot_metadata()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`clear_cache()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  plus [`dim()`](https://rdrr.io/r/base/dim.html),
  [`dimnames()`](https://rdrr.io/r/base/dimnames.html) and `[`.
- Differential expression results are standardised on assignment
  (`de_tx<-`, `de_gene<-` accept DESeq2, edgeR, limma and fishpond
  column names) and flagged with stored cutoffs
  ([`de_cutoffs()`](https://luciorq.github.io/isoformic/reference/dea.md),
  `de_cutoffs<-`,
  [`de_table()`](https://luciorq.github.io/isoformic/reference/dea.md)).
- Transcriptomic diversity summaries:
  [`isoform_fraction()`](https://luciorq.github.io/isoformic/reference/diversity.md),
  [`tx_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  (expressed isoforms, dominant isoform, effective number of isoforms),
  [`dominant_isoform()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  (per group, with dominance switches),
  [`isoform_switches()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  (transcript versus gene level calls) and
  [`exon_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  (constitutive versus alternative exons).
- [`as_summarized_experiment()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
  converts back to `SummarizedExperiment` at the transcript or gene
  level;
  [`as_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
  also accepts `tximport` lists.
- [`example_isoformic()`](https://luciorq.github.io/isoformic/reference/example_isoformic.md)
  loads a bundled experiment (preeclampsia placenta subset, GENCODE v33)
  used by all examples and tests.
- [`profile_data()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  and
  [`plot_tx_expr()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  methods for the object;
  [`plot_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  and
  [`read_gmt()`](https://luciorq.github.io/isoformic/reference/enrichment.md).
- [`summarize_to_gene()`](https://luciorq.github.io/isoformic/reference/summarize_to_gene.md)
  is exported and sums isoforms (abundance-weighted mean for `length`
  assays, as `tximport` does).
- Case study articles built from complete public datasets: chronic
  stress in the male mouse adBNST (PRJNA984605, Salmon and swish with
  batch adjustment) and a CLK inhibitor in prostate cancer cells
  sequenced with Oxford Nanopore long reads (GSE259334, Oarfish and
  swish on the bootstrap replicates)
  ([\#12](https://github.com/luciorq/isoformic/issues/12)). The
  processed bundles live in `vignettes/articles/data/` (not installed)
  and are rebuilt from the SRA reads by the scripts in
  `data-raw/case_studies/`.
- The articles cite the isoformic paper
  ([\#24](https://github.com/luciorq/isoformic/issues/24),
  [\#25](https://github.com/luciorq/isoformic/issues/25)) and
  `inst/CITATION` carries its DOI.

### Changed

- [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  is a pure `ggplot2` implementation reading the exon table of the
  annotation (one row per transcript, introns with strand arrows,
  expression labels and highlighted differentially expressed isoforms);
  it no longer needs `plotgardener`, `txdbmaker`, `TxDb`, `OrgDb` or
  `BSgenome` packages.
  [`plot_tx_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  is now an alias of its data frame method and accepts the legacy column
  names.
- [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  labels each transcript on the y axis, aligned with its track (the
  labels floated above the tracks and overlapped transcripts starting at
  the left edge), and prints expression labels without padded zeros. As
  a plain `ggplot2` object it renders the same on every graphics device,
  including Cairo and the RStudio plot panel
  ([\#21](https://github.com/luciorq/isoformic/issues/21)).
- [`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  facets several genes, orders the gene bar first and uses a common
  theme;
  [`plot_log2FC()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  is kept as a wrapper.
- [`plot_tx_expr()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  and
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  are S7 generics whose first argument is the object (or a table);
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  no longer takes the unused `lfc_cutoff` argument and expands gene sets
  with a single join.
- [`prepare_profile_data()`](https://luciorq.github.io/isoformic/reference/prepare_profile_data.md)
  returns the same long table as
  [`profile_data()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  (`feature_name`, `feature_type`, `gene_name`, `group`, `mean_expr`,
  `sd_expr`, `is_de`); gene expression is the sum of the isoforms (it
  was the mean).
- [`prepare_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_annotation.md)
  and
  [`prepare_exon_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_exon_annotation.md)
  read the Parquet annotation directory and return standard column names
  (`transcript_id`, `start_pos`, `end_pos`, …).
- Bundled example tables are now Parquet (`pe_abundance.parquet`,
  `pe_det.parquet`, `pe_deg.parquet`, `pe_col_data.parquet`).
- [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)
  discovers files by listing the release directory on the server instead
  of relying on hard-coded file names, so any current or future GENCODE
  release (and genome patch, e.g. `GRCh38.p15`) works. New `region`
  argument (`"all"`, `"chr"`, `"primary"`) selects the annotation
  flavour; versions are normalised (`49`, `"v49"`, `"38"` for mouse).
  Errors are specific: no internet, release not found (listing the
  latest releases), file not found in the release (listing the
  candidates), and failed or empty transfers (partial files are
  removed). Defaults: GENCODE release 49, `file_type = "gff"`,
  `output_path = ":cache:"`.
- Development tooling: `jarl.toml` (all jarl rules), a `jarl-check`
  GitHub workflow, lint and test-coverage workflows, and justfile
  recipes aligned with `condathis`.

### Fixed

- [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  and
  [`run_productivity_enrichment()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  run `fgsea` with
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  unless `BPPARAM` or `nproc` is given; the default socket cluster
  failed ([`serialize()`](https://rdrr.io/r/base/serialize.html) errors)
  or hung on Windows
  ([\#9](https://github.com/luciorq/isoformic/issues/9)).
- `stream_fasta()` drops the carriage returns of Windows (CRLF) FASTA
  files.
- DuckDB connections set their extension, secret and temporary
  directories under the session temporary directory, which silences the
  storage warnings of DuckDB 2.0.
- [`build_frozen_nmd_exons()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  no longer replaces the Ensembl `gene_id` with the StringTie `GeneID`
  of the source table.
- The default tests, examples and site build use only the bundled data;
  everything that downloads references or runs external tools runs only
  with `ISOFORMIC_FULL=true`.

### Removed

- `ContextData`, `create_context_data()` and `prepare_annotation_db()`.
- Dependencies `arrow`, `vroom`, `readr`, `purrr`, `stringr`, `tidyr`
  and `ggpubr`; suggested packages `plotgardener`, `txdbmaker`,
  `GenomicFeatures`, `IRanges`, `AnnotationHub`, `org.Hs.eg.db`,
  `BSgenome.Hsapiens.NCBI.GRCh38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `viridis`, `carrier` and `mirai`.

## Isoformic 0.1.2

Release Date: 2025-10-06

Development Changelog:
[0.1.2](https://github.com/luciorq/isoformic/compare/v0.1.1...v0.1.2)

### Added

- New DuckDB based `parse_annotation()`.
- New
  [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  function to visualize the genomic context of a gene with its
  transcripts.
- New
  [`tx_type_palette()`](https://luciorq.github.io/isoformic/reference/tx_type_palette.md)
  function to provide a default color palette for transcript types.

### Fixed

- In
  [`make_tx_to_gene()`](https://luciorq.github.io/isoformic/reference/make_tx_to_gene.md)
  output, wrong column name `entrez_id` replaced for `tx_length`.

## Isoformic 0.1.1

Release Date: 2024-06-18

Development Changelog:
[0.1.1](https://github.com/luciorq/isoformic/compare/v0.1.0...v0.1.1)

### Added

- [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md):
  now supports `organism = c("human", "mouse")`, with
  `organism = "human"` being the default.
- [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md):
  argument `file_type = "gtf"` is the default.
- [`prepare_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_annotation.md):
  Parse both GTF and GFF file formats into required annotation data.

### Fixed

- [`prepare_profile_data()`](https://luciorq.github.io/isoformic/reference/prepare_profile_data.md):
  accepts `matrix` and `data.frame` as input for the `txi_transcript`
  argument.

## Isoformic 0.1.0

Release Date: 2024-06-11

Development Changelog:
[0.1.0](https://github.com/luciorq/isoformic/compare/v0.0.1...v0.1.0)

### Added

- Release of the initial workflow.

## Isoformic 0.0.1

Release Date: 2024-06-08

Development Changelog:
[0.0.1](https://github.com/luciorq/isoformic/releases/tag/v0.0.1)

### Added

- Original workflow style code added.
