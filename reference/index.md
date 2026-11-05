# Package index

## The IsoformicExperiment class

Transcript, gene and exon level data for the same samples, backed by
Parquet files.

- [`IsoformicExperiment()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  : IsoformicExperiment: Transcript, Gene and Exon Level Data in One
  Object
- [`col_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`` `col_data<-`() ``](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`col_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`assays()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`assay_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`assay()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`` `assay<-`() ``](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`row_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`row_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`tx_annot()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`tx_to_gene()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`gene_assays()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`` `gene_assays<-`() ``](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`gene_assay()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`gene_names()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`gene_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`annot_data_transcripts()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`annot_data_genes()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`annot_data_exons()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`annot_data()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`annot_metadata()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`exon_table()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  [`clear_cache()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  : Accessors for IsoformicExperiment Objects
- [`example_isoformic()`](https://luciorq.github.io/isoformic/reference/example_isoformic.md)
  : Example IsoformicExperiment

## Creating experiments

- [`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
  : Read Salmon Quantifications Into an IsoformicExperiment
- [`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)
  : Read Oarfish Quantifications Into an IsoformicExperiment
- [`read_oarfish_sc_counts()`](https://luciorq.github.io/isoformic/reference/read_oarfish_sc_counts.md)
  [`read_oarfish_sc()`](https://luciorq.github.io/isoformic/reference/read_oarfish_sc_counts.md)
  : Read Oarfish single-cell quantifications
- [`as_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
  [`as_summarized_experiment()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
  : Convert Between IsoformicExperiment and SummarizedExperiment
- [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
  [`read_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
  : Save and Restore IsoformicExperiment Objects as Parquet

## Annotation

- [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)
  : Convert a GFF3/GTF Annotation to Parquet Files
- [`prepare_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_annotation.md)
  : Read a Prepared Annotation into Memory
- [`prepare_exon_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_exon_annotation.md)
  : Exon Table for Selected Genes From an Annotation
- [`make_tx_to_gene()`](https://luciorq.github.io/isoformic/reference/make_tx_to_gene.md)
  : Create a Transcript-to-Gene Mapping Table
- [`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md)
  : Merge Per-Haplotype Annotations Into One Annotation Directory
- [`add_annotation_columns()`](https://luciorq.github.io/isoformic/reference/add_annotation_columns.md)
  : Attach Extra Columns to the Annotation Tables
- [`add_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/add_fantom_cat_tables.md)
  : Attach FANTOM CAT Tables to an Annotation
- [`add_sqanti_classification()`](https://luciorq.github.io/isoformic/reference/add_sqanti_classification.md)
  : Add SQANTI3 Structural Categories to an Annotation

## Reference downloads

- [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)
  : Download Reference Annotations and Sequences
- [`download_reference_url()`](https://luciorq.github.io/isoformic/reference/download_reference_url.md)
  : Download Any Reference File by URL
- [`download_hprc_annotation()`](https://luciorq.github.io/isoformic/reference/download_hprc_annotation.md)
  : Download HPRC Per-Haplotype CAT Annotations
- [`download_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/download_fantom_cat_tables.md)
  : Download the FANTOM CAT Companion Tables
- [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md)
  : Retrieve System-Dependent Cache Path for Isoformic
- [`isoformic-options`](https://luciorq.github.io/isoformic/reference/isoformic-options.md)
  [`isoformic_options`](https://luciorq.github.io/isoformic/reference/isoformic-options.md)
  [`isoformic.frozen_mirrors`](https://luciorq.github.io/isoformic/reference/isoformic-options.md)
  : Options and environment variables used by isoformic

## Quantification uncertainty

Inferential replicates, equivalence classes and run metadata from Salmon
turned into uncertainty, ambiguity and confidence of isoform-level
results.

- [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  [`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  [`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  [`inferential_replicates_table()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  [`as_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  : Read the inferential uncertainty that Salmon writes
- [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  [`tx_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  [`isoform_ambiguity()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  [`uncertainty_report()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  [`uncertainty_methods()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  : Quantify and report the uncertainty of isoform-level results
- [`confidence_palette()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  [`plot_de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  [`plot_isoform_correlation()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  [`plot_inf_rv()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)
  : Plot the uncertainty of isoform-level results
- [`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  [`scale_counts()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
  : Put the quantification uncertainty into the counts
- [`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md)
  : Read the transcript groups that Terminus builds
- [`run_terminus()`](https://luciorq.github.io/isoformic/reference/run_terminus.md)
  : Run Terminus on Salmon output from R
- [`read_treeterminus()`](https://luciorq.github.io/isoformic/reference/treeterminus.md)
  [`cut_transcript_trees()`](https://luciorq.github.io/isoformic/reference/treeterminus.md)
  : Read and cut the transcript trees of TreeTerminus
- [`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
  : Test differential transcript usage

## Isoform structure and productivity

- [`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)
  [`splicing_events()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)
  : Structural and Coding Features of Transcripts
- [`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  [`has_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  [`transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  : Extract Transcript Sequences into an Annotation Directory
- [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  [`coding_potential()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  : Open Reading Frames of Transcripts from Their Sequences
- [`productivity_palette()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`productivity_fractions()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`run_productivity_enrichment()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`poison_exons()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`benchmark_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  : Productivity of Transcripts and Genes
- [`add_appris()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
  [`add_riboseq_orfs()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
  [`add_nmd_exons()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
  : Attach External Productivity Evidence to an Annotation
- [`domain_retention()`](https://luciorq.github.io/isoformic/reference/domain_retention.md)
  : Protein Domain Retention of Isoforms
- [`gtex_reference()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md)
  [`gtex_tissues()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md)
  : Tissue Reference Values from GTEx
- [`polya_sites()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  [`read_polya_gff()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  [`apa_usage()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  [`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
  : Poly(A) Sites and Alternative Polyadenylation of Transcripts

## Frozen evidence datasets

Reference tables pinned to a source version (APPRIS, Ribo-seq ORFs,
validated NMD exons, PolyASite, GTEx) and the coding potential model.

- [`frozen_datasets()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  [`frozen_mirrors()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  [`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  [`read_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  [`bundle_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  : Frozen Reference Datasets
- [`build_frozen_appris()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  [`build_frozen_appris_domains()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  [`build_frozen_riboseq_orfs()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  [`build_frozen_nmd_exons()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  [`build_frozen_polyasite()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  [`build_frozen_gtex()`](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  : Build the Frozen Datasets from Their Sources
- [`build_coding_model()`](https://luciorq.github.io/isoformic/reference/build_coding_model.md)
  : Train the Coding Potential Model

## Differential expression

- [`de_tx()`](https://luciorq.github.io/isoformic/reference/dea.md)
  [`de_gene()`](https://luciorq.github.io/isoformic/reference/dea.md)
  [`` `de_tx<-`() ``](https://luciorq.github.io/isoformic/reference/dea.md)
  [`` `de_gene<-`() ``](https://luciorq.github.io/isoformic/reference/dea.md)
  [`de_cutoffs()`](https://luciorq.github.io/isoformic/reference/dea.md)
  [`` `de_cutoffs<-`() ``](https://luciorq.github.io/isoformic/reference/dea.md)
  [`de_table()`](https://luciorq.github.io/isoformic/reference/dea.md) :
  Differential Expression Results in an IsoformicExperiment
- [`join_DEG_DET()`](https://luciorq.github.io/isoformic/reference/join_DEG_DET.md)
  : Merge Gene and Transcript Level Differential Expression Tables
- [`is_deg_sig()`](https://luciorq.github.io/isoformic/reference/is_deg_sig.md)
  : Flag Transcripts Whose Gene Is Differentially Expressed

## Transcriptomic diversity

- [`isoform_fraction()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  [`tx_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  [`dominant_isoform()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  [`isoform_switches()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  [`exon_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  [`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)
  : Transcriptomic Diversity Summaries

## Plots

- [`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  [`plot_log2FC()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md)
  : Plot Log2 Fold-Changes of a Gene and Its Isoforms
- [`profile_data()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  [`plot_tx_expr()`](https://luciorq.github.io/isoformic/reference/profiles.md)
  : Expression Profiles of Genes and Their Isoforms Across Groups
- [`prepare_profile_data()`](https://luciorq.github.io/isoformic/reference/prepare_profile_data.md)
  : Prepare Expression Profile Data From Loose Tables
- [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  [`plot_tx_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  : Plot the Genomic Context of a Gene's Isoforms
- [`tx_type_palette()`](https://luciorq.github.io/isoformic/reference/tx_type_palette.md)
  : Transcript Type Color Palette
- [`productivity_palette()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`productivity_fractions()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`run_productivity_enrichment()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`poison_exons()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  [`benchmark_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  : Productivity of Transcripts and Genes

## Enrichment

- [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  [`default_enrichment_types()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  [`read_gmt()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  [`plot_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
  : Transcript-Type Aware Gene Set Enrichment

## Parquet and quantification tools

- [`read_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  [`read_tabular()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  [`write_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md)
  : Read and Write Parquet Files
- [`summarize_to_gene()`](https://luciorq.github.io/isoformic/reference/summarize_to_gene.md)
  : Summarize Transcript-Level Matrices to Gene Level
- [`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md)
  : Build Salmon Index
- [`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
  : Run Salmon Quantification
