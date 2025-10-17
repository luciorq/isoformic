## Isoformic 0.1.3 (Development Version)

Release Date: Unreleased

Development Changelog: [dev](https://github.com/luciorq/isoformic/compare/v0.1.2...HEAD)

### Added

* New `output_path = ":cache:"` to store downloaded reference files in a
  dedicated cache folder.
  This is the new default behavior for `download_reference()`.
* New `get_isoformic_cache()` function to retrieve the path to the cache folder.

* New `IsoformicExperiment` class to encapsulate all workflow inputs with a
  disk-based backend and tidy interface.

* New generics and methods for the `IsoformicExperiment` class:

  * `col_data()` - Retrieve Sample Metadata
  * `row_data()` - Retrieve Metadata for Transcripts in the Assay Data
  * `annot_data()` - Retrieve Aggregated Annotation Data (transcript-centric)
  * `annot_data_transcripts()` - Retrieve Transcript Annotation Data
  * `annot_data_genes()` - Retrieve Gene Annotation Data
  * `annot_data_exons()` - Retrieve Exon Annotation Data
  * `tx_to_gene()` - Retrieve Transcript to Gene Mapping from the Annotation
  * `summarize_to_gene()` - Summarize Transcript-level Expression to Gene-level

  * `de_tx()` - Retrieve Differential Expression Results for Transcripts
  * `de_gene()` - Retrieve Differential Expression Results for Genes

* New interface for `plot_log2fc()` function to visualize log2 fold changes of
  transcripts within a gene.

* New `feature_column` argument in `plot_log2FC()` to specify the column
  name in the DE results table that contains the feature names (e.g., gene or
  transcript names).

* New `file_type = "genome_fasta"` option in `download_reference()` to
  download the genome FASTA file from GENCODE.

### Changed

* Argument `file_type` in `download_reference()` is now `"gff"` by default.
* Default GENCODE `version` in `download_reference()` is now `"49"` by default.
* Argument `output_path` in `download_reference()` is now `":cache:"`.
* Arguments `DEG_DET_table` and `selected_genes` in `plot_log2FC()` are now
  `de_data` and `feature`.

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
