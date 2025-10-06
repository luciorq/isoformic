## Isoformic [Unreleased]

### Added

* New `IsoformicExperiment` class to encapsulate transcriptomic data and annotations.

## Isoformic [0.1.2] - 2025-10-06

### Added

* New DuckDB based `parse_annotation()`.
* New `plot_genomic_context()` function to visualize the genomic context of a
  gene with its transcripts.
* New `tx_type_palette()` function to provide a color palette for transcript types.

### Fixed

* In `make_tx_to_gene()` output, wrong column name `entrez_id` replaced for `tx_length`.

## Isoformic [0.1.1] - 2024-06-18

### Added

* `download_reference()`: now supports `organism = c("human", "mouse")`, with `organism = "human"` being the default.
* `download_reference()`: argument `file_type = "gtf"` is the default.
* `prepare_annotation()`: Parse both GTF and GFF file formats into required annotation data.

### Fixed

* `prepare_profile_data()`: accepts `matrix` and `data.frame` as input for the `txi_transcript` argument.

## Isoformic [0.1.0] - 2024-06-11

### Added

* Release of the initial workflow.

## Isoformic [0.0.1] - 2024-06-08

### Added

* Original workflow style code added.

[unreleased]: https://github.com/luciorq/isoformic/compare/v0.1.2...HEAD
[0.1.2]: <https://github.com/luciorq/isoformic/compare/v0.1.1...v0.1.2>
[0.1.1]: <https://github.com/luciorq/isoformic/compare/v0.1.0...v0.1.1>
[0.1.0]: <https://github.com/luciorq/isoformic/compare/v0.0.1...v0.1.0>
[0.0.1]: <https://github.com/luciorq/isoformic/releases/tag/v0.0.1>
