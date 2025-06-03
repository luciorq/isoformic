## isoformic v0.1.2

### Features

* New `tx_type_palette()`.

### Bug Fixes

* In `make_tx_to_gene()` output, wrong column name `entrez_id` replaced for `tx_length`.

## isoformic v0.1.1

### Reference and Annotation

* `download_reference()`: now supports `organism = c("human", "mouse")`, with `organism = "human"` being the default.
* `download_reference()`: argument `file_type = "gtf"` is the default.
* `prepare_annotation()`: Parse both GTF and GFF file formats into required annotation data.

### Bug Fixes

* `prepare_profile_data()`: accepts `matrix` and `data.frame` as input for the `txi_transcript` argument.

## isoformic v0.1.0

* Release of the initial workflow.
