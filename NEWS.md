# isoformic v0.1.1

## Reference and Annotation

* `download_reference()`: now supports `organism = c("human", "mouse")`, with `organism = "human"` being the default.
* `download_reference()`: argument `file_type = "gtf"` is the default.
* `prepare_annotation()`: Parse both GTF and GFF file formats into required annotation data.

## Differential Expression Analysis

* created `run_swish_pairwise()` function.

## Bug Fixes

* `prepare_profile_data()`: accepts `matrix` and `data.frame` as input for the `txi_transcript` argument.


# isoformic v0.1.0

* Initial support for `SummarizedExperiment` and `MultiAssayExperiment`. 

# isoformic v0.0.1

* Release of the initial workflow.
