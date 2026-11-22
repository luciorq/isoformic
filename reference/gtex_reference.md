# Tissue Reference Values from GTEx

Computes, from the GTEx v8 median transcript TPM of one tissue (frozen
dataset `gtex_transcript_tpm`,
[frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md))
and the productivity calls of `x`, the reference state of every gene in
healthy tissue: its productive fraction, its dominant isoform and the
dominant isoform's productivity.
[`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
accepts the result (or a tissue name) as `reference` and reports the
difference between the samples and the tissue baseline, which separates
a disease-specific loss of productive isoforms from a tissue's normal
state.

Transcript identifiers are matched without version; GTEx v8 uses GENCODE
v26, so transcripts annotated later have no reference value.

## Usage

``` r
gtex_reference(self, tissue, gtex = NULL, include_altered = TRUE, ...)

gtex_tissues(gtex = NULL)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- tissue:

  GTEx tissue name (`SMTSD`, for example `"Whole Blood"`,
  `"Brain - Cortex"`); see `gtex_tissues()`.

- gtex:

  Data frame or Parquet path of the GTEx table; `NULL` fetches the
  frozen dataset.

- include_altered:

  Logical, count `altered_productive` isoforms as productive (default
  `TRUE`).

- ...:

  Not used; reserved for methods.

## Value

`gtex_reference()` returns a tibble with one row per gene of `x` that
has a transcript in the GTEx table: `gene_id`, `gene_name`,
`reference_tissue`, `reference_gene_tpm`,
`reference_productive_fraction` (`NA` for genes without a coding isoform
or with a total TPM of 0), `reference_dominant_transcript_id`,
`reference_dominant_productivity`, `n_reference_transcripts`.
`gtex_tissues()` returns a tibble with the `tissue` names and their
`n_samples`.

## Examples

``` r
if (FALSE) { # \dontrun{
iso <- example_isoformic()
ref <- gtex_reference(iso, tissue = "Whole Blood")
productive_fraction(iso, group = "condition", reference = ref)
} # }
```
