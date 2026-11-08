# Attach Extra Columns to the Annotation Tables

Joins columns of an external per-transcript or per-gene table into the
annotation directory (persisted in the Parquet tables, so every object
using the annotation sees them and caches refresh automatically). This
is the hook for evidence that the package does not compute itself:
FANTOM CAGE clusters, GTEx dominant isoforms per tissue, Ribo-seq
translation calls, NMD-inhibition responses, protein domain retention,
long-read full-length counts. A column named `productivity_evidence`
with values `"productive"` / `"unproductive"` overrides the call of
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md).

## Usage

``` r
add_annotation_columns(
  x,
  table,
  level = c("transcript", "gene"),
  by = NULL,
  prefix = NULL
)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- table:

  A data frame, or a path to a tabular file readable by
  [`read_tabular()`](https://luciorq.github.io/isoformic/reference/parquet.md).

- level:

  `"transcript"` (default) or `"gene"`.

- by:

  Name of the key column in `table`; defaults to
  `transcript_id`/`gene_id`, or the first column.

- prefix:

  Optional prefix added to the new column names.

## Value

The annotation directory path, invisibly.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "evidence_demo")
#> ℹ Parsing example.gff3 (GFF) into /tmp/RtmpV2WK1u/evidence_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
evidence <- data.frame(
  transcript_id = tx$transcript_id[1:2],
  riboseq_translated = c(TRUE, FALSE),
  productivity_evidence = c("productive", "unproductive")
)
add_annotation_columns(annot, evidence)
#> ✔ Added 2 columns (riboseq_translated and productivity_evidence) for 2 of 72
#>   transcripts in /tmp/RtmpV2WK1u/evidence_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
head(tx[, c("transcript_id", "riboseq_translated")], 3)
#> # A tibble: 3 × 2
#>   transcript_id     riboseq_translated
#>   <chr>             <lgl>             
#> 1 ENST00000282397.9 TRUE              
#> 2 ENST00000615611.4 FALSE             
#> 3 ENST00000617835.4 NA                
```
