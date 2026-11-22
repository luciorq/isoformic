# Attach Extra Columns to the Annotation Tables

Joins columns of an external per-transcript or per-gene table into the
annotation directory (persisted in the Parquet tables, so every object
using the annotation sees them and caches refresh automatically). This
is the hook for evidence that the package does not compute itself:
FANTOM CAGE clusters, GTEx dominant isoforms per tissue, Ribo-seq
translation calls, NMD-inhibition responses, protein domain retention,
long-read full-length counts. A column named `productivity_evidence`
holding one of the classes of
[`productivity_palette()`](https://luciorq.github.io/isoformic/reference/productivity.md)
(`"productive"`, `"altered_productive"`, `"unproductive"`,
`"noncoding"`, `"uncertain"`) overrides the call of
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md);
other values are ignored.

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

## Details

Only the first row of each key of `table` is used. An annotation column
with the name of a new column is replaced, and rows of the annotation
without a match get `NA`.

## Examples

``` r
gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "evidence_demo")
#> ℹ Parsing example.gff3.gz (GFF) into /tmp/RtmpPCTCTT/evidence_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
evidence <- data.frame(
  transcript_id = tx$transcript_id[1:2],
  riboseq_translated = c(TRUE, FALSE),
  productivity_evidence = c("productive", "unproductive")
)
add_annotation_columns(annot, evidence)
#> ✔ Added 2 columns (riboseq_translated and productivity_evidence) for 2 of 141
#>   transcripts in /tmp/RtmpPCTCTT/evidence_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
head(tx[, c("transcript_id", "riboseq_translated")], 3)
#> # A tibble: 3 × 2
#>   transcript_id     riboseq_translated
#>   <chr>             <lgl>             
#> 1 ENST00000909997.2 TRUE              
#> 2 ENST00000909998.2 FALSE             
#> 3 ENST00000909999.2 NA                
```
