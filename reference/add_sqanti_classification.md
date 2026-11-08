# Add SQANTI3 Structural Categories to an Annotation

Long-read transcript catalogs (IsoSeq, FLAIR, Bambu, StringTie) carry no
biotype for novel transcripts. SQANTI3 classifies every isoform against
a reference (`*_classification.txt`: `isoform`, `structural_category`,
`associated_gene`, `associated_transcript`, `subcategory`, ...).
`add_sqanti_classification()` writes those categories into the
`transcript_type` column of an annotation directory (and the associated
reference transcript / gene into `source_transcript_id`,
`source_gene_id` and `gene_name` when they are missing), so the SQANTI
categories behave like transcript types everywhere in the package
(plots, enrichment by category,
[`tx_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md)).
The categories are part of
[`tx_type_palette()`](https://luciorq.github.io/isoformic/reference/tx_type_palette.md).

## Usage

``` r
add_sqanti_classification(x, classification, overwrite_types = TRUE)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  (its annotation directory is modified in place and caches refresh
  automatically).

- classification:

  Path to the SQANTI3 classification file (tab separated, optionally
  gzipped) or a data frame with at least the `isoform` and
  `structural_category` columns.

- overwrite_types:

  Logical. Replace existing transcript types (default `TRUE`); with
  `FALSE` only missing types are filled.

## Value

The annotation directory path, invisibly.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "sqanti_demo")
#> ℹ Parsing example.gff3 (GFF) into /tmp/RtmpV2WK1u/sqanti_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
classification <- data.frame(
  isoform = tx$transcript_id[1:3],
  structural_category = c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog"),
  associated_gene = "RBPJ",
  associated_transcript = c(tx$transcript_id[1], "novel", "novel")
)
add_sqanti_classification(annot, classification)
#> ✔ Assigned SQANTI categories to 3 of 72 transcripts in
#>   /tmp/RtmpV2WK1u/sqanti_demo.
tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
head(tx[, c("transcript_id", "transcript_type")], 4)
#> # A tibble: 4 × 2
#>   transcript_id     transcript_type     
#>   <chr>             <chr>               
#> 1 ENST00000282397.9 full-splice_match   
#> 2 ENST00000615611.4 novel_in_catalog    
#> 3 ENST00000617835.4 novel_not_in_catalog
#> 4 ENST00000543394.2 protein_coding      
```
