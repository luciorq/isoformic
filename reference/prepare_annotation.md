# Read a Prepared Annotation into Memory

Reads the Parquet tables of an annotation directory (or converts a
GFF3/GTF file first, see
[`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md))
and returns them as a list of tibbles. This is the table-based
counterpart of the `annot_data_*()` accessors of
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Usage

``` r
prepare_annotation(
  file_path,
  file_type = c("auto", "gff", "gtf"),
  output_path = NULL
)
```

## Arguments

- file_path:

  Path to an annotation directory or to a GFF3/GTF file.

- file_type:

  Passed to
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)
  when `file_path` is a file.

- output_path:

  Passed to
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md).

## Value

A named list with tibbles `genes`, `transcripts`, `exons`, `seq_regions`
and a named list `metadata`.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
annot <- prepare_annotation(gff, output_path = tempdir())
#> ℹ Using cached annotation directory /tmp/Rtmp2Tt1au/example.
names(annot)
#> [1] "genes"       "transcripts" "exons"       "seq_regions" "metadata"   
head(annot$genes)
#> # A tibble: 3 × 14
#>   gene_id   gene_name gene_type source_gene_id gene_class gene_category dhs_type
#>   <chr>     <chr>     <chr>     <chr>          <chr>      <chr>         <chr>   
#> 1 ENSG0000… FLT1      protein_… NA             NA         NA            NA      
#> 2 ENSG0000… RBPJ      protein_… NA             NA         NA            NA      
#> 3 ENSG0000… PNCK      protein_… NA             NA         NA            NA      
#> # ℹ 7 more variables: coding_status <chr>, cumulative_support <chr>,
#> #   haplotype <chr>, seqid <chr>, start_pos <int>, end_pos <int>, strand <chr>
```
