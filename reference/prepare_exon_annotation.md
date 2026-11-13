# Exon Table for Selected Genes From an Annotation

Table-based helper returning one row per exon for the transcripts of the
requested genes, joined with transcript and gene attributes. This is the
input format of the data-frame method of
[`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md).
With an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
use
[`exon_table()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
instead.

## Usage

``` r
prepare_exon_annotation(
  gene_name,
  file_path,
  file_type = c("auto", "gff", "gtf"),
  output_path = NULL
)
```

## Arguments

- gene_name:

  Character vector of gene names (or gene IDs).

- file_path:

  Path to an annotation directory or a GFF3/GTF file.

- file_type:

  Passed to
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md).

- output_path:

  Passed to
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md).

## Value

A tibble with columns `gene_id`, `gene_name`, `transcript_id`,
`transcript_name`, `transcript_type`, `exon_id`, `exon_number`, `seqid`,
`start_pos`, `end_pos`, `strand`.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
exons <- prepare_exon_annotation("RBPJ", gff, output_path = tempdir())
#> ℹ Using cached annotation directory /tmp/RtmpmIi7db/example.
head(exons)
#> # A tibble: 6 × 12
#>   gene_id      gene_name transcript_id transcript_name transcript_type haplotype
#>   <chr>        <chr>     <chr>         <chr>           <chr>           <chr>    
#> 1 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 2 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 3 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 4 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 5 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 6 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> # ℹ 6 more variables: exon_id <chr>, exon_number <int>, seqid <chr>,
#> #   start_pos <int>, end_pos <int>, strand <chr>
```
