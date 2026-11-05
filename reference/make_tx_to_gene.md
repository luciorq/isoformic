# Create a Transcript-to-Gene Mapping Table

Builds a transcript-to-gene table from a GENCODE transcriptome FASTA
file (parsing the `|`-separated headers) or from a GFF3/GTF annotation
(through
[`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)).
Parsing is done by DuckDB, so a full GENCODE FASTA (about 250,000
headers, gzip-compressed) takes a couple of seconds.

## Usage

``` r
make_tx_to_gene(
  file_path,
  file_type = c("auto", "fasta", "gff", "gtf"),
  output_path = NULL
)
```

## Arguments

- file_path:

  Path to a FASTA (`.fa`, `.fasta`, optionally `.gz`), GFF3 or GTF file,
  or to an annotation directory.

- file_type:

  `"auto"` (default) detects the format from the file extension;
  otherwise one of `"fasta"`, `"gff"`, `"gtf"`.

- output_path:

  Passed to
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)
  for GFF3/GTF inputs.

## Value

A tibble with columns `transcript_id`, `gene_id`, `havana_gene_id`,
`havana_transcript_id`, `transcript_name`, `gene_name`, `tx_length` and
`transcript_type`. The Havana columns are `NA` for GFF3/GTF inputs.

## Details

With an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
object, use
[`tx_to_gene()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
or
[`tx_annot()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
which read the same information from the annotation directory.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
tx2gene <- make_tx_to_gene(gff, output_path = tempdir())
#> ℹ Parsing example.gff3 (GFF) into /tmp/Rtmp2Tt1au/example.
head(tx2gene)
#> # A tibble: 6 × 8
#>   transcript_id     gene_id  havana_gene_id havana_transcript_id transcript_name
#>   <chr>             <chr>    <chr>          <chr>                <chr>          
#> 1 ENST00000282397.9 ENSG000… NA             NA                   FLT1-201       
#> 2 ENST00000615611.4 ENSG000… NA             NA                   FLT1-206       
#> 3 ENST00000617835.4 ENSG000… NA             NA                   FLT1-208       
#> 4 ENST00000543394.2 ENSG000… NA             NA                   FLT1-205       
#> 5 ENST00000540678.2 ENSG000… NA             NA                   FLT1-203       
#> 6 ENST00000541932.5 ENSG000… NA             NA                   FLT1-204       
#> # ℹ 3 more variables: gene_name <chr>, tx_length <int>, transcript_type <chr>
```
