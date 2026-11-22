# Structural and Coding Features of Transcripts

Computes, for every assayed transcript, features that go beyond the
biotype label, from the exon and CDS coordinates of the annotation:

- **Ends**: `tss_pos`, `tes_pos` (genomic positions of the transcription
  start and end), `tss_group` and `tes_group` (integer labels of
  alternative promoters / 3' ends within the gene; starts or ends closer
  than `end_window` bp share a group).

- **ORF**: `has_cds`, `cds_complete` (`FALSE` when the annotation flags
  an incomplete CDS with `cds_start_NF`/`cds_end_NF`), `orf_length`
  (nt), `utr5_length`, `utr3_length`.

- **NMD rule**: `stop_to_last_junction`, the distance in transcript
  coordinates from the stop codon to the last exon-exon junction, and
  `nmd_predicted` (`TRUE` when that distance exceeds `nmd_rule` nt, the
  classical 50-55 nt rule).

- **Comparison with the canonical isoform**: `canonical_transcript_id`
  (the MANE Select / Ensembl canonical transcript of the gene, else the
  transcript with the longest CDS, else the longest transcript),
  `orf_fraction_of_canonical` and `cds_differs_from_canonical`.

- **Splicing events relative to the canonical isoform** (see
  `splicing_events()`): `n_skipped_exons`, `n_novel_exons`,
  `n_retained_introns`, `n_alt_5ss`, `n_alt_3ss`, `alt_first_exon`,
  `alt_last_exon`.

`splicing_events()` returns the event counts alone; both are cached per
object. With an annotation directory path every transcript of the
annotation is described.

## Usage

``` r
transcript_features(
  self,
  gene_name = NULL,
  nmd_rule = 50,
  end_window = 50,
  ...
)

splicing_events(self, gene_name = NULL, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  or an annotation directory path.

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- nmd_rule:

  Distance in nucleotides from the stop codon to the last junction above
  which NMD is predicted (default 50).

- end_window:

  Distance in bp within which transcription starts (or ends) are grouped
  as the same promoter (or 3' end); default 50.

- ...:

  Not used; reserved for methods.

## Value

A tibble with one row per transcript. `transcript_features()`:
`transcript_id`, `transcript_name`, `transcript_type`, `gene_id`,
`gene_name`, `strand`, `n_exons`, `tx_length`, `is_canonical` and the
feature columns listed above. `splicing_events()`: `transcript_id`,
`transcript_name`, `gene_id`, `gene_name`, `canonical_transcript_id` and
the event columns.

## Examples

``` r
iso <- example_isoformic()
feats <- transcript_features(iso, gene_name = "RBPJ")
feats[, c("transcript_name", "orf_length", "nmd_predicted", "tss_group", "n_skipped_exons")]
#> # A tibble: 65 × 5
#>    transcript_name orf_length nmd_predicted tss_group n_skipped_exons
#>    <chr>                <int> <lgl>             <int>           <int>
#>  1 RBPJ-222               283 FALSE                 1               1
#>  2 RBPJ-247               177 TRUE                  1               1
#>  3 RBPJ-243                NA FALSE                 2               1
#>  4 RBPJ-201              1503 FALSE                 3               1
#>  5 RBPJ-213               291 FALSE                 3               1
#>  6 RBPJ-223              1503 FALSE                 3               1
#>  7 RBPJ-203              1398 FALSE                 3               1
#>  8 RBPJ-206              1503 FALSE                 4               1
#>  9 RBPJ-225                NA FALSE                 5               1
#> 10 RBPJ-229               145 FALSE                 6               1
#> # ℹ 55 more rows
splicing_events(iso, gene_name = "EGFR")
#> # A tibble: 23 × 12
#>    transcript_id     transcript_name gene_id    gene_name canonical_transcript…¹
#>    <chr>             <chr>           <chr>      <chr>     <chr>                 
#>  1 ENST00001124896.1 EGFR-222        ENSG00000… EGFR      ENST00001134423.1     
#>  2 ENST00001118378.1 EGFR-221        ENSG00000… EGFR      ENST00001134423.1     
#>  3 ENST00001112636.1 EGFR-220        ENSG00000… EGFR      ENST00001134423.1     
#>  4 ENST00000420316.7 EGFR-204        ENSG00000… EGFR      ENST00001134423.1     
#>  5 ENST00000342916.8 EGFR-202        ENSG00000… EGFR      ENST00001134423.1     
#>  6 ENST00000344576.8 EGFR-203        ENSG00000… EGFR      ENST00001134423.1     
#>  7 ENST00001110465.1 EGFR-219        ENSG00000… EGFR      ENST00001134423.1     
#>  8 ENST00000455089.6 EGFR-206        ENSG00000… EGFR      ENST00001134423.1     
#>  9 ENST00000898199.2 EGFR-214        ENSG00000… EGFR      ENST00001134423.1     
#> 10 ENST00000898200.2 EGFR-215        ENSG00000… EGFR      ENST00001134423.1     
#> # ℹ 13 more rows
#> # ℹ abbreviated name: ¹​canonical_transcript_id
#> # ℹ 7 more variables: n_skipped_exons <int>, n_novel_exons <int>,
#> #   n_retained_introns <int>, n_alt_5ss <int>, n_alt_3ss <int>,
#> #   alt_first_exon <lgl>, alt_last_exon <lgl>
```
