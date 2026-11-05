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
object. The data-frame methods accept the exon table of an annotation
(see
[`exon_table()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md))
joined with the transcript columns.

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

  Passed on to methods.

## Value

A tibble with one row per transcript.

## Examples

``` r
iso <- example_isoformic()
feats <- transcript_features(iso, gene_name = "RBPJ")
feats[, c("transcript_name", "orf_length", "nmd_predicted", "tss_group", "n_skipped_exons")]
#> # A tibble: 10 × 5
#>    transcript_name orf_length nmd_predicted tss_group n_skipped_exons
#>    <chr>                <int> <lgl>             <int>           <int>
#>  1 RBPJ-222               283 FALSE                 1               2
#>  2 RBPJ-201              1503 FALSE                 2               0
#>  3 RBPJ-206              1503 FALSE                 3               0
#>  4 RBPJ-225                NA FALSE                 4               0
#>  5 RBPJ-204              1464 FALSE                 5               0
#>  6 RBPJ-205              1461 FALSE                 5               0
#>  7 RBPJ-214              1398 FALSE                 6               0
#>  8 RBPJ-208              1149 FALSE                 7               1
#>  9 RBPJ-202              1461 FALSE                 7               0
#> 10 RBPJ-210                NA FALSE                 8               0
splicing_events(iso, gene_name = "EGFR")
#> # A tibble: 7 × 12
#>   transcript_id     transcript_name gene_id     gene_name canonical_transcript…¹
#>   <chr>             <chr>           <chr>       <chr>     <chr>                 
#> 1 ENST00000275493.7 EGFR-201        ENSG000001… EGFR      ENST00000275493.7     
#> 2 ENST00000455089.5 EGFR-207        ENSG000001… EGFR      ENST00000275493.7     
#> 3 ENST00000344576.6 EGFR-203        ENSG000001… EGFR      ENST00000275493.7     
#> 4 ENST00000420316.6 EGFR-204        ENSG000001… EGFR      ENST00000275493.7     
#> 5 ENST00000459688.1 EGFR-208        ENSG000001… EGFR      ENST00000275493.7     
#> 6 ENST00000638463.1 EGFR-211        ENSG000001… EGFR      ENST00000275493.7     
#> 7 ENST00000485503.1 EGFR-210        ENSG000001… EGFR      ENST00000275493.7     
#> # ℹ abbreviated name: ¹​canonical_transcript_id
#> # ℹ 7 more variables: n_skipped_exons <int>, n_novel_exons <int>,
#> #   n_retained_introns <int>, n_alt_5ss <int>, n_alt_3ss <int>,
#> #   alt_first_exon <lgl>, alt_last_exon <lgl>
```
