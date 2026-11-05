# Poly(A) Sites and Alternative Polyadenylation of Transcripts

`polya_sites()` describes the 3' end of every assayed transcript:

- `tes_pos`, `tes_group` (alternative 3' ends of the gene, from
  [`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)),
  `tes_rank` (1 = most proximal 3' end along the direction of
  transcription), `n_tes_groups`, `is_proximal`, `is_distal`.

- `apa_type` relative to the canonical isoform: `"same_end"`,
  `"tandem_utr"` (different 3' end, same stop codon: 3' UTR shortening
  or lengthening), `"alternative_last_exon"` (different last exon),
  `"intronic"` (3' end inside an intron of the canonical isoform,
  typically a composite terminal exon) or `"other"`.

- `utr3_length` (from
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  when sequences are available, else from the annotated CDS).

- `polya_signal` and `polya_signal_distance`: the poly(A) signal hexamer
  (`AATAAA`, `ATTAAA` or one of the weaker variants) found in the last
  `signal_window` nt of the transcript sequence, and its distance to the
  3' end (canonical signals sit 10-30 nt upstream); `NA` without
  sequences.

- With `polya_atlas`, the nearest atlas site within `window` bp of the
  3' end: `pas_site_id`, `pas_class`, `pas_tpm`, `pas_distance`,
  `pas_supported`. The atlas is the frozen PolyASite 2.0 table
  (`polya_atlas = "polyasite"`, see
  [frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md)),
  a data frame with `seqid`, `cleavage_pos`, `strand` (and optionally
  `site_id`, `site_class`, `tpm`), or a GENCODE `polyAs.gff3.gz` file
  (`download_reference(file_type = "polya")`).

`apa_usage()` quantifies per gene, for genes with at least two 3' end
groups: the share of expression ending at the distal-most 3' end
(`distal_fraction`) and the abundance-weighted 3' UTR length
(`weighted_utr3_length`), per sample or per group with the change
between the last and the first group level. `test_apa()` tests the
change in distal usage between two groups with the quasi-binomial test
used by
[`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md),
on counts.

## Usage

``` r
polya_sites(
  self,
  gene_name = NULL,
  polya_atlas = NULL,
  window = 50L,
  signal_window = 50L,
  ...
)

read_polya_gff(path)

apa_usage(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_expr = 0,
  as_matrix = FALSE,
  gene_name = NULL,
  ...
)

test_apa(
  self,
  group,
  levels = NULL,
  assay = "counts",
  min_expr = 0,
  gene_name = NULL,
  ...
)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- polya_atlas:

  `NULL` (no atlas), `"polyasite"`, a data frame or a file path (Parquet
  table or GENCODE poly(A) GFF3).

- window:

  Distance in bp within which an atlas site supports a 3' end (default
  50).

- signal_window:

  Length of the 3' terminal sequence searched for a poly(A) signal
  (default 50 nt).

- ...:

  Passed on to methods.

- path:

  Path to a GENCODE `polyAs.gff3.gz` file.

- group:

  Name of a `col_data` column defining groups; `NULL` for per-sample
  values.

- levels:

  Order of the group levels; the change is `last - first`. For the tests
  ([`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md),
  `test_apa()`), the two groups to compare; samples of other groups are
  ignored.

- assay:

  Transcript assay to use (`"abundance"` for usage, `"counts"` for the
  test).

- min_expr:

  Minimum assay value for an isoform to count (default 0).

- as_matrix:

  Logical, return the genes x samples matrix of distal fractions instead
  of the summary table.

## Value

Tibbles with one row per transcript (`polya_sites()`) or per gene
(`apa_usage()`, `test_apa()`); a matrix for `as_matrix = TRUE`.

## Examples

``` r
iso <- example_isoformic()
pa <- polya_sites(iso, gene_name = "RBPJ")
#> Warning: ! 1 transcript has a sequence length differing from the exon length; its ORF is
#>   skipped.
pa[, c("transcript_name", "tes_rank", "apa_type", "utr3_length", "polya_signal")]
#> # A tibble: 10 × 5
#>    transcript_name tes_rank apa_type              utr3_length polya_signal
#>    <chr>              <int> <chr>                       <dbl> <chr>       
#>  1 RBPJ-222               1 alternative_last_exon          NA NA          
#>  2 RBPJ-201               6 same_end                     4121 NA          
#>  3 RBPJ-206               3 tandem_utr                      0 NA          
#>  4 RBPJ-225               2 alternative_last_exon         347 NA          
#>  5 RBPJ-204               4 tandem_utr                     59 AAAAAG      
#>  6 RBPJ-205               6 same_end                     4121 NA          
#>  7 RBPJ-214               5 tandem_utr                    176 NA          
#>  8 RBPJ-208               5 other                         357 NA          
#>  9 RBPJ-202               6 same_end                     4124 NA          
#> 10 RBPJ-210               5 tandem_utr                     NA NA          
head(apa_usage(iso, group = "condition"))
#> # A tibble: 6 × 13
#>   gene_id        gene_name n_tes_groups apa_types mean_gene_expr distal_fraction
#>   <chr>          <chr>            <int> <chr>              <dbl>           <dbl>
#> 1 ENSG000000094… REV3L                2 intronic           831.            0.647
#> 2 ENSG000001308… PNCK                 4 alternat…           50.6           0.594
#> 3 ENSG000001001… POLR2F               6 alternat…          350.            0.125
#> 4 ENSG000001005… PSMA3                6 alternat…          860.            0.686
#> 5 ENSG000000703… FGF22                2 intronic            48.2           0.364
#> 6 ENSG000000542… ARID4B               7 alternat…         1116.            0.703
#> # ℹ 7 more variables: weighted_utr3_length <dbl>,
#> #   distal_fraction_preeclampsia <dbl>, distal_fraction_control <dbl>,
#> #   delta_distal_fraction <dbl>, weighted_utr3_length_preeclampsia <dbl>,
#> #   weighted_utr3_length_control <dbl>, delta_weighted_utr3_length <dbl>
head(test_apa(iso, group = "condition"))
#> # A tibble: 6 × 12
#>   gene_id            gene_name distal_fraction_preeclam…¹ distal_fraction_cont…²
#>   <chr>              <chr>                          <dbl>                  <dbl>
#> 1 ENSG00000009413.16 REV3L                         0.752                  0.593 
#> 2 ENSG00000048828.17 FAM120A                       0.920                  0.961 
#> 3 ENSG00000100038.19 TOP3B                         0.925                  0.976 
#> 4 ENSG00000074319.13 TSG101                        0.0104                 0.0208
#> 5 ENSG00000146648.18 EGFR                          0.924                  0.990 
#> 6 ENSG00000100142.15 POLR2F                        0.223                  0.0834
#> # ℹ abbreviated names: ¹​distal_fraction_preeclampsia, ²​distal_fraction_control
#> # ℹ 8 more variables: delta_distal_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>
```
