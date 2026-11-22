# Transcriptomic Diversity Summaries

These generics quantify how expression is distributed across the
isoforms of each gene, which is the question an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
is built to answer. All computations use matrix operations on the cached
assays and scale to full transcriptome experiments.

- `isoform_fraction(x)`: matrix (transcripts x samples) of the fraction
  of its gene's expression carried by each transcript. `NA` where the
  gene is not expressed in a sample.

- `tx_diversity(x)`: one row per gene summarising isoform usage: number
  of annotated, assayed and expressed isoforms, number of transcript
  types, the dominant isoform and its mean fraction, and the *effective
  number of isoforms* (`exp(Shannon entropy)` of the mean isoform
  fractions, 1 = a single isoform carries all expression).

- `dominant_isoform(x, group)`: the dominant isoform of each gene,
  overall or per level of a `col_data` column, flagging genes whose
  dominant isoform changes between groups (`dominant_switch`).

- `isoform_switches(x)`: relates transcript- and gene-level differential
  expression per gene: counts of significant isoforms up and down, the
  gene-level call and a `switch_type` classification
  (`"opposite_direction"`, `"det_not_deg"`, `"consistent"`).

- `exon_diversity(x, gene_name)`: per gene, the number of distinct
  exons, how many are constitutive (shared by all isoforms) versus
  alternative, and the mean number of exons per isoform.

- `haplotype_summary(x, group)`: for haplotype-resolved annotations (see
  [`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md)),
  expression of each gene per haplotype, the haplotype fraction (allelic
  imbalance), the dominant isoform of each haplotype (compared through
  `source_transcript_id`) and whether the dominant isoform differs
  between haplotypes (`haplotype_isoform_switch`). Genes are matched
  across haplotypes by `source_gene_id`, falling back to `gene_name`.

`dominant_isoform()` also reports `dominant_is_canonical`, whether the
dominant isoform is the MANE Select / Ensembl canonical / RefSeq Select
transcript of its gene. When the object carries uncertainty tables (see
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
it adds `dominant_inf_rv`, the inferential relative variance of the
dominant isoform, and `dominant_ambiguous`, whether it cannot be told
apart from a sibling isoform. `isoform_switches()` adds
`n_det_unstable`, the number of significant isoforms whose fold-change
direction the inferential replicates do not agree on (`unstable` in
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)),
once those results are stored.

`drop_unstable = TRUE` leaves these calls out; calls whose direction was
not assessed (below `min_count`) are kept. When ambiguous isoforms were
summed (`add_uncertainty(group = , collapse_groups = TRUE)`), each
member of a group takes the direction of the group: the significant
members of an unstable group are counted in `n_det_unstable` and left
out. It is a sensitivity analysis. In a placenta cohort it would have
removed about 3% of the calls that replicated and about 10% of those
that did not. It does not remove false calls caused by biological
variability: most false calls have a stable direction.

## Usage

``` r
isoform_fraction(self, assay = "abundance", ...)

tx_diversity(
  self,
  assay = "abundance",
  min_expr = 1,
  min_samples = 1L,
  gene_name = NULL,
  ...
)

dominant_isoform(
  self,
  group = NULL,
  assay = "abundance",
  gene_name = NULL,
  min_expr = 0,
  levels = NULL,
  ...
)

isoform_switches(
  self,
  gene_name = NULL,
  confident_only = NULL,
  drop_unstable = FALSE,
  ...
)

exon_diversity(self, gene_name = NULL, ...)

haplotype_summary(
  self,
  group = NULL,
  assay = "abundance",
  gene_name = NULL,
  levels = NULL,
  ...
)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- assay:

  Name of the transcript assay to use (default `"abundance"`, falling
  back to the first assay).

- ...:

  Not used; reserved for methods.

- min_expr:

  Minimum assay value for a transcript to count as expressed in a sample
  (default `1` in `tx_diversity()`). In `dominant_isoform()` (default
  `0`, no floor), values below it are set to zero, so these isoforms
  neither dominate nor count towards the gene total.

- min_samples:

  Minimum number of samples in which a transcript must be expressed to
  count as expressed in `tx_diversity()`.

- gene_name:

  Optional character vector of gene names or IDs to restrict the result
  to.

- group:

  Name of a `col_data` column defining sample groups, or `NULL` for all
  samples together.

- levels:

  Levels of `group` to use, in this order. Samples of other groups are
  ignored, as if the object had been subset to these samples. `NULL`
  takes all the levels of `group`.

- confident_only:

  Defunct: the confidence class it filtered on was removed. Use
  `drop_unstable`, by name.

- drop_unstable:

  Logical. In `isoform_switches()`, leave out the significant isoforms
  with an unstable direction (`unstable` in
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)).

## Value

A matrix for `isoform_fraction()`, tibbles otherwise.

## Examples

``` r
iso <- example_isoformic()
isoform_fraction(iso)[1:3, 1:3]
#>                   SRR11498039 SRR11498040 SRR11498041
#> ENST00001127150.1   0.0000000   0.0000000   0.0000000
#> ENST00001127151.1   0.4706061   0.5451462   0.2777801
#> ENST00001091268.1   0.0000000   0.0000000   0.0000000
head(tx_diversity(iso))
#> # A tibble: 6 × 14
#>   gene_id         gene_name gene_type n_tx_annotated n_tx_assayed n_tx_expressed
#>   <chr>           <chr>     <chr>              <int>        <int>          <int>
#> 1 ENSG0000013082… PNCK      protein_…             59           59             15
#> 2 ENSG0000020484… ATXN2     protein_…            101          101             46
#> 3 ENSG0000015309… ACOXL     protein_…             49           49             19
#> 4 ENSG0000020458… DDR1      protein_…            135          135             24
#> 5 ENSG0000016043… RDH13     protein_…             28           27             24
#> 6 ENSG0000017697… DPP7      protein_…             65           65             23
#> # ℹ 8 more variables: n_tx_types <int>, tx_types <chr>, mean_gene_expr <dbl>,
#> #   dominant_transcript_id <chr>, dominant_transcript_name <chr>,
#> #   dominant_transcript_type <chr>, dominant_fraction <dbl>,
#> #   effective_isoforms <dbl>
dominant_isoform(iso, group = "condition")
#> # A tibble: 138 × 11
#>    gene_id            gene_name gene_type      condition  dominant_transcript_id
#>    <chr>              <chr>     <chr>          <chr>      <chr>                 
#>  1 ENSG00000092199.21 HNRNPC    protein_coding preeclamp… ENST00000553614.5     
#>  2 ENSG00000092199.21 HNRNPC    protein_coding control    ENST00000557157.5     
#>  3 ENSG00000100567.15 PSMA3     protein_coding preeclamp… ENST00000216455.9     
#>  4 ENSG00000100567.15 PSMA3     protein_coding control    ENST00000557508.6     
#>  5 ENSG00000104827.13 CGB3      protein_coding preeclamp… ENST00001140563.1     
#>  6 ENSG00000104827.13 CGB3      protein_coding control    ENST00000357383.4     
#>  7 ENSG00000105993.17 DNAJB6    protein_coding preeclamp… ENST00000488001.5     
#>  8 ENSG00000105993.17 DNAJB6    protein_coding control    ENST00000262177.9     
#>  9 ENSG00000107819.17 SFXN3     protein_coding preeclamp… ENST00000465383.5     
#> 10 ENSG00000107819.17 SFXN3     protein_coding control    ENST00000896242.2     
#> # ℹ 128 more rows
#> # ℹ 6 more variables: dominant_transcript_name <chr>,
#> #   dominant_transcript_type <chr>, dominant_is_canonical <lgl>,
#> #   dominant_fraction <dbl>, gene_expr <dbl>, dominant_switch <lgl>
isoform_switches(iso)
#> # A tibble: 60 × 12
#>    gene_id         gene_name gene_type n_tx_tested n_det_sig n_det_up n_det_down
#>    <chr>           <chr>     <chr>           <int>     <int>    <int>      <int>
#>  1 ENSG0000013311… TPT1      protein_…          47        10        5          5
#>  2 ENSG0000014664… EGFR      protein_…          20         7        5          2
#>  3 ENSG0000020484… ATXN2     protein_…          38         7        1          6
#>  4 ENSG0000010513… SYDE1     protein_…           9         6        5          1
#>  5 ENSG0000009219… HNRNPC    protein_…          28         5        4          1
#>  6 ENSG0000015442… CCSAP     protein_…           6         3        2          1
#>  7 ENSG0000010056… PSMA3     protein_…          11         2        1          1
#>  8 ENSG0000014339… PIP5K1A   protein_…          25         4        0          4
#>  9 ENSG0000012951… SNX6      protein_…          21         2        2          0
#> 10 ENSG0000012494… AHNAK     protein_…           7         1        1          0
#> # ℹ 50 more rows
#> # ℹ 5 more variables: det_sig_types <chr>, max_abs_log2FC <dbl>,
#> #   gene_log2FC <dbl>, deg_sig <chr>, switch_type <chr>
exon_diversity(iso, gene_name = c("RBPJ", "EGFR"))
#> # A tibble: 2 × 9
#>   gene_id   gene_name gene_type  n_tx mean_exons_per_tx max_exons_per_tx n_exons
#>   <chr>     <chr>     <chr>     <int>             <dbl>            <int>   <int>
#> 1 ENSG0000… RBPJ      protein_…    65              9.52               13     113
#> 2 ENSG0000… EGFR      protein_…    23             15.1                28      57
#> # ℹ 2 more variables: n_constitutive <int>, n_alternative <int>
```
