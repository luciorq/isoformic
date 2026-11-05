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
`n_det_confident`, the number of significant isoforms with high
confidence, once
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
results are stored.

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
  ...
)

isoform_switches(self, gene_name = NULL, confident_only = FALSE, ...)

exon_diversity(self, gene_name = NULL, ...)

haplotype_summary(
  self,
  group = NULL,
  assay = "abundance",
  gene_name = NULL,
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

  Passed on to methods.

- min_expr:

  Minimum assay value for a transcript to count as expressed in a
  sample.

- min_samples:

  Minimum number of samples in which a transcript must be expressed to
  count as expressed.

- gene_name:

  Optional character vector of gene names or IDs to restrict the result
  to.

- group:

  Name of a `col_data` column defining sample groups, or `NULL` for all
  samples together.

- confident_only:

  Logical. Count only the significant isoforms with high confidence (see
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
  in `isoform_switches()`.

## Value

A matrix for `isoform_fraction()`, tibbles otherwise.

## Examples

``` r
iso <- example_isoformic()
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
isoform_fraction(iso)[1:3, 1:3]
#>                   SRR11498039 SRR11498040 SRR11498041
#> ENST00000374514.8  0.96675753  0.98117418  0.98239081
#> ENST00000420982.5  0.00000000  0.00000000  0.00000000
#> ENST00000374502.7  0.02265948  0.01882582  0.01760919
head(tx_diversity(iso))
#> # A tibble: 6 × 14
#>   gene_id         gene_name gene_type n_tx_annotated n_tx_assayed n_tx_expressed
#>   <chr>           <chr>     <chr>              <int>        <int>          <int>
#> 1 ENSG0000002319… RNH1      protein_…             24           17             17
#> 2 ENSG0000010042… HDAC10    protein_…             22           16             16
#> 3 ENSG0000022980… XIST      lncRNA                30           18             18
#> 4 ENSG0000007980… DNM2      protein_…             22           15             15
#> 5 ENSG0000007231… SREBF1    protein_…             24           18             18
#> 6 ENSG0000005082… BCAR1     protein_…             22           15             15
#> # ℹ 8 more variables: n_tx_types <int>, tx_types <chr>, mean_gene_expr <dbl>,
#> #   dominant_transcript_id <chr>, dominant_transcript_name <chr>,
#> #   dominant_transcript_type <chr>, dominant_fraction <dbl>,
#> #   effective_isoforms <dbl>
dominant_isoform(iso, group = "condition")
#> # A tibble: 138 × 11
#>    gene_id            gene_name gene_type      condition  dominant_transcript_id
#>    <chr>              <chr>     <chr>          <chr>      <chr>                 
#>  1 ENSG00000009413.16 REV3L     protein_coding preeclamp… ENST00000368802.8     
#>  2 ENSG00000009413.16 REV3L     protein_coding control    ENST00000619481.1     
#>  3 ENSG00000023191.17 RNH1      protein_coding preeclamp… ENST00000356187.9     
#>  4 ENSG00000023191.17 RNH1      protein_coding control    ENST00000397604.7     
#>  5 ENSG00000072310.17 SREBF1    protein_coding preeclamp… ENST00000355815.8     
#>  6 ENSG00000072310.17 SREBF1    protein_coding control    ENST00000261646.10    
#>  7 ENSG00000082996.19 RNF13     protein_coding preeclamp… ENST00000392894.7     
#>  8 ENSG00000082996.19 RNF13     protein_coding control    ENST00000361785.10    
#>  9 ENSG00000100038.19 TOP3B     protein_coding preeclamp… ENST00000444502.5     
#> 10 ENSG00000100038.19 TOP3B     protein_coding control    ENST00000457270.5     
#> # ℹ 128 more rows
#> # ℹ 6 more variables: dominant_transcript_name <chr>,
#> #   dominant_transcript_type <chr>, dominant_is_canonical <lgl>,
#> #   dominant_fraction <dbl>, gene_expr <dbl>, dominant_switch <lgl>
isoform_switches(iso)
#> # A tibble: 34 × 12
#>    gene_id         gene_name gene_type n_tx_tested n_det_sig n_det_up n_det_down
#>    <chr>           <chr>     <chr>           <int>     <int>    <int>      <int>
#>  1 ENSG0000008299… RNF13     protein_…          11         4        1          3
#>  2 ENSG0000008915… GCN1      protein_…           6         4        3          1
#>  3 ENSG0000010056… PSMA3     protein_…           8         3        1          2
#>  4 ENSG0000006366… GPC1      protein_…           6         3        2          1
#>  5 ENSG0000014664… EGFR      protein_…           7         3        3          0
#>  6 ENSG0000004882… FAM120A   protein_…           5         3        3          0
#>  7 ENSG0000008928… FUS       protein_…          11         2        2          0
#>  8 ENSG0000003244… PNPLA6    protein_…          11         2        2          0
#>  9 ENSG0000001134… SYT7      protein_…           7         2        2          0
#> 10 ENSG0000007112… WDR1      protein_…           7         1        1          0
#> # ℹ 24 more rows
#> # ℹ 5 more variables: det_sig_types <chr>, max_abs_log2FC <dbl>,
#> #   gene_log2FC <dbl>, deg_sig <chr>, switch_type <chr>
exon_diversity(iso, gene_name = c("RBPJ", "EGFR"))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp2Tt1au/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> # A tibble: 2 × 9
#>   gene_id   gene_name gene_type  n_tx mean_exons_per_tx max_exons_per_tx n_exons
#>   <chr>     <chr>     <chr>     <int>             <dbl>            <int>   <int>
#> 1 ENSG0000… RBPJ      protein_…    31              6.52               12      65
#> 2 ENSG0000… EGFR      protein_…    11             14.4                28      46
#> # ℹ 2 more variables: n_constitutive <int>, n_alternative <int>
```
