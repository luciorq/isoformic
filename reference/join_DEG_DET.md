# Merge Gene and Transcript Level Differential Expression Tables

Table-based helper combining gene-level and transcript-level results
into the long format used by
[`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md).
With an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md),
[`de_table()`](https://luciorq.github.io/isoformic/reference/dea.md)
produces the same table from the stored results.

## Usage

``` r
join_DEG_DET(
  deg_table,
  det_table,
  logfc_cut = 1,
  pval_cut = 0.05,
  use_fdr = FALSE,
  transcript_types = default_transcript_types()
)
```

## Arguments

- deg_table:

  Gene-level results with `gene_id`, `gene_name`, `log2FC` and `pvalue`
  columns.

- det_table:

  Transcript-level results with `transcript_id`, `transcript_name`,
  `transcript_type`, `log2FC` and `pvalue` columns.

- logfc_cut:

  Absolute log2 fold-change cutoff for significance.

- pval_cut:

  P-value cutoff for significance. Like the cutoffs of
  [`de_tx()`](https://luciorq.github.io/isoformic/reference/dea.md), a
  feature is significant when `abs(log2FC) >= logfc_cut` and the p-value
  is `<= pval_cut`.

- use_fdr:

  Logical, test the `qvalue` column instead of `pvalue` (when present),
  as the stored cutoffs of an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  do.

- transcript_types:

  Transcript types to keep. `NULL` keeps all.

## Value

A tibble with `feature_id`, `feature_name`, `feature_type`, `gene_name`,
`log2FC`, `pvalue` (and `qvalue` when present) and `is_de`.

## Examples

``` r
deg_table <- data.frame(
  gene_id = c("gene1", "gene2"),
  gene_name = c("GeneA", "GeneB"),
  log2FC = c(1.5, -2.0),
  pvalue = c(0.01, 0.04)
)
det_table <- data.frame(
  transcript_id = c("tx1", "tx2", "tx3"),
  transcript_name = c("Transcript1", "Transcript2", "Transcript3"),
  transcript_type = c("protein_coding", "lncRNA", "processed_transcript"),
  gene_name = c("GeneA", "GeneA", "GeneB"),
  log2FC = c(1.2, -1.8, 0.5),
  pvalue = c(0.02, 0.03, 0.2)
)
join_DEG_DET(deg_table, det_table, logfc_cut = 1, pval_cut = 0.05)
#> # A tibble: 5 × 8
#>   feature_id feature_name feature_type     gene_id gene_name log2FC pvalue is_de
#>   <chr>      <chr>        <chr>            <chr>   <chr>      <dbl>  <dbl> <chr>
#> 1 gene1      GeneA        gene             gene1   GeneA        1.5   0.01 yes  
#> 2 gene2      GeneB        gene             gene2   GeneB       -2     0.04 yes  
#> 3 tx1        Transcript1  protein_coding   NA      GeneA        1.2   0.02 yes  
#> 4 tx2        Transcript2  lncRNA           NA      GeneA       -1.8   0.03 yes  
#> 5 tx3        Transcript3  processed_trans… NA      GeneB        0.5   0.2  no   
```
