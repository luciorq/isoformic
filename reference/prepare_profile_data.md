# Prepare Expression Profile Data From Loose Tables

Table-based counterpart of
[`profile_data()`](https://luciorq.github.io/isoformic/reference/profiles.md):
combines a transcript expression table, sample metadata, a
transcript-to-gene table and differential expression results into the
long profile table consumed by
[`plot_tx_expr()`](https://luciorq.github.io/isoformic/reference/profiles.md).
Gene expression is the sum of the isoform values unless `txi_gene` is
supplied.

## Usage

``` r
prepare_profile_data(
  txi_transcript,
  sample_metadata,
  tx_to_gene,
  de_result_gene = NULL,
  de_result_transcript = NULL,
  var,
  var_levels = NULL,
  txi_gene = NULL,
  gene_name = NULL,
  pvalue_cutoff = 0.05,
  lfc_cutoff = 1,
  use_fdr = TRUE
)
```

## Arguments

- txi_transcript:

  Transcript expression: a matrix or data frame with transcript IDs as
  row names or in a `transcript_id` column, or a `tximport` list (its
  `abundance` element is used).

- sample_metadata:

  Data frame whose first column (or `sample_id` column) holds the sample
  IDs and that contains the `var` column.

- tx_to_gene:

  Data frame with `transcript_id`, `gene_id`, `transcript_name`,
  `gene_name` and `transcript_type` columns (see
  [`make_tx_to_gene()`](https://luciorq.github.io/isoformic/reference/make_tx_to_gene.md)).

- de_result_gene:

  Gene-level results with `gene_id` or `gene_name`, `log2FC` and
  `qvalue`/`pvalue` columns. May be `NULL`.

- de_result_transcript:

  Transcript-level results with `transcript_id` or `transcript_name`,
  `log2FC` and `qvalue`/`pvalue` columns. May be `NULL`.

- var:

  Name of the grouping column in `sample_metadata`.

- var_levels:

  Order of the group levels. `NULL` uses the order of appearance.

- txi_gene:

  Optional gene expression table (same formats as `txi_transcript`, with
  gene IDs).

- gene_name:

  Optional genes to keep (names or IDs); `NULL` keeps all.

- pvalue_cutoff, lfc_cutoff, use_fdr:

  Significance cutoffs.

## Value

A tibble in the format described in
[`profile_data()`](https://luciorq.github.io/isoformic/reference/profiles.md).

## Examples

``` r
iso <- example_isoformic()
prof <- prepare_profile_data(
  txi_transcript = assay(iso, "abundance"),
  sample_metadata = col_data(iso),
  tx_to_gene = tx_annot(iso),
  de_result_gene = de_gene(iso),
  de_result_transcript = de_tx(iso),
  var = "condition",
  gene_name = "RBPJ"
)
head(prof)
#> # A tibble: 6 × 10
#>   feature_id feature_name feature_type gene_id gene_name group mean_expr sd_expr
#>   <chr>      <chr>        <chr>        <chr>   <chr>     <fct>     <dbl>   <dbl>
#> 1 ENSG00000… RBPJ         gene         ENSG00… RBPJ      pree…   1297.    427.  
#> 2 ENSG00000… RBPJ         gene         ENSG00… RBPJ      cont…   1367.    195.  
#> 3 ENST00000… RBPJ-201     protein_cod… ENSG00… RBPJ      pree…     11.2    25.1 
#> 4 ENST00000… RBPJ-201     protein_cod… ENSG00… RBPJ      cont…      3.95    8.30
#> 5 ENST00000… RBPJ-202     protein_cod… ENSG00… RBPJ      pree…    838.    303.  
#> 6 ENST00000… RBPJ-202     protein_cod… ENSG00… RBPJ      cont…   1095.    169.  
#> # ℹ 2 more variables: n_samples <int>, is_de <chr>
```
