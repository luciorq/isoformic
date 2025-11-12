# Prepare Data for Gene and Transcript Expression Profile Plot

This function processes gene and transcript-level expression data, along
with differential expression results, to prepare a tidy data frame
suitable for plotting expression profiles across different sample
groups.

## Usage

``` r
prepare_profile_data(
  txi_gene = NULL,
  txi_transcript,
  sample_metadata,
  tx_to_gene,
  de_result_gene,
  de_result_transcript,
  var,
  var_levels,
  gene_col = "gene_name",
  tx_col = "transcript_name",
  pvalue_cutoff = 0.05,
  lfc_cutoff = 1,
  use_fdr = TRUE
)
```

## Arguments

- txi_gene:

  A `tibble` or `tximport` output containing gene-level expression
  abundances. If `NULL`, gene-level abundances will be summarized from
  `txi_transcript`. Default is `NULL`.

- txi_transcript:

  A `tibble` or `tximport` output containing transcript-level expression
  abundances.

- sample_metadata:

  A `data.frame` or `tibble` containing sample metadata. The first
  column should contain sample names matching the column names in
  `txi_gene` and `txi_transcript`.

- tx_to_gene:

  A `data.frame` or `tibble` containing transcript-to-gene mapping
  information. Must include columns specified by `gene_col` and
  `tx_col`.

- de_result_gene:

  A `data.frame` or `tibble` containing differential expression results
  at the gene level. Must include `gene_name`, `log2FC`, and `qvalue`
  columns.

- de_result_transcript:

  A `data.frame` or `tibble` containing differential expression results
  at the transcript level. Must include `transcript_name`, `log2FC`, and
  `qvalue` columns.

- var:

  A string specifying the column name in `sample_metadata` that
  indicates the grouping variable (e.g., treatment, condition).

- var_levels:

  A character vector specifying the levels of `var` to include in the
  contrasts.

- gene_col:

  A string specifying the column name in `tx_to_gene` that contains gene
  names. Default is `"gene_name"`.

- tx_col:

  A string specifying the column name in `tx_to_gene` that contains
  transcript names. Default is `"transcript_name"`.

- pvalue_cutoff:

  A numeric value specifying the p-value cutoff for determining
  significant differential expression. Default is `0.05`.

- lfc_cutoff:

  A numeric value specifying the log2 fold-change cutoff for determining
  significant differential expression. Default is `1`.

- use_fdr:

  A logical value indicating whether to use the false discovery rate
  (`qvalue`) instead of p-value for significance cutoff. Default is
  `TRUE`.

## Value

A `tibble` containing processed expression data and differential
expression flags, ready for plotting.

## Details

The function combines gene and transcript expression data with
differential expression results to generate a tidy data frame. It
filters significant genes and transcripts based on specified cutoffs and
prepares the data for plotting expression profiles across specified
sample groups.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming txi_gene, txi_transcript, sample_metadata, tx_to_gene, de_result_gene,
# and de_result_transcript are pre-loaded data frames:

# Prepare data for plotting
if (FALSE) {
  expr_df <- prepare_profile_data(
    txi_gene = txi_gene,
    txi_transcript = txi_transcript,
    sample_metadata = sample_metadata,
    tx_to_gene = tx_to_gene,
    de_result_gene = de_result_gene,
    de_result_transcript = de_result_transcript,
    var = "condition",
    var_levels = c("control", "treatment"),
    gene_col = "gene_name",
    tx_col = "transcript_name",
    pvalue_cutoff = 0.05,
    lfc_cutoff = 1,
    use_fdr = TRUE
  )

  # View the prepared data
  utils::head(expr_df)

  # Plotting example (assuming ggplot2 is installed)
  library(ggplot2)
  ggplot(expr_df, aes(x = condition, y = mean_TPM, fill = DE)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ parent_gene + transcript_type)
}
} # }
```
