# Merge Gene and Transcript Level Differential Expression Tables

Combines gene-level and transcript-level differential expression results
into a single table, annotates the combined data with significance
labels based on specified cutoffs, and filters transcripts based on
their types.

## Usage

``` r
join_DEG_DET(DEG_tab, DET_final_tab, logfc_cut, pval_cut)
```

## Arguments

- DEG_tab:

  A `data.frame` or `tibble` containing gene-level differential
  expression results, including `gene_id`, `gene_name`, `log2FC`, and
  `pvalue` columns.

- DET_final_tab:

  A `data.frame` or `tibble` containing transcript-level differential
  expression results, including `transcript_id`, `transcript_name`,
  `transcript_type`, `log2FC`, and `pvalue` columns.

- logfc_cut:

  A numeric value specifying the absolute log2 fold-change cutoff for
  significance.

- pval_cut:

  A numeric value specifying the p-value cutoff for significance.

## Value

A `tibble` combining gene and transcript differential expression
results, with additional columns:

- `id`: gene or transcript ID.

- `name`: gene or transcript name.

- `transcript_type`: type of transcript or `"gene"` for gene-level
  entries.

- `abs_log2FC`: absolute value of log2 fold-change.

- `significance`: `"sig"` if significant based on cutoffs, `"not_sig"`
  otherwise.

## Examples

``` r
# Sample gene-level data
DEG_tab <- data.frame(
  gene_id = c("gene1", "gene2"),
  gene_name = c("GeneA", "GeneB"),
  log2FC = c(1.5, -2.0),
  pvalue = c(0.01, 0.04)
)

# Sample transcript-level data
DET_final_tab <- data.frame(
  transcript_id = c("tx1", "tx2", "tx3"),
  transcript_name = c("Transcript1", "Transcript2", "Transcript3"),
  transcript_type = c("protein_coding", "lncRNA", "processed_transcript"),
  log2FC = c(1.2, -1.8, 0.5),
  pvalue = c(0.02, 0.03, 0.2)
)

# Merge and annotate differential expression results
DEGs_DETs_table <- join_DEG_DET(
  DEG_tab = DEG_tab,
  DET_final_tab = DET_final_tab,
  logfc_cut = 1,
  pval_cut = 0.05
)

# View the result
print(DEGs_DETs_table)
#>   feature_id feature_name log2FC pvalue         feature_type significance is_de
#> 1      gene1        GeneA    1.5   0.01                 gene          sig   yes
#> 2      gene2        GeneB   -2.0   0.04                 gene          sig   yes
#> 3        tx1  Transcript1    1.2   0.02       protein_coding          sig   yes
#> 4        tx2  Transcript2   -1.8   0.03               lncRNA          sig   yes
#> 5        tx3  Transcript3    0.5   0.20 processed_transcript      not_sig    no
```
