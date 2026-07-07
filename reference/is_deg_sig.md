# Annotate Transcripts with Differential Gene Expression Significance

Adds a column to a transcript-level differential expression table
indicating whether each transcript originates from a gene that is
significantly differentially expressed.

## Usage

``` r
is_deg_sig(deg_sig_vector, det_table)
```

## Arguments

- deg_sig_vector:

  A character vector containing the names of transcripts from
  significantly differentially expressed genes.

- det_table:

  A `data.frame` or `tibble` containing transcript-level differential
  expression results, including a `transcript_name` column.

## Value

A `tibble` with an additional column `DEG_sig` indicating whether the
transcript is from a significantly differentially expressed gene
(`"YES"` or `"NO"`). Row order matches the input `det_table`.

## Examples

``` r
# Sample data
significant_transcripts <- c("transcript1", "transcript3")
det_table <- data.frame(
  transcript_name = c("transcript1", "transcript2", "transcript3", "transcript4"),
  log2FC = c(2.5, -1.2, 0.8, -0.5),
  pvalue = c(0.01, 0.2, 0.03, 0.6)
)

# Annotate transcripts with DEG significance
det_table_annotated <- is_deg_sig(
  deg_sig_vector = significant_transcripts,
  det_table = det_table
)

# View the result
print(det_table_annotated)
#>   transcript_name log2FC pvalue DEG_sig
#> 1     transcript1    2.5   0.01     YES
#> 2     transcript2   -1.2   0.20      NO
#> 3     transcript3    0.8   0.03     YES
#> 4     transcript4   -0.5   0.60      NO
```
