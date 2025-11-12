# Annotate Transcripts with Differential Gene Expression Significance

Adds a column to a transcript-level differential expression table
indicating whether each transcript originates from a gene that is
significantly differentially expressed.

## Usage

``` r
is_deg_sig(DegsigVector, DET_table)
```

## Arguments

- DegsigVector:

  A character vector containing the names of transcripts from
  significantly differentially expressed genes.

- DET_table:

  A `data.frame` or `tibble` containing transcript-level differential
  expression results, including a `transcript_name` column.

## Value

A `tibble` with an additional column `DEG_sig` indicating whether the
transcript is from a significantly differentially expressed gene
(`"YES"` or `"NO"`).

## Examples

``` r
# Sample data
significant_transcripts <- c("transcript1", "transcript3")
DET_table <- data.frame(
  transcript_name = c("transcript1", "transcript2", "transcript3", "transcript4"),
  log2FC = c(2.5, -1.2, 0.8, -0.5),
  pvalue = c(0.01, 0.2, 0.03, 0.6)
)

# Annotate transcripts with DEG significance
DET_table_annotated <- is_deg_sig(DegsigVector = significant_transcripts, DET_table = DET_table)

# View the result
print(DET_table_annotated)
#>   transcript_name log2FC pvalue DEG_sig
#> 1     transcript1    2.5   0.01     YES
#> 2     transcript3    0.8   0.03     YES
#> 3     transcript2   -1.2   0.20      NO
#> 4     transcript4   -0.5   0.60      NO
```
