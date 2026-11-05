# Flag Transcripts Whose Gene Is Differentially Expressed

Table-based helper adding a `DEG_sig` column (`"YES"`/`"NO"`) to a
transcript table. Row order is preserved. With an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md),
[`isoform_switches()`](https://luciorq.github.io/isoformic/reference/diversity.md)
relates transcript and gene significance directly.

## Usage

``` r
is_deg_sig(deg_sig_vector, det_table)
```

## Arguments

- deg_sig_vector:

  Character vector of transcript names (or IDs) belonging to
  significantly differentially expressed genes.

- det_table:

  Transcript-level table with a `transcript_name` (or `transcript_id`)
  column.

## Value

`det_table` with an extra `DEG_sig` column.

## Examples

``` r
det_table <- data.frame(
  transcript_name = c("tx1", "tx2", "tx3"),
  log2FC = c(2.5, -1.2, 0.8)
)
is_deg_sig(c("tx1", "tx3"), det_table)
#>   transcript_name log2FC DEG_sig
#> 1             tx1    2.5     YES
#> 2             tx2   -1.2      NO
#> 3             tx3    0.8     YES
```
