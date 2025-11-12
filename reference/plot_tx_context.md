# Plot Transcript Genomic Context

This function plots the genomic context of all transcripts of given
genes.

## Usage

``` r
plot_tx_context(exon_table, custom_colors = NULL)
```

## Arguments

- exon_table:

  a tibble with exon information. Must contain columns `tx_id`,
  `exon_left`, and `exon_right`.

- custom_colors:

  a vector of colors to use for each transcript. If not provided, the
  function will use the default colors. Actually, this argument is
  \***NOT implemented** yet.
