# Plot Transcripts Genomic Context

Generate a genomic context plot for a specified gene, displaying its
transcripts along with their types and annotations. The function
utilizes the `plotgardener` package to create a detailed visualization
of the genomic context, including an ideogram, chromosome highlight, and
transcript structures. It requires a `ContextData` object that contains
the necessary genomic information and annotations. The plot can be
customized with various parameters such as offsets, label limits, and
ideogram references.

## Usage

``` r
plot_genomic_context(
  gene_name,
  context_data,
  limit_label = TRUE,
  show_guides = FALSE,
  y_offset = 0,
  height_offset = 0,
  downstream_offset = 0,
  upstream_offset = 0,
  ideogram_reference = c("hg38", "hg19", "mm11", "mm10", "none")
)
```

## Arguments

- gene_name:

  Character string specifying the name of the gene to plot.

- context_data:

  A `ContextData` object containing genomic context information.

- limit_label:

  Logical indicating whether to limit the length of transcript labels to
  avoid overlap (default is `TRUE`).

- show_guides:

  Logical indicating whether to show guide lines on the plot (default is
  `FALSE`).

- y_offset:

  Numeric value to adjust the vertical position of the plot (default is
  `0`).

- height_offset:

  Numeric value to adjust the height of the plot (default is `0`).

- downstream_offset:

  Numeric value to extend the downstream region beyond the gene's end
  position (default is `0`).

- upstream_offset:

  Numeric value to extend the upstream region beyond the gene's start
  position (default is `0`).

- ideogram_reference:

  Character string specifying the reference genome for the ideogram.
  Options include "hg38", "hg19", "mm11", "mm10", or "none" (default is
  "hg38").

## Value

A `plotgardener` object representing the genomic context plot.
