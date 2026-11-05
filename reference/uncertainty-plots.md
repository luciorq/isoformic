# Plot the uncertainty of isoform-level results

- `plot_de_confidence(x, gene_name)` shows, for the isoforms of the
  given genes, the median log2 fold-change over the Salmon replicates
  with its 95% interval, coloured by the confidence class of
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  and with the shape marking the calls of the primary differential
  expression test. An interval that crosses zero, a faded colour or an
  ambiguous isoform tell the reader not to trust the call.

- `plot_isoform_correlation(x, gene_name)` draws the correlation between
  the replicates of every pair of isoforms of a gene. Strongly negative
  cells (outlined when ambiguous) mark isoforms whose reads cannot be
  told apart.

- `plot_inf_rv(x)` shows the distribution of the inferential relative
  variance per transcript type with the cutoff used to call isoforms
  uncertain.

- `confidence_palette()` gives the colours of the confidence classes.

## Usage

``` r
confidence_palette()

plot_de_confidence(self, gene_name, ncol = 1L, ...)

plot_isoform_correlation(self, gene_name, ncol = 1L, ...)

plot_inf_rv(self, top_types = 6L, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  with uncertainty tables (see
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)).

- gene_name:

  Gene names or IDs to plot.

- ncol:

  Number of panel columns.

- ...:

  Unused.

- top_types:

  Number of the most frequent transcript types to show.

## Value

A `ggplot` object, or a named character vector for
`confidence_palette()`.

## See also

[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
