# Plot the uncertainty of isoform-level results

- `plot_de_confidence(x, gene_name)` shows, for the isoforms of the
  given genes, the median log2 fold-change over the inferential
  replicates with its 95% interval, coloured by whether the direction of
  the fold-change is stable over the replicates (`unstable` in the table
  that `add_uncertainty(group = )` stored, see
  [`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md);
  grey for isoforms not assessed) and with the shape marking the calls
  of the primary differential expression test. An unstable direction
  marks a call whose sign the quantification does not determine; in a
  placenta cohort such calls replicated less often. A stable direction
  does not make a call right: in control-versus-control comparisons of
  that cohort, 93% of the false calls had a stable direction.

- `plot_isoform_correlation(x, gene_name)` draws the correlation between
  the replicates of every pair of isoforms of a gene. Strongly negative
  cells (outlined when ambiguous) mark isoforms whose reads cannot be
  told apart.

- `plot_inf_rv(x)` shows the distribution of the inferential relative
  variance (log scale) of the assessed transcripts per transcript type,
  with the cutoff used to call isoforms uncertain.

## Usage

``` r
plot_de_confidence(self, gene_name, ncol = 1L, ...)

plot_isoform_correlation(self, gene_name, ncol = 1L, ...)

plot_inf_rv(self, top_types = 6L, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  with uncertainty tables (see
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md));
  `plot_de_confidence()` also needs the fold-change table of
  `add_uncertainty(group = )`.

- gene_name:

  Gene names or IDs to plot.

- ncol:

  Number of panel columns.

- ...:

  Not used; reserved for methods.

- top_types:

  Number of the most frequent transcript types to show.

## Value

A `ggplot` object.

## See also

[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
