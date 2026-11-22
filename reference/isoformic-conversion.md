# Convert Between IsoformicExperiment and SummarizedExperiment

`as_isoformic()` builds an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
from a `SummarizedExperiment` (for example the output of `tximeta`),
from a `tximport` list (`abundance`, `counts`, `length`) or from a named
list of matrices. `as_summarized_experiment()` does the reverse at the
transcript or gene level, which is convenient for running Bioconductor
tools such as `DESeq2` or `fishpond` on the stored assays (a `counts`
assay is placed first, as those tools expect). The `metadata` of the
result holds `experiment_name`, `annot_path` and the object's `metadata`
without the inferential replicates.

## Usage

``` r
as_isoformic(
  x,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE
)

as_summarized_experiment(self, level = c("transcript", "gene"), ...)
```

## Arguments

- x:

  A `SummarizedExperiment`, a `tximport` list or a named list of
  matrices. An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  is returned as it is.

- annot_path:

  Annotation directory or GFF3/GTF file.

- col_data:

  Sample metadata; taken from `colData(x)` when `x` is a
  `SummarizedExperiment`.

- experiment_name:

  Optional experiment name.

- ignore_after_bar:

  Logical, truncate GENCODE `|`-separated transcript IDs.

- drop_unannotated:

  Logical, drop transcripts absent from the annotation instead of
  failing.

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- level:

  `"transcript"` (default) or `"gene"`.

- ...:

  Not used; reserved for methods.

## Value

`as_isoformic()` returns an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md);
`as_summarized_experiment()` returns a `SummarizedExperiment`.

## Examples

``` r
iso <- example_isoformic()
txi <- list(
  abundance = tx_assay(iso, "abundance"),
  counts = tx_assay(iso, "counts")
)
annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
iso2 <- as_isoformic(txi, annot_path = annot, col_data = col_data(iso))
iso2
#> <IsoformicExperiment> isoformic_f2faf1f9ff
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 25 [condition]
#>   transcripts: 1,911 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>

# Loading the SummarizedExperiment namespace takes several seconds
if (interactive() && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- as_summarized_experiment(iso, level = "gene")
  se
}
```
