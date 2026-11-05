# Example IsoformicExperiment

Loads a small
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
bundled with the package: transcript abundances (TPM) for the isoforms
of a handful of genes from a preeclampsia placenta RNA-seq study (25
samples), GENCODE v33 annotation for those genes, and transcript- and
gene-level differential expression results. It is used throughout the
documentation and is suitable for trying every method of the class in
well under a second.

## Usage

``` r
example_isoformic()
```

## Value

An
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Examples

``` r
iso <- example_isoformic()
iso
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 25 [condition]
#>   transcripts: 499 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (499 transcripts), deg (69 genes)
```
