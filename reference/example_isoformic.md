# Example IsoformicExperiment

Loads a small
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
bundled with the package: transcript abundances (TPM) and demonstration
counts (TPM x 20, rounded) for the isoforms of 69 genes from a
preeclampsia placenta RNA-seq study (25 samples), GENCODE v50 annotation
and transcript sequences for those genes, and transcript- and gene-level
differential expression results. It is used throughout the documentation
and is suitable for trying every method of the class in well under a
second.

## Usage

``` r
example_isoformic()
```

## Value

An
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
whose annotation directory is a private temporary copy of the bundled
one.

## Details

The annotation directory of the bundled experiment is copied on every
call into a new directory under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), and the returned
object points there (`iso@annot_path`). Functions that store tables in
the annotation directory, such as
[`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md),
[`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md),
[`add_appris()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md),
[`add_riboseq_orfs()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
or
[`add_nmd_exons()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md),
therefore change only that private copy: the installed package is never
modified, two calls never see each other's changes, and the copy is
removed with the temporary directory when the R session ends.
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
copies such a temporary annotation into an experiment directory outside
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), so an object saved
there can be read again in a later session.

## Examples

``` r
iso <- example_isoformic()
iso
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 25 [condition]
#>   transcripts: 1,911 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (710 transcripts), deg (69 genes)
```
