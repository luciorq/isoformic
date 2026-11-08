# Read Salmon Quantifications Into an IsoformicExperiment

Reads `quant.sf` files produced by Salmon (`salmon quant`) directly with
DuckDB, without going through `tximport`, and returns an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
with transcript-level `abundance` (TPM), `counts` (`NumReads`) and
`length` (`EffectiveLength`) assays. Gene-level assays are derived
lazily by the object (see
[`gene_assays()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)).

## Usage

``` r
read_salmon(
  paths,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE,
  inferential = FALSE
)
```

## Arguments

- paths:

  Character vector of Salmon output directories (containing `quant.sf`)
  or paths to `quant.sf` files. Names, when present, are used as sample
  IDs; otherwise the directory name is used.

- annot_path:

  Annotation directory or GFF3/GTF file (see
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)).

- col_data:

  Optional data frame of sample metadata whose first column (or a
  `sample_id` column) matches the sample IDs. When `NULL`, a minimal
  table with just `sample_id` is created.

- experiment_name:

  Optional experiment name.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- drop_unannotated:

  Logical. Drop transcripts absent from the annotation (with a message)
  instead of failing.

- inferential:

  Logical. Read the inferential replicates, equivalence classes and run
  metadata of the Salmon directories and compute the uncertainty tables
  with
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  (default settings; call that function for other cutoffs or a
  contrast). The directories must have been quantified with
  `--numGibbsSamples` or `--numBootstraps`.

## Value

An
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Details

All files are read in one DuckDB query and pivoted to wide matrices, so
reading tens of samples with 200,000 transcripts each takes seconds.
Inferential replicates (Gibbs/bootstrap samples) are read only with
`inferential = TRUE`, which also computes the uncertainty tables of
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md);
for `swish` they are imported with `tximeta`/`fishpond` and converted
with
[`as_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md).

## Examples

``` r
quant_dirs <- list.files(
  system.file("extdata", "mini_quants", package = "isoformic"),
  full.names = TRUE
)
names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
annot <- system.file("extdata", "mini_annotation", package = "isoformic")
iso <- read_salmon(quant_dirs, annot_path = annot)
iso
#> <IsoformicExperiment> isoformic_2da713f64b
#>   annotation : gencode_v49_mini_subset (15 genes, 111 transcripts, 1,071 exons)
#>   samples    : 6
#>   transcripts: 111 assayed from 15 genes
#>   assays     : abundance, counts, length
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
```
