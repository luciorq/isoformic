# Read Oarfish single-cell quantifications

With `--single-cell`, Oarfish quantifies a transcriptome-aligned BAM
whose reads carry a cell barcode (`CB:Z` tag) and writes, next to its
output prefix, `<prefix>.count.mtx` (a Matrix Market table of cells by
transcripts), `<prefix>.barcodes.txt`, `<prefix>.features.txt` and
`<prefix>.meta_info.json`. It draws no replicates in this mode.

## Usage

``` r
read_oarfish_sc_counts(path, ignore_after_bar = TRUE)

read_oarfish_sc(
  path,
  cells,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  min_cells = 1L,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE
)
```

## Arguments

- path:

  The output prefix of Oarfish (the value of `-o`), the `.count.mtx`
  file, or a directory holding a single `.count.mtx` file.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- cells:

  A data frame with a `barcode` column and a `sample_id` column.

- annot_path:

  Annotation directory or GFF3/GTF file (see
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)).

- col_data:

  Optional sample metadata with a `sample_id` column; the number of
  cells and reads of each sample are added as `n_cells` and `n_reads`.

- experiment_name:

  Optional experiment name.

- min_cells:

  Samples with fewer cells are dropped.

- drop_unannotated:

  Logical. Drop transcripts absent from the annotation (with a message)
  instead of failing.

## Value

`read_oarfish_sc_counts()`: a tibble with `barcode`, `transcript_id` and
`count`. `read_oarfish_sc()`: an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Details

- `read_oarfish_sc_counts()` returns the non-zero counts as a long
  table, one row per cell and transcript.

- `read_oarfish_sc()` sums the cells of each sample (a pseudobulk) and
  returns an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  with `counts`, `length` (from the annotation) and `abundance` (TPM
  computed from both) assays, so the isoform analyses of the package
  apply to single-cell data at the level of samples or cell types.
  `cells` says which sample each barcode belongs to; barcodes absent
  from it are dropped. To compare cell types, give each combination of
  sample and cell type its own `sample_id` and describe it in
  `col_data`.

## Examples

``` r
dir <- system.file("extdata", "oarfish_sc", package = "isoformic")
counts <- read_oarfish_sc_counts(file.path(dir, "out"))
head(counts)
#> # A tibble: 6 × 3
#>   barcode          transcript_id     count
#>   <chr>            <chr>             <dbl>
#> 1 CGGCGCGAACCACCTC ENST00000512351.5  1   
#> 2 CGGCGCGAACCACCTC ENST00000345843.7  1.26
#> 3 CGGCGCGAACCACCTC ENST00000342295.5  5.86
#> 4 CGGCGCGAACCACCTC ENST00000514380.5  1.13
#> 5 CGGCGCGAACCACCTC ENST00000355476.7  5.63
#> 6 CGGCGCGAACCACCTC ENST00000511401.5  6   
cells <- read.delim(file.path(dir, "cells.tsv"))
iso <- read_oarfish_sc(
  file.path(dir, "out"),
  cells = cells,
  annot_path = system.file(
    "extdata", "example_experiment", "annotation",
    package = "isoformic"
  )
)
col_data(iso)
#> # A tibble: 3 × 3
#>   sample_id n_cells n_reads
#>   <chr>       <int>   <dbl>
#> 1 control_1      20   3000.
#> 2 control_2      20   3000.
#> 3 treated_1      20   3000.
```
