# Download HPRC Per-Haplotype CAT Annotations

Downloads the Human Pangenome Reference Consortium year-1 CAT/Liftoff
gene annotations of an assembled sample, one GFF3 per haplotype, lifted
from GENCODE v38 onto the sample's own assembly in CHM13 or GRCh38
coordinates. The files keep GENCODE identifiers in `source_gene` /
`source_transcript`, so
[`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md)
and
[`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)
can compare haplotypes.

## Usage

``` r
download_hprc_annotation(
  sample,
  haplotype = c(1L, 2L),
  assembly = c("chm13", "hg38"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
)
```

## Arguments

- sample:

  HPRC sample ID, for example `"HG00438"` or `"HG002"`.

- haplotype:

  Haplotype numbers to download, `1`, `2` or `c(1, 2)` (default).

- assembly:

  Coordinate space of the CAT run: `"chm13"` (default) or `"hg38"`.

- output_path:

  Directory where the file is saved. Defaults to `":cache:"`, the
  package cache (see
  [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md)).

- timeout_limit:

  Maximum time in seconds for the download. Takes precedence over
  `options("timeout")`. Defaults to `3600`.

- method:

  Method used by
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).
  Defaults to `"auto"`.

## Value

A named character vector of file paths (names `hap1`, `hap2`), ready for
[`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md).

## Examples

``` r
if (FALSE) { # \dontrun{
files <- download_hprc_annotation("HG00438")
merged <- merge_annotations(files, annot_name = "HG00438_chm13")
merged
} # }
```
