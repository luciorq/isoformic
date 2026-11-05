# Download the FANTOM CAT Companion Tables

Fetches, for one FANTOM CAT level, the ID mapping table (transcript to
CAGE cluster), the gene information table and the CAGE cluster BED file,
which
[`add_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/add_fantom_cat_tables.md)
joins into an annotation directory.

## Usage

``` r
download_fantom_cat_tables(
  level = c("lv3_robust", "lv4_stringent", "lv2_permissive", "lv1_raw"),
  tables = c("id_mapping", "gene_info", "cage_clusters"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
)
```

## Arguments

- level:

  FANTOM CAT level: `"lv3_robust"` (default), `"lv1_raw"`,
  `"lv2_permissive"` or `"lv4_stringent"`.

- tables:

  Which tables to download.

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

A named character vector of file paths (`id_mapping`, `gene_info`,
`cage_clusters`).

## Examples

``` r
if (FALSE) { # \dontrun{
tables <- download_fantom_cat_tables("lv3_robust")
tables
} # }
```
