# Download Any Reference File by URL

Caches a file from an arbitrary URL in the same way as
[`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)
(no re-download when present, partial files removed on failure). Useful
for annotations without a provider layout, such as the T2T consortium
Liftoff files or CHESS releases.

## Usage

``` r
download_reference_url(
  url,
  output_path = ":cache:",
  file_name = NULL,
  timeout_limit = 3600,
  method = "auto"
)
```

## Arguments

- url:

  URL of the file.

- output_path:

  Directory where the file is saved; `":cache:"` (the default) is the
  package cache.

- file_name:

  File name to use; defaults to the last part of the URL.

- timeout_limit:

  Maximum time in seconds for the download. Takes precedence over
  `options("timeout")`. Defaults to `3600`.

- method:

  Method used by
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).
  Defaults to `"auto"`.

## Value

The full path to the file.

## Examples

``` r
if (FALSE) { # \dontrun{
liftoff <- download_reference_url(paste0(
  "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/",
  "assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
))
annot <- prepare_isoformic_annotation(liftoff)
annot
} # }
```
