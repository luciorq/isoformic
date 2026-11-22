# Retrieve System-Dependent Cache Path for Isoformic

Determines the appropriate user cache directory for the `isoformic`
package based on the operating system. On macOS, it avoids using paths
with spaces and follows the XDG base directory specification.

## Usage

``` r
get_isoformic_cache(..., ext = NULL)
```

## Arguments

- ...:

  Additional path components to append to the cache directory.

- ext:

  An optional file extension (e.g., "rds", "csv") to append to the final
  path.

## Value

A path character string representing the path to the user cache
directory for the `isoformic` package.

## Details

This function uses
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) to
determine the user cache directory, so `R_USER_CACHE_DIR` and
`XDG_CACHE_HOME` move it (see
[isoformic_options](https://luciorq.github.io/isoformic/reference/isoformic-options.md)).
The directory is not created.

## Examples

``` r
get_isoformic_cache()
#> /home/runner/.cache/R/isoformic
get_isoformic_cache("frozen")
#> /home/runner/.cache/R/isoformic/frozen
```
