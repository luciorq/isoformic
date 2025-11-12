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

This function uses the `[tools::R_user_dir()]` function to determine the
user cache directory.
