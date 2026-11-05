# Read and Write Parquet Files

Thin, dependency-light wrappers around DuckDB for reading tabular files
into tibbles and for writing data frames as Parquet. These are the
readers used throughout `isoformic`: annotation tables, experiment
caches and the bundled example data are all stored as Parquet.

## Usage

``` r
read_parquet(path, columns = NULL, n_max = Inf, ...)

read_tabular(path, columns = NULL, n_max = Inf, ...)

write_parquet(
  x,
  path,
  compression = c("zstd", "snappy", "gzip", "uncompressed"),
  ...
)
```

## Arguments

- path:

  Path to the file to read or write.

- columns:

  Optional character vector of columns to read. `NULL` reads all
  columns.

- n_max:

  Maximum number of rows to read. `Inf` (default) reads all.

- ...:

  Reserved for future use; must be empty.

- x:

  A data frame (or tibble) to write.

- compression:

  Parquet compression codec.

## Value

`read_parquet()` and `read_tabular()` return a tibble. `write_parquet()`
returns `path`, invisibly.

## Details

`read_parquet()` reads a Parquet file. `read_tabular()` reads Parquet,
CSV or TSV files (`.csv`, `.tsv`, `.tab`, `.txt`, Salmon `.sf`;
optionally gzip-compressed), auto-detecting delimiter and types.
`write_parquet()` writes any data frame to a Parquet file, creating
parent directories as needed.

Factor columns are stored as Parquet strings and come back as character.
Row names are not stored; keep identifiers in a regular column.

## Examples

``` r
tmp <- tempfile(fileext = ".parquet")
write_parquet(data.frame(id = c("a", "b"), value = c(1.5, 2)), tmp)
read_parquet(tmp)
#> # A tibble: 2 × 2
#>   id    value
#>   <chr> <dbl>
#> 1 a       1.5
#> 2 b       2  
read_parquet(tmp, columns = "id", n_max = 1)
#> # A tibble: 1 × 1
#>   id   
#>   <chr>
#> 1 a    
read_tabular(tmp)
#> # A tibble: 2 × 2
#>   id    value
#>   <chr> <dbl>
#> 1 a       1.5
#> 2 b       2  
```
