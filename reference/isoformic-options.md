# Options and environment variables used by isoformic

isoformic has two package options, `isoformic.frozen_mirrors` and
`isoformic.orf_cache`, and is otherwise configured through base R
options, the environment variables that locate the user cache and the
settings of the packages it builds on. This page lists every setting the
package reads or is affected by, what it changes and where it is used.
Settings are read at call time, so they can be changed with
[`options()`](https://rdrr.io/r/base/options.html) /
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html) in a session,
scoped with
[`withr::with_options()`](https://withr.r-lib.org/reference/with_options.html)
/
[`withr::with_envvar()`](https://withr.r-lib.org/reference/with_envvar.html),
or placed in `.Rprofile` / `.Renviron`.

## Package options

- `isoformic.frozen_mirrors`:

  Character vector of base URLs tried first by
  [`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  (and every function that loads a frozen dataset, such as
  [`add_appris()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)
  or `polya_sites(polya_atlas = "polyasite")`) before the `frozen-data`
  GitHub release and the Zenodo record returned by
  [`frozen_mirrors()`](https://luciorq.github.io/isoformic/reference/frozen_data.md).
  Each URL is joined with the file name of the dataset, for example
  `options(isoformic.frozen_mirrors = "https://zenodo.org/records/<id>/files/")`
  pins an analysis to one Zenodo deposit, and
  `"file:///data/isoformic-frozen/"` serves the files from a local or
  institutional mirror. The MD5 recorded in
  [`frozen_datasets()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  is verified whatever the mirror. Default: unset (no extra mirror).

- `isoformic.orf_cache`:

  `FALSE` stops
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  from reading and writing `annot_data_orfs.parquet`, the table of the
  ORFs computed so far that it keeps in the annotation directory. The
  ORFs are then computed in every session. Default: unset (`TRUE`).

## Base R options

- `timeout`:

  Timeout in seconds of
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html)
  and [`url()`](https://rdrr.io/r/base/connections.html) connections.
  [`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md),
  [`download_reference_url()`](https://luciorq.github.io/isoformic/reference/download_reference_url.md),
  [`download_hprc_annotation()`](https://luciorq.github.io/isoformic/reference/download_hprc_annotation.md),
  [`download_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/download_fantom_cat_tables.md),
  [`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
  and the builders in
  [build_frozen_data](https://luciorq.github.io/isoformic/reference/build_frozen_data.md)
  raise it, for the duration of the download only, to at least the
  `timeout_limit` of the call (3600 seconds by default, large reference
  files are several gigabytes) and never lower it. Set a larger value
  for slow connections.

- `download.file.method`, `download.file.extra`:

  Define what the `method = "auto"` default of the download functions
  resolves to and the extra command line flags of the `"wget"` and
  `"curl"` methods. Pass `method` explicitly to override them for one
  call.

- `rlang_interactive`:

  Overrides [`interactive()`](https://rdrr.io/r/base/interactive.html)
  for
  [`rlang::is_interactive()`](https://rlang.r-lib.org/reference/is_interactive.html).
  In a non-interactive session the download progress is silenced
  (`quiet = TRUE` in
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html))
  and missing optional packages (`fgsea`, `condathis`, `readxl`,
  `SummarizedExperiment`) raise an error instead of an installation
  prompt.

- `cli.num_colors`, `cli.unicode`, `cli.default_handler`, `cli.width`,
  `cli.progress_show_after`:

  Rendering of every message, warning and error, which the package emits
  through the `cli` package. `cli.unicode = FALSE` and
  `cli.num_colors = 1` produce plain ASCII output for logs; see
  `?cli::cli-config`.

## Environment variables

- `R_USER_CACHE_DIR`, `XDG_CACHE_HOME`:

  Root of the user cache.
  [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md)
  resolves `tools::R_user_dir("isoformic", "cache")`, which honours
  `R_USER_CACHE_DIR` first, then `XDG_CACHE_HOME` (`~/.cache` by default
  on Linux). On macOS the package sets `XDG_CACHE_HOME` to `~/.cache`
  for the lookup when it is unset, so the cache lives in
  `~/.cache/R/isoformic/` instead of
  `~/Library/Caches/org.R-project.R/R/isoformic/` and never contains a
  space. The cache holds downloaded references
  ([`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)),
  the frozen datasets (`frozen/`), their sources (`sources/<name>/`) and
  the Zenodo bundle (`zenodo_deposit/`). Point it to a shared or scratch
  disk on HPC systems.

- `LOCALAPPDATA`:

  On Windows, the root that
  [`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) uses
  when `R_USER_CACHE_DIR` is unset.

- `HOME`, `USERPROFILE`:

  Define the `~` that the defaults above and
  [`fs::path_home()`](https://fs.r-lib.org/reference/path_expand.html)
  expand to.

- `http_proxy`, `https_proxy`, `ftp_proxy`, `no_proxy`:

  Proxy settings honoured by
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html)
  and [`url()`](https://rdrr.io/r/base/connections.html) (`libcurl` and
  `internal` methods) and therefore by every download in the package,
  including the connectivity check
  [`curl::has_internet()`](https://jeroen.r-universe.dev/curl/reference/nslookup.html).

- `CURL_CA_BUNDLE`, `SSL_CERT_FILE`:

  Certificate bundle used by the `libcurl` download method, needed
  behind TLS-intercepting proxies.

## Settings of other packages

- DuckDB:

  Every Parquet read and write
  ([`read_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  [`write_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  [`read_tabular()`](https://luciorq.github.io/isoformic/reference/parquet.md))
  and the heavy joins run in a short-lived in-memory DuckDB connection
  with the DuckDB defaults: all CPU cores (`threads`) and 80% of the RAM
  as `memory_limit`. The package does not expose these; they are DuckDB
  session settings and the connections are not shared with user code.
  Extensions, secrets and the files of a query that spills to disk go to
  `isoformic-duckdb/` under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), unless the
  `duckdb.home` option or the `DUCKDB_R_HOME` environment variable is
  set: DuckDB then uses its own locations under that directory.

- `condathis.micromamba_path`, `CONDATHIS_MICROMAMBA_PATH`,
  `CONDA_PREFIX`, `XDG_DATA_HOME`:

  Read by the `condathis` package, which
  [`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md)
  and
  [`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
  use to run `salmon` in a micromamba environment: where the
  `micromamba` binary is found and where environments are created. See
  [`?condathis::condathis`](https://luciorq.github.io/condathis/reference/condathis-package.html).

- `BIOCONDUCTOR_*`, `R_LIBS*`:

  Only relevant for installing the optional packages (`fgsea`,
  `SummarizedExperiment`); the package itself does not read them.

## Development and test settings

These are not read by the package code but control its test suite and
checks (see `justfile`):

- `ISOFORMIC_FULL`:

  `"true"` runs the tests and articles that download references or run
  external tools (`just test-full`, `just render-articles-full`). Unset
  by default: the tests, examples and site then use only the data
  bundled with the package.

- `TESTTHAT_CPUS`:

  Number of parallel testthat workers used by `just test`.

- `NOT_CRAN`:

  Set to `"true"` by `devtools::test()`; tests that need optional
  packages use
  [`testthat::skip_if_not_installed()`](https://testthat.r-lib.org/reference/skip.html)
  rather than network access, so no test downloads data.

- `_R_CHECK_SYSTEM_CLOCK_`:

  Set to `"0"` by `just check` to skip the clock check of
  `R CMD check --as-cran`.

- `HOME`, `R_USER_CACHE_DIR`, `XDG_CACHE_HOME`, `R_USER_DATA_DIR`,
  `XDG_DATA_HOME`, `LOCALAPPDATA`, `APPDATA`:

  The test setup redirects all of them to temporary directories so that
  tests never touch the user cache.

## See also

[`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md),
[frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md),
[`download_reference()`](https://luciorq.github.io/isoformic/reference/download_reference.md)

## Examples

``` r
# Pin the frozen datasets to one Zenodo deposit for a reproducible run
withr::with_options(
  list(isoformic.frozen_mirrors = "https://zenodo.org/records/1234567/files/"),
  frozen_mirrors()
)
#> [1] "https://zenodo.org/records/1234567/files/"                          
#> [2] "https://github.com/luciorq/isoformic/releases/download/frozen-data/"

# Where the cache is, and how to move it
get_isoformic_cache()
#> /home/runner/.cache/R/isoformic
withr::with_envvar(
  list(R_USER_CACHE_DIR = tempdir()),
  get_isoformic_cache("frozen")
)
#> /tmp/RtmpatoNKq/R/isoformic/frozen
```
