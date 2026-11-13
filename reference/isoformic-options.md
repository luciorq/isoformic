# Options and environment variables used by isoformic

isoformic has one package option, `isoformic.frozen_mirrors`, and is
otherwise configured through base R options, the environment variables
that locate the user cache and the settings of the packages it builds
on. This page lists every setting the package reads or is affected by,
what it changes and where it is used. Settings are read at call time, so
they can be changed with
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
  with the DuckDB defaults: all CPU cores (`threads`), 80% of the RAM as
  `memory_limit` and a `.tmp` directory under the working directory when
  a query spills to disk. The package does not expose these; they are
  DuckDB session settings and the connections are not shared with user
  code.

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
#> /tmp/RtmpmIi7db/R/isoformic/frozen
```
