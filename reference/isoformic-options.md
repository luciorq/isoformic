# Options and environment variables used by isoformic

isoformic has five package options, `isoformic.duckdb_memory_limit`,
`isoformic.duckdb_threads`, `isoformic.frozen_mirrors`,
`isoformic.orf_cache` and `isoformic.workers`, and is otherwise
configured through base R options, the environment variables that locate
the user cache and the settings of the packages it builds on. This page
lists every setting the package reads or is affected by, what it changes
and where it is used. Settings are read at call time, so they can be
changed with [`options()`](https://rdrr.io/r/base/options.html) /
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html) in a session,
scoped with
[`withr::with_options()`](https://withr.r-lib.org/reference/with_options.html)
/
[`withr::with_envvar()`](https://withr.r-lib.org/reference/with_envvar.html),
or placed in `.Rprofile` / `.Renviron`.

## Package options

- `isoformic.duckdb_threads`:

  Number of threads of the DuckDB connections that read and write
  Parquet and run the annotation and uncertainty queries: a single
  positive whole number, such as `4`. Default: unset. When
  `SLURM_CPUS_PER_TASK` is set or the affinity mask of the R process
  excludes some CPUs (`taskset`, schedulers that pin cores), isoformic
  uses the smallest of `SLURM_CPUS_PER_TASK`, the CPUs of the mask and
  the cgroup CPU quota. Otherwise DuckDB uses all the CPUs of the
  machine, or the cgroup CPU quota when there is one.

- `isoformic.duckdb_memory_limit`:

  Memory limit of each DuckDB connection: a single string with a size
  and a unit, such as `"8GB"` or `"512MiB"` (units KB, MB, GB, TB, KiB,
  MiB, GiB, TiB), of at least 1 MiB. A query above it spills to disk
  when it can and fails otherwise. Default: unset. isoformic then sets
  80% of `SLURM_MEM_PER_NODE`, or of `SLURM_MEM_PER_CPU` times
  `SLURM_CPUS_PER_TASK` (or 1), when one is set, or 80% of the cgroup
  memory limit when that is smaller, and at least 64 MiB; otherwise
  DuckDB uses 80% of the RAM, or of the cgroup memory limit when there
  is one. Set it explicitly under `ulimit -v`, which DuckDB does not
  see. Both DuckDB options are read when a connection opens a database.
  A database file that is already open in the session (by a connection
  that is still open) is shared and keeps the settings it was opened
  with.

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

- `isoformic.workers`:

  Default of the `workers` argument of
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md),
  [`run_productivity_enrichment()`](https://luciorq.github.io/isoformic/reference/productivity.md)
  and
  [`run_apa_enrichment()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md):
  the number of `mirai` daemons that run their `fgsea` calls in
  parallel, a whole number of at least 1, or `"auto"` for
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html).
  The count is capped at the number of `fgsea` calls. Default: unset (1,
  serial, the results of earlier versions). Above 1 every `fgsea` call
  gets its own seed, so the sampled p-values differ slightly from a
  serial run but not between worker counts; see *Workers* in
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md).

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
  and missing optional packages (for example `fgsea`, `condathis`,
  `readxl`, `Matrix`, `SummarizedExperiment`, or `edgeR`, `limma`,
  `satuRn`, `DEXSeq` and `BiocParallel` for
  [`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md),
  `mirai` for `workers` above 1 in the enrichment functions) raise an
  error instead of an installation prompt.

- `cli.num_colors`, `cli.unicode`, `cli.default_handler`, `cli.width`,
  `cli.progress_show_after`, `crayon.enabled`, `crayon.colors`:

  Rendering of every message, warning and error, which the package emits
  through the `cli` package (`cli` also honours the older `crayon.*`
  options). `cli.unicode = FALSE` and `cli.num_colors = 1` produce plain
  ASCII output for logs; `crayon.enabled = FALSE` turns off ANSI
  colours, and `just test` and `just test-full` set it together with
  `cli.unicode = FALSE`; see `?cli::cli-config`.

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
  [`curl::has_internet()`](https://rdrr.io/pkg/curl/man/nslookup.html)
  that the download functions run first when `curl` is installed.

- `CURL_CA_BUNDLE`, `SSL_CERT_FILE`:

  Certificate bundle used by the `libcurl` download method, needed
  behind TLS-intercepting proxies.

- `SLURM_CPUS_PER_TASK`, `SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU`:

  Set by Slurm inside a job. When the `isoformic.duckdb_threads` and
  `isoformic.duckdb_memory_limit` options are unset, they bound the
  threads and set the memory limit of the DuckDB connections (memory in
  megabytes, as Slurm writes it; `0` is ignored), so that DuckDB sizes
  itself to the job and not to the node.

- `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
  `VECLIB_MAXIMUM_THREADS`, `R_DATATABLE_NUM_THREADS`, `R_LIBS`:

  Set, for the daemons only, when the enrichment functions start `mirai`
  daemons (`workers` above 1, or
  [`isoformic_workers()`](https://luciorq.github.io/isoformic/reference/isoformic_workers.md)):
  the thread variables to 1, so that each daemon uses one core, and
  `R_LIBS` to the session's
  [`.libPaths()`](https://rdrr.io/r/base/libPaths.html), so that the
  daemons find the packages of the session. The session's own values are
  not changed.

## Settings of other packages

- DuckDB:

  Every Parquet read and write
  ([`read_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  [`write_parquet()`](https://luciorq.github.io/isoformic/reference/parquet.md),
  [`read_tabular()`](https://luciorq.github.io/isoformic/reference/parquet.md))
  and the heavy joins run in a short-lived in-memory DuckDB connection,
  which is not shared with user code. Its `threads` and `memory_limit`
  come from the `isoformic.duckdb_threads` and
  `isoformic.duckdb_memory_limit` options or the Slurm and affinity
  defaults described above, and otherwise from DuckDB. Besides these and
  the storage locations below, the package only turns off the DuckDB
  progress bar (`enable_progress_bar`). Extensions, secrets and the
  files of a query that spills to disk go to `isoformic-duckdb/` under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), unless the
  `duckdb.home` option or the `DUCKDB_R_HOME` environment variable is
  set: DuckDB then uses its own locations under that directory.

- BiocParallel:

  Each `fgsea` call of the enrichment functions runs with
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  unless `BPPARAM` or `nproc` is passed through `...`, so a backend set
  with
  [`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  is not used; the calls run in parallel through `workers` instead. The
  `workers` argument of
  [`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
  is different: it builds a
  [`BiocParallel::MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  on Unix-like systems.

- mirai:

  With `workers` above 1 the enrichment functions run on the daemons of
  the `mirai` compute profile `"isoformic"`. When that profile has no
  daemons, they start them (with one thread each) and keep them for the
  session, resizing them to the `workers` of later calls and replacing
  those that died;
  [`isoformic_workers()`](https://luciorq.github.io/isoformic/reference/isoformic_workers.md)
  starts a fixed number ahead or stops them, and they are stopped when
  the session ends or isoformic is unloaded. Daemons that the user
  started on the profile (`mirai::daemons(n, .compute = "isoformic")`)
  are used as they are and left running; each call puts the session's
  library paths first on the daemon while it runs, but they keep their
  own thread settings. Starting daemons with
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  draws from the session's random number stream (isoformic's own start
  does not). Needs `mirai` 2.3.0 or later. Daemons of other profiles,
  including the default one, are not used or changed.

- parallelly:

  `workers = "auto"` (or the `isoformic.workers` option) calls
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html),
  which reads the `mc.cores` option, the `MC_CORES`,
  `SLURM_CPUS_PER_TASK` and other scheduler environment variables, the
  CPU affinity mask and the cgroup CPU quota, and its own
  `parallelly.availableCores.*` options and `R_PARALLELLY_*` variables
  (see
  [`?parallelly::availableCores`](https://parallelly.futureverse.org/reference/availableCores.html)).
  Without parallelly, `"auto"` is one worker.

- `condathis.micromamba_path`, `CONDATHIS_MICROMAMBA_PATH`,
  `CONDA_PREFIX`, `R_USER_DATA_DIR`, `XDG_DATA_HOME`:

  Read by the `condathis` package, which
  [`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md),
  [`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
  and
  [`run_terminus()`](https://luciorq.github.io/isoformic/reference/run_terminus.md)
  use to run `salmon` and Terminus in micromamba environments: where the
  `micromamba` binary is found and where environments are created. See
  [`?condathis::condathis`](https://luciorq.github.io/condathis/reference/condathis-package.html).

- `BIOCONDUCTOR_*`, `R_LIBS*`:

  Only relevant for installing the optional packages (`fgsea`,
  `SummarizedExperiment`, `edgeR`, ...); the package itself does not
  read them.

## Development and test settings

These are not read by the package code but control its test suite and
checks (see `justfile`):

- `ISOFORMIC_FULL`:

  `"true"` runs the tests and articles that download references or run
  external tools (`just test-full`, `just render-articles-full`,
  `just benchmark-memory`). Unset by default: the tests, examples and
  site then use only the data bundled with the package.

- `TESTTHAT_CPUS`:

  Number of parallel testthat workers used by `just test`; it only
  matters when parallel mode is on.

- `TESTTHAT_PARALLEL`:

  `"true"` or `"false"`; overrides `Config/testthat/parallel: true` in
  `DESCRIPTION`. `TESTTHAT_PARALLEL=false` runs the test files serially
  in one process.

- `NOT_CRAN`:

  Set to `"true"` by `devtools::test()`; tests that need optional
  packages use
  [`testthat::skip_if_not_installed()`](https://testthat.r-lib.org/reference/skip.html),
  and only the tests gated by `ISOFORMIC_FULL` download data.

- `_R_CHECK_SYSTEM_CLOCK_`:

  Set to `"0"` by `just check` to skip the clock check of
  `R CMD check --as-cran`.

- `HOME`, `USERPROFILE`, `R_USER_CACHE_DIR`, `XDG_CACHE_HOME`,
  `R_USER_DATA_DIR`, `XDG_DATA_HOME`, `LOCALAPPDATA`, `APPDATA`:

  The test setup redirects all of them to temporary directories so that
  tests never touch the user cache.

- `warnPartialMatchArgs`, `warnPartialMatchAttr`,
  `warnPartialMatchDollar`:

  Set to `TRUE` by the test setup, so that a partial match of an
  argument, attribute or `$` name is a warning.

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
#> /tmp/RtmpPCTCTT/R/isoformic/frozen
```
