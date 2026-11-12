# isoformic development rules

- In source code, meaning `*.R` files under `./R/` and `./tests/`, functions
  from dependency packages must always be called with the explicit
  `<PKG>::<FUNCTION>()` form (for example `dplyr::filter()`,
  `testthat::expect_equal()`). Only base R functions may be called unqualified.
  The `NAMESPACE` only imports the `rlang` pronouns `.data`, `.env`, `:=` and
  `%||%`.
- Classes and methods use S7. Generics are declared with explicit argument
  signatures so that the documented usage matches the code; methods are
  documented in the generic's topic, not with their own `\usage` entry.
- All on-disk data is Parquet (annotation directories, experiment
  directories, bundled example data) read and written through DuckDB
  (`read_parquet()`, `write_parquet()`, `read_tabular()`).
- Format with `air format` (see `air.toml`); lint with `lintr::lint_package()`
  (see `.lintr`).
- Anything that downloads references or runs external tools is gated by
  `ISOFORMIC_FULL=true`: tests call `skip_if_not_full()` (they live in
  `tests/testthat/test-full.R`), and articles that need GENCODE set
  `eval = full_run` in their setup chunk. The default tests, examples and
  site build use only the bundled data. `just test-full` and
  `just render-articles-full` run the full mode locally; in CI the
  `full-check` workflow and the full `pkgdown` builds (release, weekly,
  manual dispatch; the only ones that deploy) set the flag.
- Regenerate bundled data with `Rscript data-raw/prepare_extdata.R`
  (`just extdata`); the coding potential model and the frozen evidence
  datasets have their own scripts and update procedure in
  `data-raw/README.md` (`just coding-model`, `just frozen-data`).
- Case study articles (`isoformic-case-*.qmd`) run on every site build from
  the processed bundles in `vignettes/articles/data/`; rebuild them with
  `data-raw/case_studies/` (see the header of `build_case_studies.R`). Base
  their text on the rebuilt results and check every number they state.
- `pixi.toml` defines an isolated development environment (`pixi install`,
  then `pixi run <task>` or `pixi shell`) for Linux, macOS and Windows;
  `pixi.lock` is committed. On Linux and macOS the tasks delegate to the
  justfile, and the default environment gets the `Imports` from `DESCRIPTION`
  through the `[target.unix.dev]` table and the `pixi-build-r` backend,
  without installing isoformic, so use `devtools::load_all()`. The `pkg`
  environment (Linux and macOS only) also builds and installs isoformic; the
  articles need it (`pixi run render-articles`). The backend does not read
  `Suggests` and cannot map Bioconductor names: a new `Suggests` entry, dev
  tool or command line tool must be added to `pixi.toml` by hand. conda-forge
  decides the versions there, so `duckdb` is 1.5.x until conda-forge ships
  2.x. The backend builds against the newest R unless `r-base` is pinned in
  the `host-dependencies` of each package it builds (isoformic, condathis):
  change those pins together with `r-base` in `[dependencies]`. A dependency
  cannot be a `url` to a source tarball (pixi panics, prefix-dev/pixi#7230):
  use a `git` source (`github.com/cran/<pkg>` mirrors CRAN) or a
  rattler-build recipe. The activation sets `R_PROFILE_USER`,
  `R_ENVIRON_USER` and `R_MAKEVARS_USER` to empty values, so R reads no
  startup file of the user and has no CRAN mirror: a justfile recipe that
  installs packages or queries CRAN must set
  `options(repos = c(CRAN = "{{ cran_mirror }}"))`. The notes in
  `chore-pixi-build/`, next to this file, have the findings, the procedures
  and the caveats: read them before changing `pixi.toml`.
- LaTeX for the PDF manual of `R CMD check` comes from
  `.github/pixi/tinytex/recipe.yaml`, which repackages a pinned TinyTeX
  release as a conda package for Linux and macOS (its header says how to
  update it). pixi builds it on the first `pixi install` (a download of 150
  to 210 MB) and keeps about 1.3 GB of build files in `.pixi/bld`. `tlmgr` is
  not on `PATH`: a missing LaTeX package means a newer release or a change to
  the recipe. The activation sets `TEXMFVAR`, `TEXMFHOME` and `TEXMFCONFIG`
  to directories of the environment. Windows has no LaTeX, so its `check`
  task uses `--no-manual`.
- Windows (`win-64`) has no bioconda and no current conda-forge build of
  `duckdb`. `pixi.toml` lists by hand the `Imports` that conda-forge has (add
  a new one under `[target.win-64.dependencies]`), and `pixi run r-packages`
  (`.github/pixi/windows-packages.R`, which the other Windows tasks depend
  on) installs the rest of `Imports` and `Suggests` as binary packages from
  CRAN, r-universe and Bioconductor, outside `pixi.lock`. The Windows tasks
  call R directly, because the justfile needs bash. `r-base` in `pixi.toml`,
  the Bioconductor release on bioconda and `bioc_version` in that script
  change together. The `pixi` workflow runs the tests on the three platforms.
- `R/isoformic-options.R` (help page `?"isoformic-options"`) lists every
  option and environment variable the package, its tests and its articles
  read. A change that adds, removes or changes the meaning of a
  `getOption("isoformic.*")`, a `Sys.getenv()` or a setting of a dependency
  must update that page in the same change, including the count of package
  options in its description. Before a commit, compare the page with
  `grep -rnE 'getOption\(|Sys\.getenv\(' R tests vignettes justfile`.
- Research notes with backlogs live next to this file: read
  `UNCERTAINTY-LITERATURE.md` before changing `R/uncertainty.R`,
  `R/uncertainty-io.R` or the uncertainty articles, and
  `ALPHAGENOME-BENCHMARK.md` before work on productivity evidence, splice-site
  strength, the frozen APPRIS data or the GENCODE release.
