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
