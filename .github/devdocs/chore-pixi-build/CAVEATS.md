# Caveats, upstream issues and backlog

Every caveat of the pixi setup in one list. Read `README.md` first.
Evidence marks (`[run]`, `[solve]`, `[docs]`) are defined there. State on
2026-11-12, pixi 0.84.2.

## 1. Caveats

### The backend and the manifest

**C1. `pixi-build` is a preview feature.** `preview = ["pixi-build"]` is
required. The manifest format and the lock format of source packages can
change between pixi versions. Pin the pixi version in CI and test a new
pixi before using it (`README.md`, section 5).

**C2. The backend does not read `Suggests`.** Nor `Enhances`,
`Config/Needs/*` or `Remotes`. A new `Suggests` entry must be added to
`pixi.toml` by hand. Nothing detects a missing one except a failing or
skipped test.

**C3. The backend cannot name Bioconductor packages.** Every import becomes
`r-<name>`. A Bioconductor package in `Imports` or `Depends` makes the
solve fail `[run]`. There is no mapping option. isoformic is not affected
(Bioconductor in `Suggests` only).

**C4. The backend builds against the newest R.** Without `r-base` in the
`host-dependencies` of each package it builds
([prefix-dev/pixi#7214](https://github.com/prefix-dev/pixi/issues/7214)).
The pins exist in three places: `[dependencies]`,
`[package.host-dependencies]` and the inline condathis definition.

**C5. A tarball URL cannot be a source.** `url` and `sha256` make pixi
panic ([prefix-dev/pixi#7230](https://github.com/prefix-dev/pixi/issues/7230))
`[run]`. Use a git source or a recipe (`PACKAGE-SOURCES.md`).

**C6. conda-forge decides the versions.** A package can be behind CRAN by
weeks. A minimum in `DESCRIPTION` above the conda-forge version makes the
environment unsolvable. This is why `duckdb (>= 1.5.4)` stays until
conda-forge has `r-duckdb` 2.x.

**C7. Lists that repeat `DESCRIPTION`.** The `Suggests` blocks, the
`r-roxygen2` pin and the Windows imports repeat information from
`DESCRIPTION`. Nothing checks that they agree.

### The environments

**C8. isoformic is not installed in the default environment.** A tool that
calls `library(isoformic)` in a new R process fails there. Use the `pkg`
environment (`pixi run -e pkg ...`).

**C9. The `pkg` environment builds isoformic again after a source change.**
The build takes time and keeps files in `.pixi/bld/r-isoformic`.

**C10. Disk use.** `du` reports about 12 GB for `.pixi`, of which the two
environments are 3 GB each and `.pixi/bld` about 6 GB (TinyTeX 1.3 GB,
condathis 2.3 GB, isoformic 2.6 GB). Many files are hard links to the pixi
cache, so the real use is lower. `pixi clean --build` frees the build
directories.

**C11. Only two Unix platforms are locked.** `linux-64` and `osx-arm64`.
`linux-aarch64` and `osx-64` need a `platforms` entry, a TinyTeX source
and a new lock.

**C12. R has no CRAN mirror in the environment.** The activation sets
`R_PROFILE_USER`, `R_ENVIRON_USER` and `R_MAKEVARS_USER` to empty values,
so R reads no startup file of the user (`MANIFEST.md`, section 6). A mirror
that a user `.Rprofile` sets is gone with it: `getOption("repos")` is
`@CRAN@`, and `install.packages()` or `available.packages()` stop with
"trying to use CRAN without setting a mirror" `[run]`. The justfile recipes
that need a mirror set `cran_mirror` themselves. Code typed in `pixi shell`
must set `options(repos = ...)`. Two limits:

- A `.Rprofile` or `.Renviron` in the project directory is not read either.
  The repository has none.
- On Windows an empty value may unset the variable. R then reads the
  startup files of the user as it does without pixi. Not tested.

### The tasks

**C13. `just test` does not fail on a failed test.** `devtools::test()`
returns normally, and the recipe prints its own success line. CI calls
`devtools::test(stop_on_failure=TRUE)` in place of the task.

**C14. `lint` and `check` rewrite files.** `just lint` formats and applies
fixes, and `just check` depends on it. Do not run them while another
session edits the tree. For a check that does not rewrite:

```bash
pixi run Rscript -e 'rcmdcheck::rcmdcheck(args = "--as-cran")'
```

**C15. The tasks need bash.** The justfile recipes are bash scripts. The
Windows tasks are separate and fewer (`WINDOWS.md`, section 5).

### Windows and LaTeX

**C16. Windows was never run.** Solved and locked only. See `WINDOWS.md`,
section 6, for the risks (binary compatibility first).

**C17. Part of the Windows environment is not locked.** duckdb, condathis
and the Bioconductor packages come from `install.packages()`.

**C18. macOS LaTeX was never run.** The package was built for `osx-arm64`
on Linux.

**C19. `tlmgr` is not available.** A missing LaTeX package needs a change
to the recipe (`LATEX.md`, section 5).

**C20. The first install downloads TinyTeX from GitHub.** 150 to 210 MB. A
network failure there fails `pixi install`.

## 2. Upstream issues

| Issue | State | Effect here |
| --- | --- | --- |
| [prefix-dev/pixi#7230](https://github.com/prefix-dev/pixi/issues/7230) | Open | No `url` sources (C5) |
| [prefix-dev/pixi#7214](https://github.com/prefix-dev/pixi/issues/7214) | Open | Host pins of `r-base` (C4) |
| Bioconductor names in `pixi-build-r` | Not reported | C3 |
| `Suggests` in `pixi-build-r` (for example as an optional group) | Not reported | C2 |
| conda-forge `r-duckdb` for Windows | Recipe skips Windows | `WINDOWS.md` |

## 3. Backlog

Ordered by value.

1. **Read the first `pixi` workflow run.** It is the first execution on
   macOS and Windows. Expect work on Windows.
2. **A check that `pixi.toml` covers `Suggests`.** A short R script that
   compares `Suggests` in `DESCRIPTION` with the `r-*` and `bioconductor-*`
   names in `pixi.toml` and fails on a gap (C2, C7). It can run in the
   `pixi` workflow.
3. **Move to duckdb 2.x** when conda-forge has it: raise the minimum in
   `DESCRIPTION`, then `pixi lock`.
4. **Report the Bioconductor mapping gap** to the pixi build backends
   repository, with the `r-biobase` example.
5. **Use a `url` source for condathis** when #7230 is fixed, with the
   r-universe tarball and its `sha256`. This removes the build from a git
   clone.
6. **Remove the host pins** when #7214 is fixed.
7. **A `tex` feature** if the TinyTeX download becomes a problem for users
   who never build the manual (C20).
8. **Decide on the other workflows.** `R-CMD-check` and `pkgdown` use
   `r-lib/actions`, not pixi. Two ways to install the dependencies exist in
   CI. Keeping both is a test of both; it is also two things to maintain.
9. **Check the R startup variables on Windows** (C12) in the first `pixi`
   workflow run: print `Sys.getenv("R_PROFILE_USER", unset = NA)`.
