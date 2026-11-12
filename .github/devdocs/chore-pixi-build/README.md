# pixi for R package development: notes

Written 2026-11-12. Tested with pixi 0.84.2, `pixi-build-r` 0.2.x,
`pixi-build-rattler-build` 0.4.x and R 4.5.3 on Linux x86-64. Purpose:
record how the pixi development environment of isoformic is built, why it
has this shape, and what pixi can and cannot do for an R package today. The
notes are written so that the same setup can be copied to another R package.

**Where to start.** Section 1 is the conclusion. Sections 2 to 4 are the
daily use. Section 5 has the procedures for the changes that come up
(a new import, a new R release, a new TinyTeX). The other files hold the
detail:

| File | Content |
| --- | --- |
| [`MANIFEST.md`](MANIFEST.md) | How `pixi.toml` is put together, table by table, and a template for another package |
| [`PACKAGE-SOURCES.md`](PACKAGE-SOURCES.md) | Where an R package can come from: conda channels, git, CRAN, r-universe, Bioconductor, a recipe |
| [`WINDOWS.md`](WINDOWS.md) | Why Windows differs and how its packages are installed |
| [`LATEX.md`](LATEX.md) | The TinyTeX conda package for the PDF manual |
| [`CAVEATS.md`](CAVEATS.md) | All caveats in one list, the upstream issues and the backlog |

Evidence level is marked where it matters: `[run]` the thing was executed
on Linux x86-64, `[solve]` pixi solved and locked it for that platform but
nothing was executed there, `[docs]` the statement comes from the pixi
documentation and was not tested.

## 1. Main conclusion

pixi gives an R package a locked, isolated development environment with one
command (`pixi install`), and the `pixi-build-r` backend lets `DESCRIPTION`
stay the only place where `Imports` and their minimum versions are written.
The pieces that make this work:

1. **`[package]` with the `pixi-build-r` backend** describes the R package.
   The backend reads `DESCRIPTION`.
2. **The `[dev]` table** installs the dependencies of that package without
   building or installing the package. The loop is `devtools::load_all()`.
3. **A second environment (`pkg`)** adds the package itself, built as a
   conda package, for the tools that need `library(isoformic)` in a fresh R
   process (Quarto articles).
4. **Everything the backend does not read is listed by hand**: `Suggests`,
   Bioconductor packages, development packages and command line tools.

The limits that shaped the manifest:

- `pixi-build` is a preview feature. The manifest can need changes when pixi
  is upgraded.
- The backend ignores `Suggests` and cannot name Bioconductor packages.
- conda-forge decides the versions. A minimum version in `DESCRIPTION` that
  conda-forge does not have makes the environment unsolvable (`duckdb` 2.x).
- Windows has no bioconda and no current `r-duckdb`, so on Windows a script
  installs part of the packages outside the lock file.
- A source tarball URL (CRAN, r-universe) cannot be a dependency: pixi
  panics. A git repository or a recipe can.
- conda-forge has no LaTeX distribution, so a recipe in this repository
  repackages TinyTeX.

## 2. Quick start

```bash
pixi install            # default environment, from pixi.lock
pixi install --all      # also the `pkg` environment
pixi run test           # a task
pixi shell              # a shell with R and the tools on PATH
pixi run -e pkg render-articles
```

On Windows run `pixi run r-packages` once after `pixi install` (the other
Windows tasks do it for you).

Inside the environment isoformic is not installed. Use
`devtools::load_all()`, which is what the justfile recipes do.

The first `pixi install` builds three packages from source: condathis
(about 15 s), TinyTeX (about 1 min, a download of 150 to 210 MB) and, for
`pkg`, isoformic. With a warm cache `pixi install --all` takes about 1 min
`[run]`.

## 3. Files

| File | Role |
| --- | --- |
| `pixi.toml` | The manifest: package, environments, dependencies, tasks |
| `pixi.lock` | Exact packages for the three platforms. Committed |
| `.gitattributes` | Marks `pixi.lock` as generated and binary for merges |
| `.github/pixi/tinytex/recipe.yaml` | Recipe of the TinyTeX conda package |
| `.github/pixi/windows-packages.R` | Installs on Windows what conda lacks |
| `.github/workflows/pixi.yaml` | Installs the lock and runs the tests on Linux, macOS and Windows |
| `.gitignore`, `.Rbuildignore` | Keep `.pixi/` out of git, and the pixi files out of the R package tarball |

## 4. Environments and tasks

| Environment | Platforms | Content |
| --- | --- | --- |
| `default` | linux-64, osx-arm64, win-64 | R 4.5, `Imports`, `Suggests`, development packages, command line tools, LaTeX (not on Windows). isoformic is not installed |
| `pkg` | linux-64, osx-arm64 | `default` plus isoformic built and installed |

| Task | Linux and macOS | Windows |
| --- | --- | --- |
| `test`, `document`, `examples`, `format`, `check` | `just <task>` | R or `air` called directly |
| `lint`, `test-full` | `just <task>` | none |
| `render-articles`, `render-articles-full` | `just <task>`, in `pkg` | none |
| `r-packages` | none | `.github/pixi/windows-packages.R` |

`pixi run lint` and `pixi run check` rewrite source files (`check` depends
on `lint`, which runs `air format` and the fixers). Do not run them while
another editor or session is saving the same files.

## 5. Procedures

**A new `Imports` entry.** Edit `DESCRIPTION`, then `pixi lock`. On Linux
and macOS nothing else is needed. For Windows also add `r-<name>` under
`[target.win-64.dependencies]` if conda-forge has it; if not, the Windows
script installs it. A Bioconductor package in `Imports` does not work with
the backend at all: see `PACKAGE-SOURCES.md`, section 5.

**A new `Suggests` entry.** Edit `DESCRIPTION` and add the package to
`pixi.toml`: `r-<name>` in `[dependencies]` for a CRAN package,
`bioconductor-<name>` (lower case) in `[target.unix.dependencies]` for a
Bioconductor package. Then `pixi lock`. The Windows script reads
`DESCRIPTION` and needs no change.

**A development package or a command line tool.** Add it to
`[dependencies]`, then `pixi lock`.

**A package that is on no conda channel.** Use an inline git package or a
recipe: `PACKAGE-SOURCES.md`.

**A new commit of condathis.** Read `RemoteSha` from
<https://luciorq.r-universe.dev/api/packages/condathis>, put it in `rev` of
`r-condathis`, then `pixi lock`.

**Newer versions of everything.** `pixi update`, then run the tests, then
commit `pixi.lock`.

**A new R release.** These change together:

1. `r-base` in `[dependencies]`.
2. `r-base` in `[package.host-dependencies]` and in the `host-dependencies`
   of every inline `pixi-build-r` package (condathis).
3. The Bioconductor release: bioconda builds one release per R version, and
   `bioc_version` in `.github/pixi/windows-packages.R` must name the same
   one.

**A new TinyTeX.** See the header of `.github/pixi/tinytex/recipe.yaml` and
`LATEX.md`, section 5.

**A new pixi.** Run `pixi install --locked --all` with it. If the lock file
is rejected, run `pixi lock` and state the new minimum in `requires-pixi`.
Update `pixi-version` in `.github/workflows/pixi.yaml`.

**Free disk space.** `pixi clean --build` removes the build directories
(`.pixi/bld`, several GB by `du`) and the built packages; the next
`pixi install` builds them again. `pixi clean` removes the environments.

## 6. What was tested

| Check | Result |
| --- | --- |
| `pixi install --all` on Linux | `[run]` both environments install |
| Test suite in `default` | `[run]` passes; isoformic is not in the library and `.libPaths()` has one entry |
| `R CMD check --as-cran` with the PDF manual | `[run]` 0 errors, 0 warnings, 1 note (new submission); `pdflatex`, `makeindex`, `qpdf` and `gs` come from the environment |
| An article that runs code | `[run]` fails in `default`, renders in `pkg` |
| Lock installs with older pixi | `[run]` 0.82.1, 0.83.0 and 0.84.1 install it unchanged; 0.82.0 and older reject it |
| macOS (osx-arm64) | `[solve]` only. The TinyTeX package was cross-built on Linux, never run |
| Windows (win-64) | `[solve]` only. The package script was tested with a simulation, never on Windows |

macOS and Windows get their first real run when the `pixi` workflow runs
after a push.
