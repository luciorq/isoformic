# Windows

Why the Windows environment differs from Linux and macOS, and how its R
packages are installed. Read `README.md` first. Evidence marks (`[run]`,
`[solve]`, `[docs]`) are defined there.

**Nothing in this file was executed on Windows.** pixi solved and locked
the `win-64` platform on Linux `[solve]`, and the package script was tested
on Linux with a simulation of its inputs. The first real run is the
`windows-latest` job of the `pixi` workflow.

## 1. What conda lacks on Windows

| Package | State on 2026-11-11 |
| --- | --- |
| `r-duckdb` | The conda-forge recipe skips Windows. The newest `win-64` build is 1.1.2; isoformic needs 1.5.4 |
| Every `bioconductor-*` package | bioconda builds for Linux and macOS only |
| `r-condathis` | On no channel; the inline git build would need the R tool chain on Windows. Not tried |
| LaTeX | The TinyTeX recipe has Linux and macOS sources only |

Two results follow:

1. The `[dev]` table cannot be solved for `win-64`, because it asks for
   `r-duckdb >=1.5.4` `[solve]`. So it sits under `[target.unix.dev]`, and
   Windows does not get its dependencies from `DESCRIPTION`.
2. The Bioconductor `Suggests`, condathis and TinyTeX sit under
   `[target.unix.dependencies]`.

## 2. The design

The Windows environment has two layers:

| Layer | Source | Locked |
| --- | --- | --- |
| R, the `Imports` that conda-forge has, the CRAN `Suggests`, the development packages, the command line tools | conda-forge, through `pixi.toml` | Yes |
| duckdb, condathis, the Bioconductor packages and the dependencies they still need | Binary packages from CRAN, r-universe and Bioconductor, through `.github/pixi/windows-packages.R` | No |

`[target.win-64.dependencies]` lists by hand the 12 `Imports` that
conda-forge has for Windows. This list repeats `DESCRIPTION`: a new import
must be added there too, or the script installs it from CRAN (which works,
but outside the lock file).

pixi has no way to say "take this package from a CRAN-like repository" in
the manifest. A feature or a target can change the conda channels and the
conda packages, not the kind of package. So the second layer is a task.

## 3. The script

`pixi run r-packages` runs `.github/pixi/windows-packages.R`. It:

1. Reads `Imports` and `Suggests`, with the minimum versions, from
   `DESCRIPTION`.
2. Finds the packages that the library of the environment lacks or has in
   an older version than the minimum.
3. Installs them with `install.packages(type = "win.binary")` into the
   library of the environment (`.Library`), with the dependencies they
   still need.
4. Checks that each one is installed and loads. `install.packages()` only
   warns on a failure, so the script stops with an error instead.

The repositories, in this order:

| Repository | For |
| --- | --- |
| `https://bioconductor.org/packages/3.22/{bioc,data/annotation,data/experiment}` | Bioconductor |
| `https://luciorq.r-universe.dev` | condathis |
| `https://cloud.r-project.org` | duckdb and the other dependencies |

Properties:

- The second run installs nothing: it only reads the library. The other
  Windows tasks (`test`, `check`, `document`, `examples`) depend on
  `r-packages`, so they always run it first.
- `ISOFORMIC_BIOC_REPOS` (comma separated) replaces the Bioconductor
  repositories, for example with `https://bioc-release.r-universe.dev`.
- The script stops with an error on Linux and macOS.
- In the simulation it selected 14 packages, which need 80 more
  dependencies.

## 4. Why Bioconductor 3.22

bioconda builds one Bioconductor release for each R version: 3.22 for
R 4.5. The script names the same release (`bioc_version <- "3.22"`), so the
three platforms use the same Bioconductor versions. `r-base` in
`pixi.toml`, the release on bioconda and `bioc_version` change together.

A Bioconductor release does not change after the next one is out, so the
unlocked layer is stable for Bioconductor packages once 3.23 exists. It is
not stable for CRAN and r-universe packages: they are the newest binary of
the day.

## 5. Tasks on Windows

The justfile recipes are bash scripts, so the Windows tasks call R and
`air` directly:

| Task | Command |
| --- | --- |
| `r-packages` | `Rscript .github/pixi/windows-packages.R` |
| `test` | `devtools::test(stop_on_failure=TRUE)` |
| `check` | `rcmdcheck::rcmdcheck(args='--no-manual', error_on='warning')` |
| `document` | `devtools::document()` |
| `examples` | `devtools::run_examples()` |
| `format` | `air format ./R/ ./tests/ ./data-raw/` |

No `lint`, `test-full` or article tasks, and no `pkg` environment.

## 6. Risks

- **Binary compatibility.** The CRAN and Bioconductor binaries are built
  for the R of CRAN (Rtools). The R of conda-forge is built with a
  different tool chain. Packages with compiled code can fail to load. This
  is why the script loads every package it installs. Not tested.
- **Not locked.** Two installs on different days can differ.
- **`pixi install` can undo it.** When pixi updates the environment it can
  remove or replace files that it does not know. The tasks run the script
  again, so the packages come back.
- **A repeated list.** `[target.win-64.dependencies]` repeats part of
  `Imports` without their minimum versions.

## 7. Alternatives that were set aside

| Alternative | Why not |
| --- | --- |
| Drop Windows from `platforms` | The user asked for a Windows route |
| Build duckdb and the Bioconductor packages from source with recipes | Long builds, a compiler tool chain on Windows, one recipe per package |
| WSL | Works today with the `linux-64` environment and needs nothing from this file. It is the simplest route for a Windows user |
| Wait for conda-forge `r-duckdb` on Windows | Then `[dev]` could move back to the top level, and only Bioconductor would need the script |
