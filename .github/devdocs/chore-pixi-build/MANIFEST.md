# The manifest, table by table

How `pixi.toml` is put together and why. Read `README.md` first. Evidence
marks (`[run]`, `[solve]`, `[docs]`) are defined there.

## 1. A template for another R package

This is the smallest manifest that gives a development environment for an R
package whose `Imports` are all on conda-forge `[run]`:

```toml
[workspace]
name = "mypackage"
channels = ["conda-forge"]
platforms = ["linux-64", "osx-arm64"]
preview = ["pixi-build"]
requires-pixi = ">=0.82.1"

# The R package in this directory. Name and version come from DESCRIPTION
# when they are left out; the conda name is `r-<lower case name>`.
[package]
name = "r-mypackage"

[package.build.backend]
name = "pixi-build-r"
version = "0.2.*"
channels = ["conda-forge"]

# Build against the R of the workspace (see section 4).
[package.host-dependencies]
r-base = "4.5.*"

# Install the dependencies of the package, not the package.
[dev]
r-mypackage = { path = "." }

[dependencies]
r-base = "4.5.*"
r-devtools = "*"
r-testthat = ">=3.0.0"

# Keep a user library and the R startup files of the user out.
[activation.env]
R_ENVIRON_USER = ""
R_LIBS_USER = "$CONDA_PREFIX/lib/R/library"
R_MAKEVARS_USER = ""
R_PROFILE_USER = ""

[tasks]
test = "Rscript -e \"devtools::test(stop_on_failure=TRUE)\""
```

Add to `.gitignore`:

```text
.pixi/*
!.pixi/config.toml
```

Add to `.Rbuildignore`: `^pixi\.toml$`, `^pixi\.lock$`, `^\.pixi$`,
`^\.gitattributes$`. Add to `.gitattributes`:

```text
pixi.lock merge=binary linguist-language=YAML linguist-generated=true -diff
```

Commit `pixi.lock`.

## 2. `[package]` and the `pixi-build-r` backend

The backend turns the R package into a conda package with
`R CMD INSTALL`. From `DESCRIPTION` it reads `[docs]`, and for the marked
items `[run]`:

| Field | Use |
| --- | --- |
| `Package`, `Version` | Name and version, when the manifest does not give them `[run]` |
| `Imports`, `Depends` | Host and run dependencies, as `r-<lower case name>`, with the minimum version (`(>= 1.5)` becomes `>=1.5`) `[run]` |
| `LinkingTo` | Host dependencies |
| `SystemRequirements` | Adds `rust`, `cmake` or `pkg-config` to the build when it names them |
| A `src/` directory or `LinkingTo` | Adds the C, C++ and Fortran compilers |
| `License` | Mapped to SPDX, best effort |

Base packages (`stats`, `utils`, `grid`, `parallel`, `tools`) are dropped,
because they come with `r-base`. Recommended packages (`Matrix`, `MASS`)
are ordinary conda packages.

What it does not read:

- `Suggests` and `Enhances`.
- `Config/Needs/*` fields.
- `Remotes`.
- Bioconductor: every name becomes `r-<name>`, and Bioconductor packages
  are `bioconductor-<name>` on bioconda. A Bioconductor package in `Imports`
  makes the solve fail with `No candidates were found for r-biobase *`
  `[run]`. The backend has no option to change the mapping.

The options under `[package.build.config]` are `extra-args` (for
`R CMD INSTALL`), `env`, `extra-input-globs`, `compilers` and `channels`
`[docs]`. There are two channel settings and they are not the same thing:

- `[package.build.backend] channels`: where pixi gets the backend program.
- `[package.build.config] channels`: where the backend looks for the R
  dependencies. The default is conda-forge.

## 3. Three ways to use the package of the workspace

| Form | What pixi does | Use |
| --- | --- | --- |
| `[dev]` with `r-pkg = { path = "." }` | Installs the build, host and run dependencies of the package. Does not build or install the package | Daily development with `devtools::load_all()` |
| `[dependencies]` (or a feature) with `r-pkg = { path = "." }` | Builds the package in a separate build environment and installs the result. Builds again when a source file changes | Tools that call `library(pkg)` in a new R process |
| `pixi build` | Writes a `.conda` file | Distribution |

isoformic uses the first form in the default environment and the second in
the `pkg` feature:

```toml
[target.unix.dev]
r-isoformic = { path = "." }

[feature.pkg]
platforms = ["linux-64", "osx-arm64"]

[feature.pkg.dependencies]
r-isoformic = { path = "." }

[environments]
pkg = ["pkg"]
```

Findings:

- With `[dev]`, `DESCRIPTION` is the only place that lists `Imports`. A
  change there needs `pixi lock` and nothing else `[run]`.
- With `[dev]`, a change to a source file does not start a build `[run]`.
  The first prototype had the path dependency in `[dependencies]` and built
  the package on every change.
- An article that runs code failed to render in the default environment and
  rendered in `pkg` `[run]`: Quarto starts a new R process that calls
  `library(isoformic)`.
- `[dev]` can sit under a target (`[target.unix.dev]`). isoformic needs
  this, because the dependencies of the package cannot be solved on Windows
  (`WINDOWS.md`).
- A feature can restrict the platforms, so `pkg` does not exist on Windows.

## 4. Pin R in the host dependencies

The backend solves the host environment of each package it builds without
the `r-base` of the workspace, so it takes the newest R
([prefix-dev/pixi#7214](https://github.com/prefix-dev/pixi/issues/7214),
open). A test package built against R 4.6 and then could not be installed
next to `r-base = "4.5.*"` `[run]`. isoformic did not fail without the pin,
but the pin is there as a guard:

```toml
[package.host-dependencies]
r-base = "4.5.*"
```

An inline package needs the same pin in its own definition
(`PACKAGE-SOURCES.md`, section 2). After the pin, the built packages depend
on `r-base >=4.5,<4.6.0a0` `[run]`.

## 5. What is listed by hand

| Group | Table | Reason |
| --- | --- | --- |
| `Suggests` on CRAN | `[dependencies]` | The backend does not read `Suggests` |
| `Suggests` on Bioconductor | `[target.unix.dependencies]` | bioconda names, and bioconda has no Windows builds |
| condathis | `[target.unix.dependencies]` | On no conda channel: an inline git package |
| `Config/Needs/dev`, `Config/Needs/website` | `[dependencies]` | The backend does not read them |
| Command line tools of the justfile | `[dependencies]` | `just`, `air`, `jarl`, `lychee`, `quarto`, `qpdf`, `ghostscript`, `rumdl`, `markdownlint-cli2` |
| TinyTeX | `[target.unix.dependencies]` | `LATEX.md` |
| `Imports` on Windows | `[target.win-64.dependencies]` | `[dev]` is not used there |

Two pins tie the manifest to other files: `r-roxygen2 = "8.1.*"` follows
`Config/roxygen2/version` in `DESCRIPTION`, and `r-base = "4.5.*"` follows
the Bioconductor release on bioconda (3.22 for R 4.5).

The first test run in the environment found a dependency that the global
library had hidden: one test needs `readr` and `jsonlite`, which were not
in `Suggests`. An isolated environment shows this kind of gap; a shared
library does not.

## 6. Activation

```toml
[target.unix.activation.env]
R_ENVIRON_USER = ""
R_LIBS_USER = "$CONDA_PREFIX/lib/R/library"
R_MAKEVARS_USER = ""
R_PROFILE_USER = ""
TEXMFCONFIG = "$CONDA_PREFIX/opt/tinytex/texmf-config"
TEXMFHOME = "$CONDA_PREFIX/opt/tinytex/texmf-local"
TEXMFVAR = "$CONDA_PREFIX/opt/tinytex/texmf-cache"

[target.win-64.activation.env]
R_ENVIRON_USER = ""
R_LIBS_USER = "%CONDA_PREFIX%\\Lib\\R\\library"
R_MAKEVARS_USER = ""
R_PROFILE_USER = ""
```

- `R_LIBS_USER` set to the library of the environment keeps a user library
  (`~/R/...`) out of `.libPaths()`. Result: one entry `[run]`.
- The three empty `R_*` variables keep the R startup files of the user out.
  R treats a variable that is set and empty as "no file": it does not fall
  back to `~/.Rprofile`, `~/.Renviron` or `~/.R/Makevars`, and it does not
  read a `.Rprofile` or `.Renviron` of the project directory. `[activation.env]`
  accepts `""` and exports an empty variable `[run]`. Result in a shell that
  exports the three variables: they are empty inside the environment,
  `getOption("repos")` is `@CRAN@` and `tools::makevars_user()` is empty
  `[run]`. The cost is in `CAVEATS.md`, C12: no CRAN mirror by default.
- The `TEXMF*` variables keep TeX away from trees in the home directory
  (`LATEX.md`, section 4).

`pixi run` and `pixi shell` pass the other variables of the calling shell
through.

The justfile has a `cran_mirror` variable, and every recipe that installs
packages or queries CRAN sets `options(repos = c(CRAN = ...))` from it
(`build-readme`, `build-vignettes`, `check`, `install-deps`,
`install-local`, `install-dev-deps`, `build-pkgdown-website`,
`pre-release`). A new recipe that needs CRAN must do the same.

## 7. Tasks

On Linux and macOS each task is `just <name>`, so the justfile stays the
one place with the commands, and `just` comes from the environment. On
Windows the tasks call R directly, because the justfile recipes are bash
scripts. A task under a feature (`[feature.pkg.tasks]`) runs in an
environment that has the feature: `pixi run -e pkg render-articles`.

## 8. The lock file and the pixi version

- `pixi.lock` records every package for the three platforms, the commit of
  each git package and the path of each source package.
- `pixi install --locked` fails when the lock file does not match the
  manifest or `DESCRIPTION`. CI uses it (`locked: true`).
- `pixi lock --check` tells whether the lock file is current.
- `.gitattributes` treats the file as generated: no text diff, no 3-way
  merge. After a merge conflict take one side and run `pixi lock`.

Which pixi is needed:

| Feature | First pixi version |
| --- | --- |
| `package.name` optional | 0.54.0 `[docs]` |
| `[dev]` table | 0.62.0 `[docs]` |
| Inline package definitions | 0.72.0, with fixes in 0.72.1 and 0.74.0 `[docs]` |
| `pixi-build-r` 0.2.0 | 0.82.1 `[docs]` |
| This lock file | 0.82.1 `[run]` |

pixi 0.82.0, 0.81.0, 0.80.0, 0.78.0, 0.74.0 and 0.72.1 reject the lock file
with "lock file not up-to-date with the workspace" `[run]`. So
`requires-pixi = ">=0.82.1"`.

## 9. Continuous integration

`.github/workflows/pixi.yaml`:

- `prefix-dev/setup-pixi` with a fixed `pixi-version`, `locked: true` and
  `cache: true`, on Ubuntu, macOS and Windows.
- On Windows, `pixi run r-packages`.
- The tests with `devtools::test(stop_on_failure=TRUE)`, not `just test`
  (`CAVEATS.md`, C13).
- On Linux and macOS, `R CMD Rd2pdf`, to test the LaTeX package.
- It runs when `pixi.toml`, `pixi.lock`, `DESCRIPTION`, `.github/pixi/**`
  or the workflow change, and on manual dispatch.

The other workflows of the repository (`R-CMD-check`, `pkgdown`,
`full-check`) do not use pixi.
