# LaTeX for the PDF manual

How the pixi environment gets the LaTeX that `R CMD check` needs. Read
`README.md` first. Evidence marks (`[run]`, `[solve]`, `[docs]`) are
defined there.

## 1. The problem

`R CMD check` builds the PDF manual with `pdflatex` and `makeindex`, and
`--as-cran` also uses `qpdf` and Ghostscript. Without LaTeX the check must
run with `--no-manual`. An isolated environment that takes `pdflatex` from
the home directory is not isolated.

## 2. What was tried

| Option | Result |
| --- | --- |
| conda-forge `texlive-core` | The programs alone: no LaTeX packages (no `article.cls`), and its `tlmgr` does not work `[run]`. Not usable |
| `tinytex::install_tinytex()` in a task | Installs outside the environment or outside the lock file. Not locked, not reproducible |
| TinyTeX bundle `TinyTeX-1` | Has no `makeindex`, which the manual needs `[run]` |
| TinyTeX bundle `TinyTeX` (the full one), release v2026.11 | Builds the manual `[run]`. Chosen |
| `qpdf`, `ghostscript` | On conda-forge. Plain dependencies |

## 3. The recipe

`.github/pixi/tinytex/recipe.yaml` is a rattler-build recipe. pixi uses it
through a path dependency:

<!-- markdownlint-disable MD013 -->

```toml
tinytex = { path = ".github/pixi/tinytex", package.build.backend = { name = "pixi-build-rattler-build", version = "0.4.*" } }
```

<!-- markdownlint-enable MD013 -->

It downloads a pinned release bundle from
<https://github.com/rstudio/tinytex-releases>, checks the `sha256`, and
repackages the files as a conda package. Nothing is compiled. Details that
needed work:

- **`target_directory: tinytex`.** The work directory of the build also
  holds files of the build tool. Without a target directory the TeX tree
  mixes with them.
- **`binary_relocation: false`.** The TeX programs find their tree relative
  to their own location (`SELFAUTOPARENT`). The files must stay as
  released.
- **`TEXMFVAR` in `texmf.cnf`.** TeX writes the fonts and formats that it
  generates to `TEXMFVAR`. The recipe points it at
  `opt/tinytex/texmf-cache`, a directory that is not part of the package.
  pixi links the package files from its cache, so a program that changed
  them would change the cache.
- **Symbolic links in `bin/`.** Every program of the bundle gets a link in
  `$PREFIX/bin`, with one exception: `tlmgr`. It would change the package
  files.

The tree is at `$CONDA_PREFIX/opt/tinytex`.

Costs `[run]`:

| Item | Value |
| --- | --- |
| Download on the first `pixi install` | 150 to 210 MB, by platform |
| Build time | About 1 min |
| Build files kept in `.pixi/bld/tinytex` | About 1.3 GB |

## 4. Isolation from the home directory

A shell that exports `TEXMFVAR`, `TEXMFHOME` or `TEXMFCONFIG` (an XDG
setup does) overrides `texmf.cnf`. In the first test the environment wrote
a font file to `~/.cache/texlive` `[run]`. The activation now sets the
three variables to directories of the environment:

```toml
[target.unix.activation.env]
TEXMFCONFIG = "$CONDA_PREFIX/opt/tinytex/texmf-config"
TEXMFHOME = "$CONDA_PREFIX/opt/tinytex/texmf-local"
TEXMFVAR = "$CONDA_PREFIX/opt/tinytex/texmf-cache"
```

After the change `kpsewhich` reports directories of the environment for the
three, the generated fonts are under `opt/tinytex/texmf-cache`, and a check
writes no file to the TeX directories of the home directory `[run]`.

## 5. Procedures

**Update TinyTeX.** In the recipe set `version` to a release of
`tinytex-releases` and replace each `sha256` with the digest that the
release page shows for the file. Then `pixi lock`, `pixi install` and a
check with the manual.

**A LaTeX package is missing.** `tlmgr` is not on `PATH`, on purpose. The
routes are a newer release, or a change to the recipe that adds the package
files to the tree. Do not install into the environment by hand: the next
`pixi install` does not keep it.

**Test the manual alone.**

```bash
pixi run R CMD Rd2pdf --no-preview --force --output=/tmp/isoformic.pdf .
```

86 pages on 2026-11-12 `[run]`. The `pixi` workflow runs the same command
on Linux and macOS.

## 6. Limits

- **macOS.** The recipe has a `darwin` source and pixi built the package
  for `osx-arm64` on Linux, but no program of it was executed `[solve]`.
- **Windows.** No LaTeX. The Windows `check` task uses `--no-manual`.
- **Linux x86-64 only** among the Linux platforms. `linux-aarch64` needs a
  third source entry.
- **A failed download breaks `pixi install`.** TinyTeX is in the default
  environment, so the GitHub release must be reachable on the first
  install. The alternative is a feature of its own (a `tex` environment),
  at the cost of a third environment. Not done.
- **License.** The package is a set of TeX Live components with their own
  licenses; the recipe says `LicenseRef-TeX-Live`. The package is built
  locally and is not published.
