# Where an R package can come from

The routes that were tested for a package in a pixi environment, with the
manifest text that worked. Read `README.md` first. Evidence marks (`[run]`,
`[solve]`, `[docs]`) are defined there.

## 1. Which route

| The package is | Route | Locked |
| --- | --- | --- |
| On conda-forge (most of CRAN) | `r-<name> = "*"` | Yes |
| On bioconda (Bioconductor), Linux or macOS | `bioconductor-<name> = "*"` | Yes |
| On CRAN, not on conda-forge, CRAN dependencies only | Inline git package from `github.com/cran/<name>` (section 2) | Yes, by commit |
| On GitHub or r-universe, CRAN dependencies only | Inline git package (sections 2 and 3) | Yes, by commit |
| Any package with a Bioconductor dependency | A rattler-build recipe (section 5) | Yes |
| A source tarball by URL | Does not work in the manifest (section 4); use a recipe | |
| Anything on Windows that conda lacks | `install.packages()` in a task (`WINDOWS.md`) | No |

conda-forge can be behind CRAN. On 2026-11-11 CRAN had `duckdb` 2.0.0 and
conda-forge had `r-duckdb` 1.5.5. A minimum version in `DESCRIPTION` above
the conda-forge version makes the environment unsolvable, which is why
isoformic keeps `duckdb (>= 1.5.4)`.

## 2. An inline git package

A dependency can carry its own package definition. pixi clones the
repository, the backend reads its `DESCRIPTION` and builds it `[run]`:

<!-- markdownlint-disable MD013 -->

```toml
r-condathis = { git = "https://github.com/luciorq/condathis", rev = "ecc46447b094aba91e5fc1cd61493fa97f0194a0", package = { build.backend = { name = "pixi-build-r", version = "0.2.*" }, host-dependencies = { r-base = "4.5.*" } } }
```

<!-- markdownlint-enable MD013 -->

- The conda name comes from the key (`r-condathis`) and the version from
  `DESCRIPTION`; `version` can be left out.
- `rev`, `tag` or `branch` select the commit. The lock file records the
  commit in every case.
- condathis builds in about 15 s.
- The `host-dependencies` pin is needed: see `MANIFEST.md`, section 4.
- The dependencies of the package must be on the channels or be another
  inline entry. A Bioconductor dependency fails (section 5).

CRAN has a git mirror with one tag per release, so the same form installs a
CRAN package that conda-forge lacks `[run]`:

<!-- markdownlint-disable MD013 -->

```toml
r-fortunes = { git = "https://github.com/cran/fortunes", tag = "1.5-5", package = { build.backend = { name = "pixi-build-r", version = "0.2.*" }, host-dependencies = { r-base = "4.5.*" } } }
```

<!-- markdownlint-enable MD013 -->

## 3. r-universe

pixi cannot use an r-universe repository as a channel. But every r-universe
package is built from a git commit, and the API gives it:

```bash
curl -s https://luciorq.r-universe.dev/api/packages/condathis |
  jq -r '.RemoteUrl, .RemoteSha'
```

Use the two values as `git` and `rev` of an inline package (section 2).
This is how condathis is installed.

## 4. A tarball URL does not work

pixi accepts `url` and `sha256` in a source dependency, and the first steps
work: pixi downloads the tarball and the backend returns the metadata of
the package. Then pixi panics `[run]`:

<!-- markdownlint-disable MD013 -->

```toml
r-fortunes = { url = "https://cran.r-project.org/src/contrib/fortunes_1.5-5.tar.gz", sha256 = "f8c103b826dd626d99179fde869f1997882a71d414b3155682951e4cd2000a2e", package = { build.backend = { name = "pixi-build-r", version = "0.2.*" }, host-dependencies = { r-base = "4.5.*" } } }
```

<!-- markdownlint-enable MD013 -->

```text
crates/pixi_record/src/pinned_source.rs:595:9: not implemented
```

- The limit is in pixi, not in `pixi-build-r`: the code that writes a URL
  source to the lock file is `unimplemented!()`, at tag `v0.84.2` and on
  `main` (read on 2026-11-12).
- Upstream issue:
  [prefix-dev/pixi#7230](https://github.com/prefix-dev/pixi/issues/7230),
  open, reported with the same CRAN example on pixi 0.81.0.
- `[package.build.source]` accepts `git` and `path`, not `url`.
- The same panic came from an r-universe tarball
  (`.../src/contrib/condathis_0.2.0.9000.tar.gz`).

When #7230 is fixed, a CRAN or r-universe tarball with its `sha256` should
become the simplest route for a package that conda-forge lacks. Until then
use git (section 2) or a recipe (section 5).

## 5. A rattler-build recipe

A recipe gives full control: any source (a URL too), any dependency names.
It is the only route for a package with Bioconductor dependencies, because
the requirements are written by hand. Tested with CMScaller, which imports
`Biobase` `[run]`, 24 s:

<!-- markdownlint-disable MD013 -->

```toml
r-cmscaller = { path = "recipes/r-cmscaller", package.build.backend = { name = "pixi-build-rattler-build", version = "0.4.*" } }
```

<!-- markdownlint-enable MD013 -->

```yaml
# recipes/r-cmscaller/recipe.yaml
package:
  name: r-cmscaller
  version: 2.0.1
source:
  git: https://github.com/Lothelab/CMScaller.git
  rev: f1541730b56c2f37c123037476eed0e6de0ee6f0
build:
  number: 0
  script: R CMD INSTALL --build .
requirements:
  host: [r-base, bioconductor-biobase]
  run: [r-base, bioconductor-biobase]
```

The same with a CRAN tarball `[run]`, 8 s:

```yaml
# recipes/r-fortunes/recipe.yaml
package:
  name: r-fortunes
  version: 1.5.5
source:
  url: https://cran.r-project.org/src/contrib/fortunes_1.5-5.tar.gz
  sha256: f8c103b826dd626d99179fde869f1997882a71d414b3155682951e4cd2000a2e
build:
  number: 0
  script: R CMD INSTALL --build .
requirements:
  host: [r-base 4.5.*]
  run: [r-base]
```

Notes:

- With this backend the binary dependencies go in the recipe, not in
  `pixi.toml` `[docs]`.
- The requirements do not follow `DESCRIPTION`. Update them by hand when
  the package changes its imports.
- The work directory of the build also holds files of the build tool. Give
  a source that is a whole directory tree a `target_directory`
  (`LATEX.md`, section 3).

## 6. Bioconductor

- As a dependency of the environment: `bioconductor-<name>` from bioconda,
  Linux and macOS only. bioconda builds one Bioconductor release per R
  version (3.22 for R 4.5).
- As an `Imports` entry of a package that `pixi-build-r` builds: not
  supported (`MANIFEST.md`, section 2). For the package of the workspace
  this would mean a recipe in place of the backend, and the loss of
  `DESCRIPTION` as the source of the dependencies. isoformic has
  Bioconductor packages in `Suggests` only, so it does not have the
  problem.
- From the Bioconductor git server (`git.bioconductor.org/packages/<pkg>`):
  possible as a recipe source. Not tested.

## 7. Not locked: `install.packages()` in a task

A task can install into the library of the environment with
`install.packages()`. The lock file does not record the result, and
`pixi install` can remove or replace such packages when it updates the
environment. isoformic does this on Windows only, where no locked route
exists (`WINDOWS.md`).
