#!/usr/bin/env just
# shellcheck shell=bash

package_name := 'isoformic'

github_org := 'luciorq'

@default:
  just --choose

# >>> rstats-package-dev-tasks >>>

# =============================================================================
# General R Package Development Tasks
# =============================================================================

# Update Package Documentation
@document:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);usethis::use_tidy_description();';
  R -q -s -e 'devtools::load_all(quiet=TRUE);devtools::document(quiet=TRUE);';
  \builtin echo "Documentation updated!";

# Lint R Package Code and Documentation
@lint:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  air format ./R/ ./tests/ ./data-raw/;
  # jarl rules are configured in `jarl.toml`; `--fix` applies the safe fixes.
  jarl check . --fix --allow-dirty || jarl check . --statistics;
  R -q -s -e 'lints <- lintr::lint_package(); print(lints); if (length(lints) > 0) quit(status = 1)';
  just document;
  # Fail on non-ASCII characters (em dashes, curly quotes, ...) in package
  # + sources. Non-ASCII outside comments breaks `R CMD check --as-cran`
  # + ("Portable packages must use only ASCII characters"). README.md and the
  # + articles are excluded on purpose: they contain author names and
  # + rendered output.
  if \grep -rnI -P '[^\x00-\x7F]' ./R/ ./tests/ ./man/ NAMESPACE NEWS.md; then
    \builtin echo "ERROR: non-ASCII characters found (listed above). Replace them with ASCII equivalents.";
    \builtin exit 1;
  fi
  # Remove duplicate lines from .gitignore and .Rbuildignore without changing their order
  awk '!seen[$0]++' .gitignore > .gitignore.tmp && \mv .gitignore.tmp .gitignore;
  awk '!seen[$0]++' .Rbuildignore > .Rbuildignore.tmp && \mv .Rbuildignore.tmp .Rbuildignore;
  lychee --no-progress . || \builtin echo 'lychee not installed or links failed.';
  \builtin echo "Linting done!";

# Run All Unit Tests
@test:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'withr::with_envvar(new=list(`TESTTHAT_CPUS`=4),code={withr::with_options(list(crayon.enabled=FALSE,cli.unicode=FALSE),code={devtools::load_all(quiet=TRUE);devtools::test(pkg=".");})});';
  \builtin echo "All tests passed!";

# Run All Unit Tests Including The Full-Data Tests (downloads GENCODE, runs Salmon)
@test-full:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  ISOFORMIC_FULL=true R -q -s -e 'withr::with_options(list(crayon.enabled=FALSE,cli.unicode=FALSE),code={devtools::load_all(quiet=TRUE);devtools::test(pkg=".");})';
  \builtin echo "All tests passed!";

# Run Documentation Examples
@examples:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);devtools::run_examples();';

# Format R Code With air
@format:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  air format ./R/ ./tests/ ./data-raw/;

# Regenerate Bundled Example Data Under inst/extdata
@extdata:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  Rscript data-raw/prepare_extdata.R;

# Retrain The Coding Potential Model Under inst/extdata/coding_model
@coding-model:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  Rscript data-raw/build_coding_model.R;

# Rebuild The Frozen Evidence Datasets (all, or `just frozen-data appris`)
@frozen-data *ARGS:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  Rscript data-raw/build_frozen_data.R {{ARGS}};

# Assemble The Zenodo Deposit Bundle Of Frozen Files (`just zenodo-deposit [dir]`)
@zenodo-deposit *ARGS:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  Rscript data-raw/zenodo_deposit.R {{ARGS}};

# Render Articles Locally With Quarto
@render-articles:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);pkgdown::build_articles();';

# Render Articles Locally Running The GENCODE-Based Code Too
@render-articles-full:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  ISOFORMIC_FULL=true R -q -s -e 'devtools::load_all(quiet=TRUE);pkgdown::build_articles();';

# Build and Lint README File
@build-readme: lint
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  # Lint markdown files
  [[ -f ./README.Rmd ]] && cat ./README.Rmd | rumdl check --stdin --disable 'MD046' || true;
  [[ -f ./README.qmd ]] && cat ./README.qmd | rumdl check --stdin --disable 'MD046' || true;
  just install-deps;
  R -q -s -e 'pak::local_install(upgrade=TRUE, dependencies=TRUE);';
  # R -q -s -e 'devtools::install(pkg = ".", build_vignettes = TRUE, dependencies = c("Imports", "Suggests", "Depends"), upgrade = "always");';
  [[ -f ./README.Rmd ]] && R -q -s -e 'devtools::load_all(quiet=TRUE);if(file.exists("README.Rmd"))rmarkdown::render("README.Rmd", encoding = "UTF-8")' || true;
  [[ -f ./README.qmd ]] && quarto render README.qmd --to gfm || true;
  # Lint Final README.md
  sed -i -e 's|[[:blank:]]*$||g' README.md;
  sed -i -e '/./,$!d' README.md;
  sed -i '/<!-- badges: start -->/{n; /^\s*$/d}' README.md;
  rumdl check README.md || true;
  markdownlint README.md || true;
  markdownlint-cli2 --no-globs README.md || true;
  \builtin echo "README built and linted!";

# Run All Examples in the Documentation Including `dontrun`
@test-all-examples: document
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);devtools::document();devtools::run_examples(run_dontrun = TRUE, run_donttest = TRUE);';

# Run Tests from a Specific Test File
@test-file file_name:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);devtools::test_active_file("tests/testthat/test-{{ file_name }}.R")';

# Run R CMD Check on the Package With Cran Like Checks
@check: lint test examples
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  if [[ $(\builtin command -v qpdf) ]]; then
    R -q -s -e 'withr::with_envvar(new=list(`_R_CHECK_SYSTEM_CLOCK_`="0"),code={rcmdcheck::rcmdcheck(args = c("--as-cran"), repos = c(CRAN = "https://cloud.r-project.org"));});';
  else
    \builtin echo -ne '\n\nNOTE: LaTeX and PDF handling tools missing, skipping the manual.\n';
    \builtin echo -ne '  - Try `R -q -s -e "tinytex::install_tinytex()"` and `pixi global install qpdf ghostscript`.\n\n';
    R -q -s -e 'withr::with_envvar(new=list(`_R_CHECK_SYSTEM_CLOCK_`="0"),code={rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"), repos = c(CRAN = "https://cloud.r-project.org"));});';
  fi

# Force GitHub Actions Checks to start for the main branch
@check-gha-trigger:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  gh workflow run "r-cmd-check" --ref main;

# Print latest GitHub Actions Checks results for the main branch
@monitor-gha:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  gh run list;
  latest_job_id="$(gh run list -w "r-cmd-check" --json databaseId --jq '.[0].databaseId')";
  gh run view "${latest_job_id}";

# Use R package version on the DESCRIPTION file to tag latest commit of the git repo
@git-tag:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  git pull origin --tags || true;
  git pull upstream --tags || true;
  __r_pkg_version="$(R -q --no-echo --silent -e 'suppressMessages({pkgload::load_all(quiet=TRUE)});cat(as.character(utils::packageVersion("{{ package_name }}")));')";
  \builtin echo -ne "Tagging version: ${__r_pkg_version}\n";
  git tag -a "v${__r_pkg_version}" HEAD -m "Version ${__r_pkg_version} released";
  # git push --tags;
  # git pull upstream --tags;
  # git push upstream --tags;

# Build Vignettes
@build-vignettes:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'devtools::load_all(quiet=TRUE);devtools::document();';
  just install-deps;
  just install-local;
  # TODO: `devtools::install()` changed completely its signature by using
  # + `pak::local_install_deps()` for managing dependencies.
  # R -q -s -e 'devtools::install(pkg = ".", build_vignettes = TRUE, dependencies = c("Imports", "Suggests", "Depends"), upgrade = "always");';
  R -q -s -e 'devtools::install(pkg = ".", build_vignettes = TRUE, dependencies = TRUE, upgrade = TRUE);';
  R -q -s -e 'print(vignette(package = "{{ package_name }}"));';

# Install Package Development Dependencies Including Suggests
@install-deps:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -e 'if(!requireNamespace("pak", quietly=TRUE)) {install.packages("pak")};';
  R -q -e 'pak::local_install_dev_deps(upgrade=TRUE, dependencies=TRUE);';

# Install Local Version of the Package With All Dependencies
@install-local:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'if(!requireNamespace("pak", quietly=TRUE)) {install.packages("pak")};';
  R -q -s -e 'pak::local_install(upgrade=TRUE, dependencies=TRUE);';

# Install Dev Environment Dependencies
@install-dev-deps: install-deps
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'pak::local_install_dev_deps(dependencies="Config/Needs/dev", upgrade=TRUE)';

# Build the pkgdown Website
@build-pkgdown-website: install-deps
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  R -q -s -e 'pak::pak("pkgdown", upgrade=TRUE, dependencies=TRUE);';
  R -q -s -e 'pkgdown::build_favicons(overwrite=FALSE);';
  R -q -s -e 'devtools::document(quiet=TRUE);devtools::load_all(quiet=TRUE);pkgdown::build_site();';
  # Steps for manually deploying the pkgdown website,
  # + not necessary if using GitHub Actions.
  # git add docs/;
  # git add pkgdown/;
  # git add _pkgdown.yml;
  # git commit -m "chore: update pkgdown website";
  # git push;

# Release a New Version on GitHub Releases Using Content From NEWS.md
@release-github:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  # Using `gh` cli to create a new release
  # gh release create v0.1.0 --title "v0.1.0" --notes "First Zenodo archiving release"
  # Using `usethis` to create a new release
  # First, change the NEWS.md header temporarily to make usethis happy
  # + Check: <https://github.com/r-lib/pkgdown/issues/2897>
  # + <https://github.com/r-lib/usethis/issues/2130>
  # + <https://github.com/tidyverse/style/issues/245#issuecomment-2959459898>
  sed -i -e "s|^## {{ package_name }}|# {{ package_name }}|g" NEWS.md;
  R -q -s -e 'devtools::load_all(quiet=TRUE);usethis::use_github_release();';
  sed -i -e "s|^# {{ package_name }}|## {{ package_name }}|g" NEWS.md;
  \builtin echo "Check the GH Releases!";

# Steps to Run Before Releasing a New Version
@pre-release:
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  # TODO: Update version on NEWS.md and DESCRIPTION files.
  # TODO: If using changelog links, also perform a git push and tag,
  # + so `url_check()` do not fail.
  # git commit -m "chore: prepare for vX.Y.Z release";
  # git push;
  # just git-tag;
  # git push --tags;
  # git push upstream --tags;
  R -q -e 'urlchecker::url_check()';
  # R -q -e 'devtools::build_readme()';
  just build-readme;
  R -q -e 'withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org")), {devtools::check(remote = TRUE, manual = TRUE)})';
  R -q -e 'devtools::check_win_devel()';
  # R -q -s -e 'if(!requireNamespace("revdepcheck", quietly=TRUE)) pak::pak("r-lib/revdepcheck");';
  # R -q -e 'revdepcheck::revdep_reset();';
  # R -q -e 'revdepcheck::revdep_check(num_workers = 4);';
  # TODO: Update CRAN comments with output from `check`.
  # cat ./cran-comments.md;
  # usethis::use_version('patch')
  # devtools::build_rmd("vignettes/my-vignette.Rmd")
  # just build-vignettes;
  \builtin echo "Pre-release checks done!";
  # TODO: This step needs to be run manually, from R.
  # + Wait for previous `check_win_devel()` to finish.
  # devtools::submit_cran();
  # TODO: Check your email! Click the link, and check all boxes!
  # TODO: When accepted run:
  # just release-github;
  # TODO: Update Zenodo archive with the new release.
  # + Update the DOI badge in the README file.
  \builtin echo "Release successful!";

# <<< rstats-package-dev-tasks <<<

# =============================================================================
# Isoformic Specific Tasks
# =============================================================================

# Check if package can be installed on a conda environment
@check-install-conda tag_version='main':
  #!/usr/bin/env bash
  \builtin set -euo pipefail;
  # micromamba create -n isoformic-env -y --override-channels -c bioconda -c conda-forge r-base r-devtools r-readr r-rlang r-dplyr r-ggplot2 r-biocmanager;
  # micrmmamba run -n isoformic-env R -q -e 'pak::pkg_install("github::{{ github_org }}/{{ package_name }}@{{ tag_version }},ask=FALSE")';
  # micromamba run -n isoformic-env R -q -e 'utils::packageVersion("{{ package_name }}")';
  \builtin echo "Not implemented yet";

