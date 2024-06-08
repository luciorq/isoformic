#!/usr/bin/env just
# shellcheck shell=bash

@default:
  just --choose

test:
  #!/usr/bin/env bash
  set -euxo pipefail;
  R -q -e 'devtools::load_all();';
  R -q -e 'devtools::document();';
  R -q -e 'devtools::test();';

check:
  #!/usr/bin/env bash
  set -euxo pipefail;
  R -q -e 'rcmdcheck::rcmdcheck();';

build-pkgdown-website:
  #!/usr/bin/env bash
  set -euxo pipefail;
  R -q -e 'pkgdown::build_site();';
  git add docs/;
  git commit -m "chore: update pkgdown website";

# Check if package can be installed on a conda environment
check-install-conda:
  #!/usr/bin/env bash -i
  set -euxo pipefail;
  conda create -n isoformic-env -y -c bioconda -c conda-forge r-base r-devtools r-readr r-rlang r-dplyr r-ggplot2 r-biocmanager r-pak bioconductor-summarizedexperiment bioconductor-multiassayexperiment;
  conda run -n isoformic-env R -q -e 'pak::pkg_install("github::luciorq/isoformic@v0.0.1")';
  conda run -n isoformic-env R -q -e 'utils::packageVersion("isoformic")';

# Use R package version on the Description file to tag latest commit of the git repo
git-tag:
  #!/usr/bin/env bash
  # set -euxo pipefail;
  set -euo pipefail;
  __r_pkg_version="$(R -q --no-echo --silent -e 'suppressMessages({pkgload::load_all()});cat(as.character(utils::packageVersion("isoformic")));')";
  builtin echo -ne "Tagging version: ${__r_pkg_version}\n";
  git tag -a "v${__package_version}" HEAD -m "Version ${_r_pkg_version} released";
  git push --tags;
