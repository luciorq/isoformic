#!/usr/bin/env just
# shellcheck shell=bash

package_name := 'isoformic'

github_org := 'luciorq'

@default:
  just --choose

@test:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();styler::style_pkg();';
  R -q -e 'devtools::load_all();devtools::document();';
  R -q -e 'devtools::load_all();devtools::test();';

@check:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'rcmdcheck::rcmdcheck(args="--as-cran");';

@build-pkgdown-website:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'devtools::document();pkgdown::build_site();';
  git add docs/;
  git commit -m "chore: update pkgdown website";
  git push;

# Check if package can be installed on a conda environment
@check-install-conda tag_version='main':
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  conda create -n isoformic-env -y --override-channels -c bioconda -c conda-forge r-base r-devtools r-readr r-rlang r-dplyr r-ggplot2 r-biocmanager r-pak bioconductor-summarizedexperiment bioconductor-multiassayexperiment;
  conda run -n isoformic-env R -q -e 'pak::pkg_install("github::{{ github_org }}/{{ package_name }}@{{ tag_version }},ask=FALSE")';
  conda run -n isoformic-env R -q -e 'utils::packageVersion("{{ package_name }}")';

# Use R package version on the DESCRIPTION file to tag latest commit of the git repo
@git-tag:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  __r_pkg_version="$(R -q --no-echo --silent -e 'suppressMessages({pkgload::load_all()});cat(as.character(utils::packageVersion("{{ package_name }}")));')";
  \builtin echo -ne "Tagging version: ${__r_pkg_version}\n";
  git tag -a "v${__r_pkg_version}" HEAD -m "Version ${__r_pkg_version} released";
  git push --tags;
