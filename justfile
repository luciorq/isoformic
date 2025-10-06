#!/usr/bin/env just
# shellcheck shell=bash

package_name := 'isoformic'

github_org := 'luciorq'

@default:
  just --choose

@test:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();styler::style_pkg();';
  R -q -e 'devtools::load_all();usethis::use_tidy_description();';
  R -q -e 'devtools::load_all();devtools::document();';
  R -q -e 'devtools::load_all();devtools::run_examples();';
  R -q -e 'devtools::load_all();devtools::test();';
  # R -q -e 'devtools::load_all();rmarkdown::render("README.Rmd", encoding = "UTF-8")';

@test-all-examples:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();devtools::run_examples(run_dontrun = TRUE, run_donttest = TRUE);';

@check:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'rcmdcheck::rcmdcheck(args = c("--as-cran"), repos = c(CRAN = "https://cloud.r-project.org"));';

# Check if package can be installed on a conda environment
@check-install-conda tag_version='main':
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  # conda create -n isoformic-env -y --override-channels -c bioconda -c conda-forge r-base r-devtools r-readr r-rlang r-dplyr r-ggplot2 r-biocmanager;
  # conda run -n isoformic-env R -q -e 'pak::pkg_install("github::{{ github_org }}/{{ package_name }}@{{ tag_version }},ask=FALSE")';
  # conda run -n isoformic-env R -q -e 'utils::packageVersion("{{ package_name }}")';
  \builtin echo "Not implemented yet";

# Use R package version on the DESCRIPTION file to tag latest commit of the git repo
@git-tag:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  __r_pkg_version="$(R -q --no-echo --silent -e 'suppressMessages({pkgload::load_all()});cat(as.character(utils::packageVersion("{{ package_name }}")));')";
  \builtin echo -ne "Tagging version: ${__r_pkg_version}\n";
  git tag -a "v${__r_pkg_version}" HEAD -m "Version ${__r_pkg_version} released";
  git push --tags;

@pre-release:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'urlchecker::url_check()';
  R -q -e 'devtools::build_readme()';
  R -q -e 'withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org")), {devtools::check(remote = TRUE, manual = TRUE)})';
  R -q -e 'devtools::check_win_devel()';
  # revdepcheck::revdep_check(num_workers = 4)
  # Update CRAN comments
  # usethis::use_version('patch')
  # devtools::build_rmd("vignettes/my-vignette.Rmd")
  # devtools::submit_cran()


@build-vignettes:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();devtools::document();';
  R -q -e 'devtools::install(pkg = ".", build_vignettes = TRUE, dependencies = c("Imports", "Suggests", "Depends"), upgrade = "always");';
  R -q -e 'print(vignette(package = "{{ package_name }}"));';

@build-pkgdown-website:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();devtools::document();pkgdown::build_site();';
  git add docs/;
  git commit -m "chore: update pkgdown website";
  git push;

@release-github:
  #!/usr/bin/env bash
  \builtin set -euxo pipefail;
  # gh release create v0.1.2 --title "v0.1.2 (beta)" --notes "First Zenodo archiving release"
  \builtin echo "Not implemented yet";
