#!/usr/bin/env just
# shellcheck shell=bash

@default:
  just --choose

test:
  #!/usr/bin/env bash
  R -q -e 'devtools::load_all();';
  R -q -e 'devtools::document();';
  R -q -e 'devtools::test();';

check:
  #!/usr/bin/env bash
  R -q -e 'rcmdcheck::rcmdcheck();';

check:
  #!/usr/bin/env bash -i
  conda create -n isoformic-env -c bioconda -c conda-forge r-base r-devtools bioconductor-summarizedexperiment bioconductor-multiassayexperiment
