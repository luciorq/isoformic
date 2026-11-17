#!/usr/bin/env bash
# Simulate long reads and quantify them with Oarfish (30 bootstrap replicates).
#
# Usage: data-raw/oarfish/run_oarfish.sh /path/to/workdir
#
# Needs micromamba and Rscript. Installs Oarfish 0.10.3 into WORKDIR/envs,
# simulates three samples from inst/extdata/example_transcripts.fa.gz, runs
# Oarfish on them and copies the output to inst/extdata/oarfish/. Run it from
# the package root.
set -euo pipefail
work=${1:?usage: run_oarfish.sh WORKDIR}
mkdir -p "$work"/{envs,sim,out}
if [ ! -x "$work/envs/oarfish/bin/oarfish" ]; then
  micromamba create -y -p "$work/envs/oarfish" -c conda-forge -c bioconda oarfish=0.10.3
fi
Rscript data-raw/oarfish/simulate_long_reads.R "$work/sim"
zcat inst/extdata/example_transcripts.fa.gz >"$work/sim/transcripts.fa"
for s in control_1 control_2 treated_1; do
  "$work/envs/oarfish/bin/oarfish" --reads "$work/sim/$s.fastq.gz" \
    --annotated "$work/sim/transcripts.fa" --seq-tech ont-cdna \
    --num-bootstraps 30 -j 4 -o "$work/out/$s"
done
mkdir -p inst/extdata/oarfish
cp "$work"/out/*.quant "$work"/out/*.ambig_info.tsv "$work"/out/*.meta_info.json \
  "$work"/out/*.infreps.pq inst/extdata/oarfish/
