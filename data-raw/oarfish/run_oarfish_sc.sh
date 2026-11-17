#!/usr/bin/env bash
# Simulate barcoded long reads and quantify them with Oarfish --single-cell.
#
# Usage: data-raw/oarfish/run_oarfish_sc.sh /path/to/workdir
#
# Needs micromamba and Rscript; installs Oarfish 0.10.3, minimap2 and samtools
# into WORKDIR/envs. The reads are aligned to the bundled transcripts with
# minimap2, the barcode of each read name is moved to the CB:Z tag, and the
# output is copied to inst/extdata/oarfish_sc/ with the table of cells. Run it
# from the package root.
set -euo pipefail
work=${1:?usage: run_oarfish_sc.sh WORKDIR}
mkdir -p "$work"/{envs,sc}
env="$work/envs/oarfish"
if [ ! -x "$env/bin/oarfish" ]; then
  micromamba create -y -p "$env" -c conda-forge -c bioconda oarfish=0.10.3 minimap2 samtools
fi
Rscript data-raw/oarfish/simulate_single_cell.R "$work/sc"
zcat inst/extdata/example_transcripts.fa.gz >"$work/sc/transcripts.fa"
"$env/bin/minimap2" -t 4 -ax map-ont -N 100 "$work/sc/transcripts.fa" "$work/sc/cells.fastq.gz" |
  awk 'BEGIN { OFS = "\t" } /^@/ { print; next } { split($1, a, "_"); print $0, "CB:Z:" a[1] }' |
  "$env/bin/samtools" view -b -o "$work/sc/cells.bam" -
"$env/bin/oarfish" --single-cell --alignments "$work/sc/cells.bam" -o "$work/sc/out/out" -j 4
mkdir -p inst/extdata/oarfish_sc
cp "$work"/sc/out/out.* "$work/sc/cells.tsv" inst/extdata/oarfish_sc/
