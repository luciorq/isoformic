#!/usr/bin/env bash
# Split-half experiment: quantify the first and second 6M read pairs of each
# GSE129540 run (needs the output of quantify_knockdown.sh in WORKDIR).
# Usage: data-raw/knockdown/split_halves.sh /path/to/workdir
# GENCODE_RELEASE selects the Salmon index of quantify_knockdown.sh (default 50).
set -euo pipefail
work=${1:?usage: split_halves.sh WORKDIR}
release=${GENCODE_RELEASE:-50}
cd "$work"
mkdir -p halves quants_half
for run in SRR8870780 SRR8870781 SRR8870782 SRR8870783 SRR8870784 SRR8870785; do
  for mate in 1 2; do
    f=sub/GSE129540_${run}_${mate}.fastq.gz
    [ -s halves/${run}_h1_${mate}.fastq.gz ] || pigz -dc $f | head -n 24000000 | pigz -p 4 > halves/${run}_h1_${mate}.fastq.gz
    [ -s halves/${run}_h2_${mate}.fastq.gz ] || pigz -dc $f | tail -n +24000001 | pigz -p 4 > halves/${run}_h2_${mate}.fastq.gz
  done
  for h in h1 h2; do
    [ -s quants_half/${run}_$h/quant.sf ] && continue
    envs/salmon/bin/salmon quant -i ref/salmon_index_v$release -l A -1 halves/${run}_${h}_1.fastq.gz -2 halves/${run}_${h}_2.fastq.gz \
      -o quants_half/${run}_$h -p 16 --numGibbsSamples 30 --dumpEq --seqBias --gcBias > quants_half/${run}_$h.log 2>&1
  done
done
echo HALVES_DONE
