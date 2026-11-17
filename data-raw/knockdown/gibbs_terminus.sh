#!/usr/bin/env bash
# Salmon runs for the Terminus comparisons and the quantification stability
# article, from the reads of quantify_knockdown.sh.
#
# Usage: data-raw/knockdown/gibbs_terminus.sh /path/to/workdir
#
# Environment variables:
#   GENCODE_RELEASE  release of the Salmon index of quantify_knockdown.sh
#                    (default 50)
#   THREADS          threads (default 16)
#   TERMINUS         terminus binary (default WORKDIR/envs/terminus/bin/terminus,
#                    Bioconda package terminus)
#   TREETERMINUS     treeterminus binary (default
#                    WORKDIR/tools/TreeTerminus/target/release/treeterminus,
#                    built from source with Rust and GSL)
#
# Writes, in the work directory:
#   quants_terminus/<run>  the six GSE129540 runs, 30 Gibbs draws, with
#                          --dumpEqWeights (Terminus needs the weights)
#   terminus_all/out/      terminus group on every run, then collapse
#   g100/<run>             SRR8870780, SRR8870781 and SRR8870783 with 100
#                          Gibbs draws
#   tt_single/<run>_<n>    TreeTerminus per sample, at 30 and 100 draws
#   tt_mean/               TreeTerminus on the six runs (--mean_inf true,
#                          about 18 GB of memory)
# terminus group writes a 3.8 GB delta.log per sample, deleted once it is no
# longer needed.
set -euo pipefail
work=${1:?usage: gibbs_terminus.sh WORKDIR}
release=${GENCODE_RELEASE:-50}
threads=${THREADS:-16}
cd "$work"
salmon=envs/salmon/bin/salmon
terminus=${TERMINUS:-envs/terminus/bin/terminus}
treeterminus=${TREETERMINUS:-tools/TreeTerminus/target/release/treeterminus}
index=ref/salmon_index_v$release
runs="SRR8870780 SRR8870781 SRR8870782 SRR8870783 SRR8870784 SRR8870785"
quant() { # run draws outdir
  [ -s "$3/quant.sf" ] && return 0
  "$salmon" quant -i "$index" -l A \
    -1 "sub/GSE129540_$1_1.fastq.gz" -2 "sub/GSE129540_$1_2.fastq.gz" \
    -o "$3" -p "$threads" --numGibbsSamples "$2" --dumpEq --dumpEqWeights \
    --seqBias --gcBias >"$3.log" 2>&1
}
mkdir -p quants_terminus g100 tt_single
for r in $runs; do quant "$r" 30 "quants_terminus/$r"; done
for r in SRR8870780 SRR8870781 SRR8870783; do quant "$r" 100 "g100/$r"; done

# Terminus on all six runs
if [ ! -s terminus_all/collapse.done ]; then
  rm -rf terminus_all
  mkdir -p terminus_all
  for r in $runs; do
    "$terminus" group -d "quants_terminus/$r" -o terminus_all/out >"terminus_all/group_$r.log" 2>&1
  done
  dirs=""
  for r in $runs; do dirs="$dirs quants_terminus/$r"; done
  # shellcheck disable=SC2086
  "$terminus" collapse -d $dirs -o terminus_all/out -t "$threads" >terminus_all/collapse.log 2>&1
  for r in $runs; do rm -f "terminus_all/out/$r/delta.log"; done
  touch terminus_all/collapse.done
fi

# TreeTerminus per sample at 30 and 100 draws
for r in SRR8870780 SRR8870781 SRR8870783; do
  for n in 30 100; do
    d=$([ "$n" = 30 ] && echo "quants_terminus/$r" || echo "g100/$r")
    out=tt_single/${r}_$n
    [ -s "$(find "$out" -name group_nwk.txt 2>/dev/null | head -1)" ] && continue
    rm -rf "$out"
    "$treeterminus" group -d "$d" -o "$out" --mean_inf false >"$out.log" 2>&1
    find "$out" -name delta.log -delete
  done
done

# TreeTerminus on the six runs (mean inferential variance)
if [ ! -s tt_mean/cluster_nwk.txt ]; then
  rm -rf tt_samples tt_mean
  mkdir -p tt_samples
  for r in $runs; do ln -s "../quants_terminus/$r" "tt_samples/$r"; done
  "$treeterminus" group -d tt_samples -o tt_mean --mean_inf true >tt_mean.log 2>&1
  find tt_mean -name delta.log -delete
fi
echo GIBBS_TERMINUS_DONE
