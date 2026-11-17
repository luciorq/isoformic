#!/usr/bin/env bash
# Quantify the preeclampsia placenta dataset (BioProject PRJNA622166) with
# Salmon, for the bundled example data.
#
# Usage: data-raw/placenta/quantify_placenta.sh /path/to/workdir
#
# Bulk RNA-seq (Illumina HiSeq 2500, paired end) of placentas from 8
# pregnancies with preeclampsia and 17 controls. Every run is quantified with
# 30 Gibbs draws against the decoy-aware GENCODE index of
# data-raw/knockdown/quantify_knockdown.sh (complete quantifications in
# quants/). The six runs of the mini example are also quantified against an
# index of the transcripts of the genes in mini_genes.tsv only, with 20 Gibbs
# draws (mini_quants/): three preeclampsia and three control placentas.
#
# Environment variables:
#   GENCODE_RELEASE  GENCODE human release (default 50)
#   INDEX            decoy-aware Salmon index of that release (default
#                    WORKDIR/ref/salmon_index_v<release>, built by
#                    quantify_knockdown.sh in the same work directory)
#   RUNS             space-separated run accessions to process (default: all
#                    25)
#   THREADS          threads for Salmon and fasterq-dump (default 16)
#
# Reads come from NCBI SRA (data-raw/case_studies/fetch_sra.sh), one run at a
# time (about 15 GB of FASTQ). Needs micromamba, curl and pigz.
# data-raw/prepare_extdata.R then reads quants/ and mini_quants/ from the
# work directory (ISOFORMIC_PLACENTA_DIR).
set -euo pipefail

work=${1:?usage: quantify_placenta.sh WORKDIR}
release=${GENCODE_RELEASE:-50}
threads=${THREADS:-16}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
index=${INDEX:-$work/ref/salmon_index_v$release}
root="$work/placenta"
mkdir -p "$work/envs" "$root"/{ref,reads,quants,mini_quants}

# 1. Tools ---------------------------------------------------------------------
if [ ! -x "$work/envs/salmon/bin/salmon" ]; then
  micromamba create -y -p "$work/envs/salmon" -c conda-forge -c bioconda salmon=2.8.0
fi
if [ ! -x "$work/envs/sra-tools/bin/fasterq-dump" ]; then
  micromamba create -y -p "$work/envs/sra-tools" -c conda-forge -c bioconda sra-tools=3.4.1
fi
salmon="$work/envs/salmon/bin/salmon"
fasterq_dump="$work/envs/sra-tools/bin/fasterq-dump"
# shellcheck source=../case_studies/fetch_sra.sh
. "$script_dir/../case_studies/fetch_sra.sh"
# Salmon 2 writes info.json, Salmon 1 versionInfo.json
[ -s "$index/info.json" ] || [ -s "$index/versionInfo.json" ] || {
  echo "No Salmon index at $index: run quantify_knockdown.sh first" >&2
  exit 1
}

# 2. Index of the transcripts of the mini example genes -------------------------
gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release
tx_fa="$work/ref/transcripts.fa.gz"
[ -s "$tx_fa" ] ||
  curl -sSL -o "$tx_fa" "$gencode/gencode.v$release.transcripts.fa.gz"
mini_index="$root/ref/mini_index_v$release"
if [ ! -s "$mini_index/.done" ]; then
  genes=$(awk -F'\t' 'NR > 1 { print $1 }' "$script_dir/mini_genes.tsv" | paste -sd'|')
  zcat "$tx_fa" | awk -v genes="$genes" '
    /^>/ { split($0, f, "|"); g = f[2]; sub(/\..*/, "", g); keep = (g ~ "^(" genes ")$") }
    keep' >"$root/ref/mini_transcripts.fa"
  "$salmon" index -t "$root/ref/mini_transcripts.fa" -i "$mini_index" -p "$threads" --gencode
  touch "$mini_index/.done"
fi

# 3. Reads and quantification, one run at a time ------------------------------
mini_runs="SRR11498039 SRR11498040 SRR11498041 SRR11498077 SRR11498078 SRR11498079"
all_runs=$(awk -F'\t' 'NR > 1 { print $1 }' "$script_dir/placenta_samples.tsv")
for run in ${RUNS:-$all_runs}; do
  is_mini=$(grep -qw -- "$run" <<<"$mini_runs" && echo 1 || echo 0)
  in_all=$(grep -qw -- "$run" <<<"$all_runs" && echo 1 || echo 0)
  if { [ "$in_all" = 0 ] || [ -s "$root/quants/$run/quant.sf" ]; } &&
    { [ "$is_mini" = 0 ] || [ -s "$root/mini_quants/${run}_quant/quant.sf" ]; }; then
    continue
  fi
  fetch_sra "$run" "$root/reads/$run" "$threads"
  r1="$root/reads/$run/${run}_1.fastq"
  r2="$root/reads/$run/${run}_2.fastq"
  if [ "$in_all" = 1 ] && [ ! -s "$root/quants/$run/quant.sf" ]; then
    "$salmon" quant -i "$index" -l A -1 "$r1" -2 "$r2" \
      -o "$root/quants/$run.tmp" -p "$threads" \
      --numGibbsSamples 30 --dumpEq --seqBias --gcBias
    mv "$root/quants/$run.tmp" "$root/quants/$run"
  fi
  if [ "$is_mini" = 1 ] && [ ! -s "$root/mini_quants/${run}_quant/quant.sf" ]; then
    "$salmon" quant -i "$mini_index" -l A -1 "$r1" -2 "$r2" \
      -o "$root/mini_quants/${run}_quant.tmp" -p "$threads" \
      --numGibbsSamples 20
    mv "$root/mini_quants/${run}_quant.tmp" "$root/mini_quants/${run}_quant"
  fi
  rm -rf "$root/reads/$run"
done
echo "Done: ISOFORMIC_PLACENTA_DIR=$work Rscript data-raw/prepare_extdata.R"
