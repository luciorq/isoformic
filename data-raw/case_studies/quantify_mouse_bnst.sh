#!/usr/bin/env bash
# Quantify the mouse chronic stress case study (BioProject PRJNA984605) with
# Salmon.
#
# Usage: data-raw/case_studies/quantify_mouse_bnst.sh /path/to/workdir
#
# Bulk RNA-seq (paired end) of the anterodorsal bed nucleus of the stria
# terminalis (adBNST) of control and chronically stressed mice (Degroat et
# al. 2024, Psychoneuroendocrinology). The article uses the 19 male samples
# (10 control, 9 stressed).
#
# Environment variables:
#   RUNS     space-separated run accessions to process (default: the 19 males)
#   THREADS  threads for Salmon and fasterq-dump (default 16)
#   GENCODE_RELEASE      GENCODE mouse release (default M39)
#   GENCODE_TRANSCRIPTS  local gencode.v<release>.transcripts.fa.gz (default:
#                        download)
#
# Reads come from NCBI SRA (see fetch_sra.sh). Needs micromamba and curl and
# about 60 GB of free space in the work directory, so keep it off small home
# drives (the decoy-aware Salmon index takes 20 GB, the reads of one run at a
# time up to 40 GB, the quantifications with the Gibbs draws about 3 GB). data-raw/case_studies/build_case_studies.R then
# reads mouse_bnst/quants/ from the work directory.
set -euo pipefail

work=${1:?usage: quantify_mouse_bnst.sh WORKDIR}
threads=${THREADS:-16}
release=${GENCODE_RELEASE:-M39}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root="$work/mouse_bnst"
mkdir -p "$work/envs" "$root"/{ref,reads,quants}

# 1. Tools ---------------------------------------------------------------------
if [ ! -x "$work/envs/salmon/bin/salmon" ]; then
  micromamba create -y -p "$work/envs/salmon" -c conda-forge -c bioconda salmon=2.11.0
fi
if [ ! -x "$work/envs/sra-tools/bin/fasterq-dump" ]; then
  micromamba create -y -p "$work/envs/sra-tools" -c conda-forge -c bioconda sra-tools=3.4.1
fi
salmon="$work/envs/salmon/bin/salmon"
fasterq_dump="$work/envs/sra-tools/bin/fasterq-dump"
# shellcheck source=fetch_sra.sh
. "$script_dir/fetch_sra.sh"

# 2. Decoy-aware index of GENCODE transcripts + GRCm39 -------------------------
# The genome only provides the decoy sequences, so it comes from UCSC (mm39,
# same chromosome names as GENCODE), which is much faster than the EBI FTP.
# GENCODE_TRANSCRIPTS can point to a local copy of the transcript FASTA.
gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_$release
index="$root/ref/salmon_index_v$release"
if [ ! -s "$index/.done" ]; then
  tx_fa="$root/ref/gencode.v$release.transcripts.fa.gz"
  if [ -n "${GENCODE_TRANSCRIPTS:-}" ]; then
    cp "$GENCODE_TRANSCRIPTS" "$tx_fa"
  elif [ ! -s "$tx_fa" ]; then
    curl -sSL -o "$tx_fa.tmp" "$gencode/gencode.v$release.transcripts.fa.gz"
    mv "$tx_fa.tmp" "$tx_fa"
  fi
  curl -sSL -o "$root/ref/genome.fa.gz" \
    https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
  zcat "$root/ref/genome.fa.gz" | sed -n 's/^>\([^ ]*\).*/\1/p' >"$root/ref/decoys.txt"
  cat "$tx_fa" "$root/ref/genome.fa.gz" >"$root/ref/gentrome.fa.gz"
  "$salmon" index -t "$root/ref/gentrome.fa.gz" -d "$root/ref/decoys.txt" \
    -i "$index" -p "$threads" --gencode
  rm -f "$root/ref/gentrome.fa.gz" "$root/ref/genome.fa.gz"
  touch "$index/.done"
fi

# 3. Reads and quantification, one run at a time ------------------------------
# mouse_bnst_samples.tsv lists every run of the BioProject with the sex,
# condition and replicate parsed from the sample titles.
runs=${RUNS:-$(awk -F'\t' 'NR > 1 && $2 == "male" { print $1 }' "$script_dir/mouse_bnst_samples.tsv")}
for run in $runs; do
  [ -s "$root/quants/$run/quant.sf" ] && continue
  fetch_sra "$run" "$root/reads/$run" "$threads"
  "$salmon" quant -i "$index" -l A \
    -1 "$root/reads/$run/${run}_1.fastq" -2 "$root/reads/$run/${run}_2.fastq" \
    -o "$root/quants/$run.tmp" -p "$threads" \
    --numGibbsSamples 30 --dumpEq --posBias --seqBias --gcBias
  mv "$root/quants/$run.tmp" "$root/quants/$run"
  rm -rf "$root/reads/$run"
done
echo "Done: ISOFORMIC_CASE_STUDIES_DIR=$work Rscript data-raw/case_studies/build_case_studies.R"
