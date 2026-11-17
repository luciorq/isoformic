#!/usr/bin/env bash
# Quantify the long-read case study (GEO GSE259334, BioProject PRJNA1081609)
# with Oarfish.
#
# Usage: data-raw/case_studies/quantify_prostate_ont.sh /path/to/workdir
#
# Oxford Nanopore PCR-cDNA sequencing (SQK-PCS109, PromethION) of 22Rv1
# prostate cancer cells treated for 72 h with 50 nM lorecivivint, a CLK
# inhibitor, or with DMSO (3 + 3 samples; Van Goubergen et al. 2024, Mol
# Oncol). Oarfish aligns the reads to the GENCODE transcripts with its
# built-in minimap2 and writes 30 bootstrap replicates per sample.
#
# Environment variables:
#   RUNS     space-separated run accessions to process (default: all six)
#   THREADS  threads for Oarfish and fasterq-dump (default 16)
#   GENCODE_RELEASE      GENCODE human release (default 50)
#   GENCODE_TRANSCRIPTS  local gencode.v<release>.transcripts.fa.gz (default:
#                        download)
#
# Reads come from NCBI SRA (see fetch_sra.sh). Needs micromamba and curl and about 45 GB of free space in the work
# directory (the reads of one run at a time; the SRA file is 9 GB and the
# FASTQ about 25 GB). data-raw/case_studies/build_case_studies.R then reads
# prostate_ont/quants/ from the work directory.
set -euo pipefail

work=${1:?usage: quantify_prostate_ont.sh WORKDIR}
threads=${THREADS:-16}
release=${GENCODE_RELEASE:-50}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root="$work/prostate_ont"
mkdir -p "$work/envs" "$root"/{ref,reads,quants}

# 1. Tools ---------------------------------------------------------------------
if [ ! -x "$work/envs/oarfish/bin/oarfish" ]; then
  micromamba create -y -p "$work/envs/oarfish" -c conda-forge -c bioconda oarfish=0.11.0
fi
if [ ! -x "$work/envs/sra-tools/bin/fasterq-dump" ]; then
  micromamba create -y -p "$work/envs/sra-tools" -c conda-forge -c bioconda sra-tools=3.4.1
fi
oarfish="$work/envs/oarfish/bin/oarfish"
fasterq_dump="$work/envs/sra-tools/bin/fasterq-dump"
# shellcheck source=fetch_sra.sh
. "$script_dir/fetch_sra.sh"

# 2. GENCODE transcripts -------------------------------------------------------
gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release
transcripts="$root/ref/gencode.v$release.transcripts.fa"
if [ ! -s "$transcripts" ]; then
  # Keep only the transcript ID of the GENCODE header, as `salmon --gencode`
  # does, so that the Oarfish output matches the isoformic annotation.
  if [ -n "${GENCODE_TRANSCRIPTS:-}" ]; then
    src_cmd=(zcat "$GENCODE_TRANSCRIPTS")
  else
    src_cmd=(sh -c "curl -sSL '$gencode/gencode.v$release.transcripts.fa.gz' | zcat")
  fi
  "${src_cmd[@]}" |
    sed 's/^>\([^|]*\)|.*/>\1/' >"$transcripts.tmp"
  mv "$transcripts.tmp" "$transcripts"
fi

# 3. Reads and quantification, one run at a time ------------------------------
runs=${RUNS:-$(awk -F'\t' 'NR > 1 { print $1 }' "$script_dir/prostate_ont_samples.tsv")}
for run in $runs; do
  [ -s "$root/quants/$run.quant" ] && continue
  fetch_sra "$run" "$root/reads/$run" "$threads"
  "$oarfish" --reads "$root/reads/$run/$run.fastq" \
    --annotated "$transcripts" --seq-tech ont-cdna \
    --num-bootstraps 30 -j "$threads" -o "$root/quants/$run.tmp"
  for f in "$root/quants/$run.tmp".*; do
    mv "$f" "$root/quants/$run${f#"$root/quants/$run.tmp"}"
  done
  rm -rf "$root/reads/$run"
done
echo "Done: ISOFORMIC_CASE_STUDIES_DIR=$work Rscript data-raw/case_studies/build_case_studies.R"
