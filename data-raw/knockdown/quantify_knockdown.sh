#!/usr/bin/env bash
# Quantify the ESRP knockdown examples (GSE75489 and GSE129540) with Salmon.
#
# Usage: data-raw/knockdown/quantify_knockdown.sh /path/to/workdir
#
# Environment variables:
#   GENCODE_RELEASE  GENCODE human release of the index (default 50)
#   N_PAIRS          read pairs kept per run (default 12000000)
#   RUNS             space-separated run accessions to process (default: all
#                    twelve)
#   THREADS          threads for Salmon and fasterq-dump (default 16)
#
# Needs micromamba, curl, pigz and about 60 GB of free space in the work
# directory (the decoy-aware Salmon index takes 20 GB; the subsampled reads
# 18 GB; the complete reads of one run at a time up to 20 GB; the
# quantifications with the Gibbs draws and the equivalence classes 1.1 GB).
# Everything is written to the work directory;
# data-raw/build_knockdown_example.R then reads quants_sub/ from it. A
# ref/transcripts.fa.gz or ref/genome.fa.gz already in the work directory is
# used instead of a download.
#
# Only the first 12 million read pairs of each run are used. Reads come from
# NCBI SRA (data-raw/case_studies/fetch_sra.sh); fasterq-dump writes them in
# the order of the submitted files, so the first 12 million pairs are those of
# the FASTQ files on ENA.
set -euo pipefail

work=${1:?usage: quantify_knockdown.sh WORKDIR}
release=${GENCODE_RELEASE:-50}
n_pairs=${N_PAIRS:-12000000}
n_lines=$((4 * n_pairs))
only_runs=${RUNS:-}
threads=${THREADS:-16}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$work"/{envs,ref,sub,reads,quants_sub}
cd "$work"

# 1. Salmon 2.8.0 --------------------------------------------------------------
if [ ! -x envs/salmon/bin/salmon ]; then
  micromamba create -y -p envs/salmon -c conda-forge -c bioconda salmon=2.8.0
fi
if [ ! -x envs/sra-tools/bin/fasterq-dump ]; then
  micromamba create -y -p envs/sra-tools -c conda-forge -c bioconda sra-tools=3.4.1
fi
salmon=envs/salmon/bin/salmon
fasterq_dump=envs/sra-tools/bin/fasterq-dump
# shellcheck source=../case_studies/fetch_sra.sh
. "$script_dir/../case_studies/fetch_sra.sh"

# 2. Decoy-aware index of GENCODE transcripts + GRCh38 -------------------------
gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release
index=ref/salmon_index_v$release
if [ ! -s "$index/.done" ]; then
  [ -s ref/transcripts.fa.gz ] ||
    curl -sSL -o ref/transcripts.fa.gz "$gencode/gencode.v$release.transcripts.fa.gz"
  [ -s ref/genome.fa.gz ] ||
    curl -sSL -o ref/genome.fa.gz "$gencode/GRCh38.primary_assembly.genome.fa.gz"
  zcat ref/genome.fa.gz | sed -n 's/^>\([^ ]*\).*/\1/p' > ref/decoys.txt
  cat ref/transcripts.fa.gz ref/genome.fa.gz > ref/gentrome.fa.gz
  $salmon index -t ref/gentrome.fa.gz -d ref/decoys.txt -i "$index" \
    -p "$threads" --gencode
  rm -f ref/gentrome.fa.gz
  touch "$index/.done"
fi

# 3. Reads: first 12M pairs of every run ----------------------------------------
runs=(
  GSE75489_SRR2966412 GSE75489_SRR2966413 GSE75489_SRR2966414
  GSE75489_SRR2966415 GSE75489_SRR2966416 GSE75489_SRR2966417
  GSE129540_SRR8870780 GSE129540_SRR8870781 GSE129540_SRR8870782
  GSE129540_SRR8870783 GSE129540_SRR8870784 GSE129540_SRR8870785
)
for entry in "${runs[@]}"; do
  run=${entry#*_}
  if [ -n "$only_runs" ] && ! grep -qw -- "$run" <<<"$only_runs"; then
    continue
  fi
  [ -s "sub/${entry}_1.fastq.gz" ] && [ -s "sub/${entry}_2.fastq.gz" ] && continue
  fetch_sra "$run" "reads/$run" "$threads"
  for mate in 1 2; do
    head -n "$n_lines" "reads/$run/${run}_$mate.fastq" | pigz -p 4 >"sub/${entry}_$mate.fastq.gz.tmp"
    [ "$(pigz -dc "sub/${entry}_$mate.fastq.gz.tmp" | wc -l)" -eq "$n_lines" ]
    mv "sub/${entry}_$mate.fastq.gz.tmp" "sub/${entry}_$mate.fastq.gz"
  done
  rm -rf "reads/$run"
done

# 4. Quantify every run with 30 Gibbs draws and the equivalence classes --------
for r1 in sub/*_1.fastq.gz; do
  base=$(basename "$r1" _1.fastq.gz)
  [ -s "quants_sub/$base/quant.sf" ] && continue
  $salmon quant -i "$index" -l A \
    -1 "$r1" -2 "sub/${base}_2.fastq.gz" -o "quants_sub/$base" -p "$threads" \
    --numGibbsSamples 30 --dumpEq --seqBias --gcBias
done
echo "Done: run ISOFORMIC_KNOCKDOWN_DIR=$work Rscript data-raw/build_knockdown_example.R"
