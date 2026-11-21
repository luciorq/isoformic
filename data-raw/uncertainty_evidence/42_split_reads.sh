#!/usr/bin/env bash
# 42_split_reads.sh: random read halves of the 12 ESRP knockdown runs (WP-C2).
#
# Purpose: split each 12 M read-pair FASTQ pair of GSE129540 (SRR8870780-85)
# and GSE75489 (SRR2966412-17) into two random, disjoint halves A and B. Each
# 4-line record is assigned with gawk's rand() after srand(20261117); the two
# mates are processed by the same program with the same seed, so record i of
# mate 1 and record i of mate 2 receive the same draw and the pair stays
# together. The same seed is used for every run (record i goes to the same
# half in every run; records of different runs are independent, so this only
# fixes the pattern). This is a technical replicate design: the two halves of
# one library differ only by read sampling and by the EM, the quantity that
# InfRV and the Gibbs draws estimate. It is not proof of quantification error.
#
# Inputs (environment variables):
#   KD_SUB_DIR   default $ISOFORMIC_KNOCKDOWN_DIR/sub, the 12 M read-pair
#                subsamples of data-raw/knockdown/
#                (GSE*_SRR*_{1,2}.fastq.gz, 12,000,000 pairs each)
#   RHALF_DIR    default $ISOFORMIC_UE_OUT/rhalf
#   One of KD_SUB_DIR or ISOFORMIC_KNOCKDOWN_DIR, and one of RHALF_DIR or
#   ISOFORMIC_UE_OUT, must be set.
#   SPLIT_SEED   default 20261117
#   RUNS         optional space-separated subset of the 12 run entries
#
# Outputs (under RHALF_DIR):
#   fastq/<entry>_{A,B}_{1,2}.fastq.gz   about 18 GB in all (estimate)
#   split_check.tsv  one row per run: records per file, pairwise name
#                    agreement, disjointness and coverage of the halves
#   split_reads.log  timings and tool versions
#
# Compute: one run at a time; pigz -dc, gawk and two single-thread pigz
# writers, so about 4 cores; no lock needed (under the 4-thread limit). The
# check pass runs 4 decompressors at a time. Memory is small (streams;
# well under 1 GB). Measured: 5 to 9 min per run to split and about 1 min
# to check, 1.5 h for the 12 runs. Needs pigz and gawk.
#
# Run from anywhere:
#   bash data-raw/uncertainty_evidence/42_split_reads.sh
set -euo pipefail

if [ -z "${KD_SUB_DIR:-}" ] && [ -z "${ISOFORMIC_KNOCKDOWN_DIR:-}" ]; then
  echo "Set KD_SUB_DIR or ISOFORMIC_KNOCKDOWN_DIR (see README.md)" >&2
  exit 1
fi
if [ -z "${RHALF_DIR:-}" ] && [ -z "${ISOFORMIC_UE_OUT:-}" ]; then
  echo "Set RHALF_DIR or ISOFORMIC_UE_OUT (see README.md)" >&2
  exit 1
fi
sub=${KD_SUB_DIR:-$ISOFORMIC_KNOCKDOWN_DIR/sub}
out=${RHALF_DIR:-$ISOFORMIC_UE_OUT/rhalf}
seed=${SPLIT_SEED:-20261117}
runs=${RUNS:-"GSE129540_SRR8870780 GSE129540_SRR8870781 GSE129540_SRR8870782
GSE129540_SRR8870783 GSE129540_SRR8870784 GSE129540_SRR8870785
GSE75489_SRR2966412 GSE75489_SRR2966413 GSE75489_SRR2966414
GSE75489_SRR2966415 GSE75489_SRR2966416 GSE75489_SRR2966417"}
mkdir -p "$out/fastq"
log=$out/split_reads.log
{
  echo "date: $(date -Is)"
  echo "seed: $seed"
  gawk --version | head -1
  pigz --version 2>&1
} >>"$log"

split_one() {
  local entry=$1 mate=$2
  local in=$sub/${entry}_${mate}.fastq.gz
  local a=$out/fastq/${entry}_A_${mate}.fastq.gz
  local b=$out/fastq/${entry}_B_${mate}.fastq.gz
  if [ -s "$a" ] && [ -s "$b" ]; then return 0; fi
  pigz -dc -p 1 "$in" |
    gawk -v seed="$seed" -v fa="pigz -p 1 > '$a.tmp'" -v fb="pigz -p 1 > '$b.tmp'" '
      BEGIN { srand(seed) }
      NR % 4 == 1 { h = (rand() < 0.5) ? 1 : 2 }
      { if (h == 1) print | fa; else print | fb }
      END { close(fa); close(fb) }'
  mv "$a.tmp" "$a"
  mv "$b.tmp" "$b"
}

# Spot number of each record header (first token, after the last dot).
names() { pigz -dc -p 1 "$1" | gawk 'NR % 4 == 1 { sub(/^@/, "", $1); n = split($1, p, "."); print p[n] }'; }

check_one() {
  local entry=$1 tmp
  tmp=$(mktemp -d "$out/check.XXXXXX")
  names "$sub/${entry}_1.fastq.gz" >"$tmp/o1" &
  names "$out/fastq/${entry}_A_1.fastq.gz" >"$tmp/a1" &
  names "$out/fastq/${entry}_B_1.fastq.gz" >"$tmp/b1" &
  names "$out/fastq/${entry}_A_2.fastq.gz" >"$tmp/a2" &
  wait
  names "$out/fastq/${entry}_B_2.fastq.gz" >"$tmp/b2"
  local n_o n_a1 n_b1 n_a2 n_b2 pairA pairB inc dup cover
  n_o=$(wc -l <"$tmp/o1")
  n_a1=$(wc -l <"$tmp/a1"); n_b1=$(wc -l <"$tmp/b1")
  n_a2=$(wc -l <"$tmp/a2"); n_b2=$(wc -l <"$tmp/b2")
  pairA=$(cmp -s "$tmp/a1" "$tmp/a2" && echo TRUE || echo FALSE)
  pairB=$(cmp -s "$tmp/b1" "$tmp/b2" && echo TRUE || echo FALSE)
  # Original spot numbers strictly increasing (so a numeric merge is valid).
  inc=$(gawk 'NR > 1 && $1 <= prev { bad++ } { prev = $1 } END { print bad + 0 }' "$tmp/o1")
  # Disjoint: no spot in both halves; cover: merge of halves equals original.
  dup=$(sort -m -n "$tmp/a1" "$tmp/b1" | uniq -d | wc -l)
  cover=$(sort -m -n "$tmp/a1" "$tmp/b1" | cmp -s - "$tmp/o1" && echo TRUE || echo FALSE)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$entry" "$n_o" "$n_a1" "$n_b1" \
    "$n_a2" "$n_b2" "$pairA" "$pairB" "$inc" "$dup" "$cover" >>"$out/split_check.tsv"
  rm -rf "$tmp"
}

printf 'entry\tn_original\tn_A_1\tn_B_1\tn_A_2\tn_B_2\tnames_agree_A\tnames_agree_B\tnon_increasing_original\tspots_in_both_halves\thalves_cover_original\n' \
  >"$out/split_check.tsv"
for entry in $runs; do
  t0=$(date +%s)
  split_one "$entry" 1
  split_one "$entry" 2
  t1=$(date +%s)
  check_one "$entry"
  t2=$(date +%s)
  echo "$entry split_s=$((t1 - t0)) check_s=$((t2 - t1))" >>"$log"
done
echo "SPLIT_DONE $(date -Is)" >>"$log"
