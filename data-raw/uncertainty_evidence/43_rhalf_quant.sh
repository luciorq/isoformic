#!/usr/bin/env bash
# 43_rhalf_quant.sh: Salmon on the random read halves of 42_split_reads.sh.
#
# Purpose: quantify each half (A, B) of the 12 knockdown runs with the
# settings of data-raw/knockdown/split_halves.sh and quantify_knockdown.sh:
# salmon 2.8.0, decoy-aware GENCODE 50 index, -l A, --numGibbsSamples 30,
# --dumpEq, --seqBias, --gcBias. 24 runs.
#
# Inputs (environment variables):
#   RHALF_DIR    default $ISOFORMIC_UE_OUT/rhalf (one of the two must be set)
#                (reads in fastq/ from 42_split_reads.sh; the script waits for
#                each run's halves to appear, so it can start with 42)
#   ISOFORMIC_V50_DIR  work directory with human/envs/salmon/bin/salmon and
#                human/ref/salmon_index_v50 (required)
#   THREADS      Salmon threads, default 10
#   MEMCAP_WRAPPER  optional command that runs a command under a memory cap,
#                called as `$MEMCAP_WRAPPER <MEM_GB> <command...>` (for
#                example a small script around `systemd-run --user --scope
#                -p MemoryMax=<GB>G`); default: none, Salmon runs directly
#   MEM_GB       memory cap of each Salmon run, default 20 (the 6.2 GB
#                decoy-aware index plus mapping buffers; estimate)
#
# Outputs (under RHALF_DIR):
#   quants/<entry>_<A|B>/        Salmon output with 30 Gibbs draws and eq classes
#   quants/<entry>_<A|B>.log     Salmon stderr plus /usr/bin/time -v
#   quant_times.tsv              one row per attempt: wall time (s, from date),
#                                max RSS (kB), lock wait (s), start, exit status
#                                (0 = verified run), attempt
#   MAX_ATTEMPTS (env, default 8) tries per run before the script stops
#   MIN_AVAIL_GB (env, default 19) effective available memory (MemAvailable
#                plus reclaimable ZFS ARC) needed before each attempt
#
# Measured on 2026-11-18: 533 to 1112 s per half with 10 threads, peak RSS
# about 11.7 GB; 6.5 h in all including waits for other heavy jobs. Run from
# anywhere; the driver loop itself is light (well under 1 GB):
#   bash data-raw/uncertainty_evidence/43_rhalf_quant.sh
set -euo pipefail

if [ -z "${RHALF_DIR:-}" ] && [ -z "${ISOFORMIC_UE_OUT:-}" ]; then
  echo "Set RHALF_DIR or ISOFORMIC_UE_OUT (see README.md)" >&2
  exit 1
fi
if [ -z "${ISOFORMIC_V50_DIR:-}" ]; then
  echo "Set ISOFORMIC_V50_DIR (see README.md)" >&2
  exit 1
fi
out=${RHALF_DIR:-$ISOFORMIC_UE_OUT/rhalf}
v50=$ISOFORMIC_V50_DIR
threads=${THREADS:-10}
mem_gb=${MEM_GB:-20}
wrap=()
if [ -n "${MEMCAP_WRAPPER:-}" ]; then
  wrap=("$MEMCAP_WRAPPER" "$mem_gb")
fi
min_avail=${MIN_AVAIL_GB:-19}
salmon=$v50/human/envs/salmon/bin/salmon
index=$v50/human/ref/salmon_index_v50
runs="GSE129540_SRR8870780 GSE129540_SRR8870781 GSE129540_SRR8870782
GSE129540_SRR8870783 GSE129540_SRR8870784 GSE129540_SRR8870785
GSE75489_SRR2966412 GSE75489_SRR2966413 GSE75489_SRR2966414
GSE75489_SRR2966415 GSE75489_SRR2966416 GSE75489_SRR2966417"
eff_avail_gb() {
  local a z=0
  a=$(gawk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
  [ -r /proc/spl/kstat/zfs/arcstats ] &&
    z=$(gawk '$1 == "size" {s = $3} $1 == "c_min" {m = $3} END {v = int((s - m) / 1048576); print (v > 0 ? v : 0)}' /proc/spl/kstat/zfs/arcstats)
  echo $(((a + z) / 1024))
}
mkdir -p "$out/quants"
times=$out/quant_times.tsv
[ -s "$times" ] || printf 'run\twall_s\tmax_rss_kb\tlock_wait_s\tstart\texit_status\tattempt\n' >"$times"

for entry in $runs; do
  for h in A B; do
    q=$out/quants/${entry}_$h
    # Done only when the run has its row in quant_times.tsv (written after
    # Salmon exits 0 with its Gibbs draws); a run killed part way is redone.
    if gawk -v r="${entry}_$h" '$1 == r && $6 == "0" { f = 1 } END { exit !f }' "$times" &&
      [ -s "$q/aux_info/bootstrap/bootstraps.gz" ]; then continue; fi
    rm -rf "$q"
    r1=$out/fastq/${entry}_${h}_1.fastq.gz
    r2=$out/fastq/${entry}_${h}_2.fastq.gz
    # Wait for 42_split_reads.sh (files appear only when complete).
    while [ ! -s "$r1" ] || [ ! -s "$r2" ]; do sleep 60; done
    # Up to MAX_ATTEMPTS tries: a memory watchdog or the cap may kill
    # Salmon; a run counts as done only
    # when Salmon exits 0 and wrote quant.sf and its 30 Gibbs draws.
    attempt=0
    while :; do
      attempt=$((attempt + 1))
      rm -rf "$q"
      t_wait=$(date +%s)
      # Memory gate: Salmon peaks at about 11 GB (measured), so wait until
      # MIN_AVAIL_GB of effective available memory (MemAvailable plus the ZFS
      # ARC above its minimum, where ZFS is used) are free.
      while [ "$(eff_avail_gb)" -lt "$min_avail" ]; do
        sleep 60
      done
      # The optional wrapper runs Salmon under a hard memory cap of MEM_GB
      # (exit 137 if exceeded).
      st=0
      ${wrap[@]+"${wrap[@]}"} bash -c '
        t0=$(date +%s)
        echo "$t0" > "$5.t0"
        /usr/bin/time -v -o "$5.time" "$1" quant -i "$2" -l A -1 "$3" -2 "$4" -o "$5" -p "$6" \
          --numGibbsSamples 30 --dumpEq --seqBias --gcBias > "$5.log" 2>&1
        st=$?
        echo $(( $(date +%s) - t0 )) > "$5.wall"
        exit $st
      ' _ "$salmon" "$index" "$r1" "$r2" "$q" "$threads" || st=$?
      cat "$q.time" >>"$q.log" 2>/dev/null || true
      t0=$(cat "$q.t0" 2>/dev/null || echo "$t_wait")
      rss=$(grep 'Maximum resident' "$q.time" 2>/dev/null | gawk '{print $NF}')
      ok=0
      if [ "$st" -eq 0 ] && [ -s "$q/quant.sf" ] && [ -s "$q/aux_info/bootstrap/bootstraps.gz" ] &&
        grep -q '"num_bootstraps": 30' "$q/aux_info/meta_info.json"; then ok=1; fi
      printf '%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$entry" "$h" "$(cat "$q.wall" 2>/dev/null || echo NA)" "${rss:-NA}" \
        "$((t0 - t_wait))" "$(date -d "@$t0" -Is)" "$( [ $ok = 1 ] && echo 0 || echo "${st}_failed")" "$attempt" >>"$times"
      cp "$q.log" "$q.attempt$attempt.log" 2>/dev/null || true
      rm -f "$q.t0" "$q.wall" "$q.time"
      [ $ok = 1 ] && break
      if [ "$attempt" -ge "${MAX_ATTEMPTS:-8}" ]; then
        echo "FAILED $entry $h after $attempt attempts" >&2
        exit 1
      fi
      sleep 120
    done
  done
done
echo "QUANT_DONE $(date -Is)" >>"$times.done"
