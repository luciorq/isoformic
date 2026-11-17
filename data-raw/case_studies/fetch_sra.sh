#!/usr/bin/env bash
# Helpers shared by the case study quantification scripts.
#
# Reads come from the NCBI SRA Open Data mirror on AWS S3
# (sra-pub-run-odp), which is much faster from the US than the ENA FTP
# mirror. `fetch_sra RUN DIR` downloads the SRA file of RUN into DIR and
# converts it to FASTQ with fasterq-dump (sra-tools), so that DIR then
# holds RUN.fastq (single end) or RUN_1.fastq and RUN_2.fastq (paired end).
# Callers delete DIR after the quantification to keep the disk use to one
# run at a time.

# fetch_sra RUN DIR THREADS
fetch_sra() {
  local run=$1 dir=$2 threads=${3:-8}
  local url="https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}"
  mkdir -p "$dir"
  if [ ! -s "$dir/$run.sra.ok" ]; then
    curl -sS --retry 5 --retry-all-errors -C - -o "$dir/$run.sra" "$url"
    # The size of the download must match the object size on S3.
    local want have
    want=$(curl -sSI "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" { print $2 }')
    have=$(stat -c %s "$dir/$run.sra")
    if [ "$want" != "$have" ]; then
      echo "Incomplete download of $run ($have of $want bytes)" >&2
      return 1
    fi
    echo "$have" >"$dir/$run.sra.ok"
  fi
  "$fasterq_dump" --split-3 --threads "$threads" --temp "$dir" \
    --outdir "$dir" "$dir/$run.sra"
  rm -f "$dir/$run.sra"
}
