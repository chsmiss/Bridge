#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 R1.fastq.gz R2.fastq.gz OUTDIR [REFERENCES.fasta] [MAX_PAIRS]" >&2
  exit 2
fi

R1=$1
R2=$2
OUT=$3
REFS=${4:-}
MAX_PAIRS=${5:-}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"

cargo build --release --manifest-path "$ROOT/Cargo.toml"
CMD=("$ROOT/target/release/bridgeasm" assemble -1 "$R1" -2 "$R2" -o "$OUT/bridgeasm" -k 31 \
     --min-count 2 --mercy-max-kmers 16 --mercy-min-support 1 \
     --min-read-support 2 --min-contig-length 200 --threads "${THREADS:-8}")
if [[ -n "$MAX_PAIRS" ]]; then
  CMD+=(--max-pairs "$MAX_PAIRS")
fi
/usr/bin/time -v -o "$OUT/bridgeasm.time.txt" "${CMD[@]}"

if command -v megahit >/dev/null 2>&1; then
  MEGAHIT_CMD=(megahit -1 "$R1" -2 "$R2" -o "$OUT/megahit" -t "${THREADS:-8}")
  if [[ -n "$MAX_PAIRS" ]]; then
    echo "MEGAHIT subset comparison requires physically subset FASTQs; skipping." >&2
  else
    /usr/bin/time -v -o "$OUT/megahit.time.txt" "${MEGAHIT_CMD[@]}"
  fi
fi

if [[ -n "$REFS" ]] && command -v metaquast.py >/dev/null 2>&1; then
  ASSEMBLIES=("$OUT/bridgeasm/primary_contigs.fasta")
  [[ -s "$OUT/megahit/final.contigs.fa" ]] && ASSEMBLIES+=("$OUT/megahit/final.contigs.fa")
  metaquast.py "${ASSEMBLIES[@]}" -r "$REFS" -o "$OUT/metaquast" --threads "${THREADS:-8}"
fi
