#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=${1:-"$ROOT/benchmark_synthetic"}
rm -rf "$WORK"
mkdir -p "$WORK"

python3 "$ROOT/scripts/simulate_metagenome.py" --output "$WORK/data" --seed 43
cargo run --release --manifest-path "$ROOT/Cargo.toml" -- assemble \
  -1 "$WORK/data/reads_R1.fastq.gz" \
  -2 "$WORK/data/reads_R2.fastq.gz" \
  -o "$WORK/bridgeasm" \
  -k 31 --min-count 2 --mercy-max-kmers 16 --mercy-min-support 1 \
  --min-read-support 2 --min-primary-support 5 --primary-dominance 0.75 \
  --min-contig-length 200 --threads 4
python3 "$ROOT/scripts/evaluate_synthetic.py" \
  --major "$WORK/data/major.fasta" \
  --minor "$WORK/data/minor.fasta" \
  --variants "$WORK/data/minor_variants.tsv" \
  --contigs "$WORK/bridgeasm/primary_contigs.fasta" \
  --additional "$WORK/bridgeasm/variants.fasta" \
  --additional "$WORK/bridgeasm/haplotigs.fasta" | tee "$WORK/evaluation.tsv"
