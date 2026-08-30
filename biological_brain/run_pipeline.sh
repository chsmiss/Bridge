#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  biological_brain/run_pipeline.sh \
    --gfa graph.gfa \
    (--proteins assembled.faa | --reads reads.fastq[,...]) \
    --output-dir results [options]

Options:
  --threads N             Plass threads (default: 4)
  --mode MODE             conservative|balanced|exploratory (default: balanced)
  --junction-nt N         nucleotide flank per side (default: 450)
  --with-esm              run the optional ESM scorer (requires torch/transformers)
  --esm-model NAME        default: facebook/esm2_t6_8M_UR50D
  --release               use cargo --release

The wrapper never creates graph edges or generated nucleotide sequence.  It scores only
links already present in the input GFA.
EOF
}

GFA=""
PROTEINS=""
READS=""
OUTPUT_DIR=""
THREADS=4
MODE=balanced
JUNCTION_NT=450
WITH_ESM=0
ESM_MODEL=facebook/esm2_t6_8M_UR50D
RELEASE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gfa) GFA=${2:?}; shift 2 ;;
    --proteins) PROTEINS=${2:?}; shift 2 ;;
    --reads) READS=${2:?}; shift 2 ;;
    --output-dir) OUTPUT_DIR=${2:?}; shift 2 ;;
    --threads) THREADS=${2:?}; shift 2 ;;
    --mode) MODE=${2:?}; shift 2 ;;
    --junction-nt) JUNCTION_NT=${2:?}; shift 2 ;;
    --with-esm) WITH_ESM=1; shift ;;
    --esm-model) ESM_MODEL=${2:?}; shift 2 ;;
    --release) RELEASE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "${GFA}" || -z "${OUTPUT_DIR}" ]]; then
  usage >&2
  exit 2
fi
if [[ -n "${PROTEINS}" && -n "${READS}" ]]; then
  echo "choose either --proteins or --reads, not both" >&2
  exit 2
fi
if [[ -z "${PROTEINS}" && -z "${READS}" ]]; then
  echo "one of --proteins or --reads is required" >&2
  exit 2
fi
if [[ ! -s "${GFA}" ]]; then
  echo "missing or empty GFA: ${GFA}" >&2
  exit 2
fi
if [[ ! "${MODE}" =~ ^(conservative|balanced|exploratory)$ ]]; then
  echo "invalid mode: ${MODE}" >&2
  exit 2
fi

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "${OUTPUT_DIR}"

if [[ -n "${READS}" ]]; then
  command -v plass >/dev/null 2>&1 || {
    echo "Plass is required when --reads is used" >&2
    exit 2
  }
  combined="${OUTPUT_DIR}/protein_input_reads.fastq"
  : > "${combined}"
  IFS=',' read -r -a read_files <<< "${READS}"
  for read_file in "${read_files[@]}"; do
    if [[ ! -s "${read_file}" ]]; then
      echo "missing or empty reads file: ${read_file}" >&2
      exit 2
    fi
    cat "${read_file}" >> "${combined}"
  done
  PROTEINS="${OUTPUT_DIR}/plass_proteins.faa"
  mkdir -p "${OUTPUT_DIR}/plass_tmp"
  plass assemble "${combined}" "${PROTEINS}" "${OUTPUT_DIR}/plass_tmp" --threads "${THREADS}"
fi

if [[ ! -s "${PROTEINS}" ]]; then
  echo "missing or empty protein assembly: ${PROTEINS}" >&2
  exit 2
fi

python3 "${ROOT}/biological_brain/protein_bridge_evidence.py" \
  --gfa "${GFA}" \
  --proteins "${PROTEINS}" \
  --output "${OUTPUT_DIR}/protein_evidence.tsv" \
  --junction-nt "${JUNCTION_NT}"

python3 "${ROOT}/biological_brain/kmer_junction_lm.py" \
  --gfa "${GFA}" \
  --output "${OUTPUT_DIR}/dna_lm.tsv"

ESM_ARGS=()
if [[ "${WITH_ESM}" == "1" ]]; then
  python3 "${ROOT}/biological_brain/esm_breakpoint_score.py" \
    --candidates "${OUTPUT_DIR}/protein_evidence.tsv" \
    --output "${OUTPUT_DIR}/esm_scores.tsv" \
    --model "${ESM_MODEL}"
  ESM_ARGS=(--esm-scores "${OUTPUT_DIR}/esm_scores.tsv")
fi

CARGO_ARGS=(run --quiet --bin bridgeasm-evidence-path)
if [[ "${RELEASE}" == "1" ]]; then
  CARGO_ARGS=(run --release --quiet --bin bridgeasm-evidence-path)
fi

(
  cd "${ROOT}"
  cargo "${CARGO_ARGS[@]:1}" -- \
    --gfa "${GFA}" \
    --edge-evidence "${OUTPUT_DIR}/protein_evidence.tsv" \
    --dna-lm-scores "${OUTPUT_DIR}/dna_lm.tsv" \
    "${ESM_ARGS[@]}" \
    --output "${OUTPUT_DIR}/evidence_contigs.fasta" \
    --report "${OUTPUT_DIR}/evidence_path_report.tsv" \
    --mode "${MODE}"
)

python3 "${ROOT}/biological_brain/assembly_stats.py" \
  --input "${OUTPUT_DIR}/evidence_contigs.fasta" \
  --format fasta \
  --minimum 200 \
  --evidence "${OUTPUT_DIR}/protein_evidence.tsv" \
  --label biological_brain \
  --output "${OUTPUT_DIR}/assembly_stats.json"

cat "${OUTPUT_DIR}/assembly_stats.json"
