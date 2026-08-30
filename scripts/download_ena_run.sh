#!/usr/bin/env bash
set -euo pipefail

RUN=${1:-ERR2935805}
OUT=${2:-data/${RUN}}
mkdir -p "$OUT"
REPORT=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=fastq_ftp,fastq_bytes,read_count,base_count&format=tsv")
printf '%s\n' "$REPORT" | tee "$OUT/ena_report.tsv"
URLS=$(printf '%s\n' "$REPORT" | awk -F '\t' 'NR==2 {print $2}')
if [[ -z "$URLS" ]]; then
  echo "No FASTQ URLs found for ${RUN}" >&2
  exit 1
fi
IFS=';' read -r -a FILES <<< "$URLS"
for item in "${FILES[@]}"; do
  url="https://${item}"
  target="$OUT/$(basename "$item")"
  if [[ -s "$target" ]]; then
    echo "exists: $target"
  else
    echo "downloading: $url"
    curl -fL --retry 5 --continue-at - "$url" -o "$target"
  fi
done
