#!/usr/bin/env python3
"""Stream the first N paired FASTQ records from an ENA run.

The script stops reading the remote gzip streams after the requested number of
records, so a real-data smoke test does not require downloading the full run.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import urllib.request
from pathlib import Path


def ena_metadata(run: str) -> dict[str, str]:
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run}&result=read_run"
        "&fields=fastq_ftp,fastq_bytes,read_count,base_count,library_layout"
        "&format=tsv"
    )
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one ENA row for {run}, received {len(rows)}")
    return dict(rows[0])


def stream_fastq_subset(url: str, output: Path, records: int) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    request = urllib.request.Request(url, headers={"User-Agent": "bridgeasm-benchmark/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        with gzip.GzipFile(fileobj=response, mode="rb") as source:
            with gzip.open(output, "wb", compresslevel=6) as target:
                for _ in range(records):
                    lines = [source.readline() for _ in range(4)]
                    if not lines[0]:
                        break
                    if any(not line for line in lines[1:]):
                        raise RuntimeError(f"truncated FASTQ record from {url}")
                    if not lines[0].startswith(b"@") or not lines[2].startswith(b"+"):
                        raise RuntimeError(f"malformed FASTQ record from {url}")
                    target.writelines(lines)
                    written += 1
    return written


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", type=int, default=20_000)
    args = parser.parse_args()

    metadata = ena_metadata(args.run)
    ftp_urls = [item for item in metadata.get("fastq_ftp", "").split(";") if item]
    if metadata.get("library_layout") != "PAIRED" or len(ftp_urls) != 2:
        raise RuntimeError(f"{args.run} is not represented by two paired FASTQ files")

    args.output.mkdir(parents=True, exist_ok=True)
    counts: list[int] = []
    files: list[str] = []
    for mate, ftp_url in enumerate(ftp_urls, start=1):
        url = ftp_url if ftp_url.startswith("http") else f"https://{ftp_url}"
        destination = args.output / f"{args.run}_{mate}.fastq.gz"
        counts.append(stream_fastq_subset(url, destination, args.pairs))
        files.append(str(destination))
    if counts[0] != counts[1]:
        raise RuntimeError(f"paired subset count mismatch: {counts}")

    provenance = {
        "run": args.run,
        "requested_pairs": args.pairs,
        "written_pairs": counts[0],
        "files": files,
        "ena_metadata": metadata,
    }
    (args.output / "subset_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
