#!/usr/bin/env python3
"""Download the complete paired FASTQ files for one ENA run with provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import shutil
import urllib.request
from pathlib import Path


def ena_metadata(run: str) -> dict[str, str]:
    fields = "fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count,library_layout"
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run}&result=read_run&fields={fields}&format=tsv"
    )
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one ENA row for {run}, received {len(rows)}")
    return dict(rows[0])


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "bridgeasm-stage34/0.1"})
    with urllib.request.urlopen(request, timeout=180) as response:
        with destination.open("wb") as target:
            shutil.copyfileobj(response, target, length=8 << 20)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("output", type=Path)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()

    metadata = ena_metadata(args.run)
    ftp_urls = [item for item in metadata.get("fastq_ftp", "").split(";") if item]
    expected_bytes = [int(item) for item in metadata.get("fastq_bytes", "").split(";") if item]
    expected_md5 = [item for item in metadata.get("fastq_md5", "").split(";") if item]
    if metadata.get("library_layout") != "PAIRED" or len(ftp_urls) != 2:
        raise RuntimeError(f"{args.run} is not represented by two paired FASTQ files")

    args.output.mkdir(parents=True, exist_ok=True)
    files: list[str] = []
    if not args.metadata_only:
        for mate, ftp_url in enumerate(ftp_urls, start=1):
            url = ftp_url if ftp_url.startswith("http") else f"https://{ftp_url}"
            destination = args.output / f"{args.run}_{mate}.fastq.gz"
            download(url, destination)
            files.append(str(destination))
            if len(expected_bytes) == 2 and destination.stat().st_size != expected_bytes[mate - 1]:
                raise RuntimeError(
                    f"mate {mate} byte mismatch: {destination.stat().st_size} != {expected_bytes[mate - 1]}"
                )
            if len(expected_md5) == 2:
                observed = md5_file(destination)
                if observed != expected_md5[mate - 1]:
                    raise RuntimeError(f"mate {mate} md5 mismatch: {observed} != {expected_md5[mate - 1]}")

    provenance = {
        "run": args.run,
        "complete_run": not args.metadata_only,
        "files": files,
        "expected_compressed_bytes": expected_bytes,
        "expected_compressed_bytes_total": sum(expected_bytes),
        "ena_metadata": metadata,
    }
    path = args.output / "full_provenance.json"
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps(provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
