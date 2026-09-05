#!/usr/bin/env python3
"""Resolve the paired FASTQ source used by an existing benchmark workflow.

The workflow may contain explicit FASTQ URLs or may call ``download_ena_subset.py RUN``.
The latter form is used by current Zymo benchmarks, so this resolver queries ENA for that
run and emits the exact mate URLs. This keeps protein-guide evaluation matched to the read
set that produced the selected GFA artifact.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import urllib.request
from pathlib import Path
from typing import Optional, Sequence

URL_RE = re.compile(r"https?://[^\s'\"<>]+?(?:\.fastq\.gz|\.fq\.gz)")
ENA_RUN_RE = re.compile(r"download_ena_subset\.py\s+([DES]RR\d+)", re.IGNORECASE)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow", type=Path, required=True)
    parser.add_argument("--github-env", type=Path)
    return parser.parse_args(argv)


def ena_fastq_urls(run: str) -> list[str]:
    endpoint = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run}&result=read_run&fields=fastq_ftp,library_layout&format=tsv"
    )
    request = urllib.request.Request(
        endpoint, headers={"User-Agent": "BridgeAsm-protein-pilot/0.1"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one ENA row for {run}, received {len(rows)}")
    row = rows[0]
    if row.get("library_layout") != "PAIRED":
        raise RuntimeError(f"{run} is not paired-end according to ENA")
    ftp_urls = [value for value in row.get("fastq_ftp", "").split(";") if value]
    if len(ftp_urls) != 2:
        raise RuntimeError(f"expected two ENA FASTQ files for {run}, received {len(ftp_urls)}")
    return [value if value.startswith("http") else f"https://{value}" for value in ftp_urls]


def resolve_urls(text: str, workflow: Path) -> tuple[list[str], str | None]:
    urls: list[str] = []
    for match in URL_RE.finditer(text):
        url = match.group(0).rstrip("),]")
        if url not in urls:
            urls.append(url)
    if len(urls) >= 2:
        return urls, None

    runs: list[str] = []
    for match in ENA_RUN_RE.finditer(text):
        run = match.group(1).upper()
        if run not in runs:
            runs.append(run)
    if len(runs) != 1:
        raise SystemExit(
            f"could not resolve one paired FASTQ source in {workflow}; "
            f"explicit_urls={len(urls)} ena_runs={runs}"
        )
    return ena_fastq_urls(runs[0]), runs[0]


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    text = args.workflow.read_text(encoding="utf-8")
    urls, ena_run = resolve_urls(text, args.workflow)

    # Prefer explicit R1/R2 or _1/_2 pairing when auxiliary URLs are present.
    first = next(
        (url for url in urls if re.search(r"(?:R?1|_1)\.f(?:ast)?q\.gz$", url, re.I)),
        urls[0],
    )
    second = next(
        (
            url
            for url in urls
            if url != first and re.search(r"(?:R?2|_2)\.f(?:ast)?q\.gz$", url, re.I)
        ),
        next(url for url in urls if url != first),
    )
    lines = [f"READ1_URL={first}", f"READ2_URL={second}"]
    if ena_run is not None:
        lines.append(f"ENA_RUN={ena_run}")
    if args.github_env:
        with args.github_env.open("a", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
    else:
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
