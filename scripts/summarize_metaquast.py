#!/usr/bin/env python3
"""Extract a compact, machine-readable summary from MetaQUAST TSV output."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

DEFAULT_FIELDS = [
    "# contigs",
    "Largest contig",
    "Total length",
    "N50",
    "# misassemblies",
    "# local misassemblies",
    "Genome fraction (%)",
    "Duplication ratio",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp",
    "Total aligned length",
    "NA50",
]


def parse_scalar(value: str) -> Any:
    text = value.strip()
    if text == "-":
        return None
    try:
        if any(character in text for character in ".eE"):
            return float(text)
        return int(text)
    except ValueError:
        return text


def read_report(path: Path, fields: list[str]) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in ["Assembly", *fields] if field not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"missing MetaQUAST columns: {', '.join(missing)}")
        rows = []
        for row in reader:
            summary: dict[str, Any] = {"Assembly": row["Assembly"]}
            for field in fields:
                summary[field] = parse_scalar(row[field])
            rows.append(summary)
        return rows


def write_tsv(rows: list[dict[str, Any]], fields: list[str]) -> None:
    writer = csv.DictWriter(
        __import__("sys").stdout,
        fieldnames=["Assembly", *fields],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path, help="MetaQUAST transposed_report.tsv")
    parser.add_argument("--json", type=Path, help="Optional JSON output path")
    parser.add_argument(
        "--field",
        action="append",
        dest="fields",
        help="Metric to retain; repeat to override the default field set",
    )
    args = parser.parse_args()

    fields = args.fields or DEFAULT_FIELDS
    rows = read_report(args.report, fields)
    write_tsv(rows, fields)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
