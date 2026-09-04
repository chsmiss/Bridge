#!/usr/bin/env python3
"""Outer-join BridgeBin Biological Brain feature TSVs by contig.

Foundation-model inference intentionally lives in separate adapters (DNA LM, GENERanno,
ESM-C, protein repertoire). The pair head, however, consumes one row per contig. This
utility performs a strict, auditable join without silently overwriting conflicting values.

Examples:
  python scripts/bridgebin_join_features.py \
    --input bio.tsv \
    --input generanno.architecture.tsv \
    --input protein_repertoire.tsv \
    --output biobrain_features.tsv

Rules:
- ``contig`` / ``contig_id`` / ``sequence`` identifies the row.
- all non-ID columns are preserved in first-seen order;
- empty/./NA values do not overwrite a populated value;
- two different populated values for the same contig+column are an error by default;
- ``--prefer-last`` explicitly allows the later input to win.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

MISSING = {"", ".", "NA", "na", "NaN", "nan"}
ID_COLUMNS = ("contig", "contig_id", "sequence")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--prefer-last",
        action="store_true",
        help="explicitly let later inputs replace conflicting non-missing values",
    )
    return parser.parse_args(argv)


def normalize(value: Optional[str]) -> str:
    if value is None:
        return "."
    value = value.strip()
    return "." if value in MISSING else value


def read_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        id_column = next((column for column in ID_COLUMNS if column in reader.fieldnames), None)
        if id_column is None:
            raise ValueError(f"{path}: needs one of ID columns {ID_COLUMNS}")
        payload_columns = [column for column in reader.fieldnames if column != id_column]
        rows: List[Dict[str, str]] = []
        seen = set()
        for line_number, raw in enumerate(reader, start=2):
            contig = normalize(raw.get(id_column))
            if contig == ".":
                continue
            if contig in seen:
                raise ValueError(f"{path}:{line_number}: duplicate contig {contig!r}")
            seen.add(contig)
            row = {"contig": contig}
            for column in payload_columns:
                row[column] = normalize(raw.get(column))
            rows.append(row)
        return payload_columns, rows


def merge_value(
    contig: str,
    column: str,
    old: str,
    new: str,
    source: Path,
    prefer_last: bool,
) -> str:
    if old == ".":
        return new
    if new == "." or new == old:
        return old
    if prefer_last:
        return new
    raise ValueError(
        f"conflicting feature for contig={contig!r} column={column!r}: "
        f"existing={old!r}, new={new!r} from {source}; use --prefer-last only if intentional"
    )


def join(inputs: Sequence[Path], prefer_last: bool) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    columns: List[str] = []
    data: Dict[str, Dict[str, str]] = {}
    for path in inputs:
        table_columns, table_rows = read_table(path)
        for column in table_columns:
            if column not in columns:
                columns.append(column)
        for row in table_rows:
            contig = row["contig"]
            target = data.setdefault(contig, {})
            for column in table_columns:
                new = row.get(column, ".")
                old = target.get(column, ".")
                target[column] = merge_value(
                    contig, column, old, new, path, prefer_last
                )
    return columns, data


def write_output(path: Path, columns: Sequence[str], data: Dict[str, Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", *columns])
        for contig in sorted(data):
            row = data[contig]
            writer.writerow([contig, *[row.get(column, ".") for column in columns]])


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    columns, data = join(args.input, args.prefer_last)
    write_output(args.output, columns, data)
    populated = {
        column: sum(data[contig].get(column, ".") != "." for contig in data)
        for column in columns
    }
    summary = " ".join(f"{column}={count}" for column, count in populated.items())
    print(
        f"bridgebin-join-features: inputs={len(args.input)} contigs={len(data)} "
        f"columns={len(columns)} {summary}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
