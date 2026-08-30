#!/usr/bin/env python3
"""Export exact protein-guided candidate junction windows for external scorers.

The exporter reconstructs each eligible ``bridgeasm-proteinguide`` candidate
from the immutable backbone GFA, the PenguiN guide FASTA, the original PAF, and
the edge report.  It writes joined DNA windows without inventing sequence:
positive-gap bases come only from the guide, and overlap joins reuse only the
validated overlap reported by the Rust candidate generator.

The resulting FASTA/manifest are suitable for:

* DNA-language-model joined-vs-alternative likelihood scoring, and
* Prodigal translation followed by ESM coding-continuity scoring.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Anchor:
    query: str
    target: str
    reverse: bool
    full_start: int
    full_end: int
    identity: float
    query_fraction: float
    mapq: int
    alignment: int

    @property
    def score(self) -> float:
        return (
            self.identity * 1_000_000.0
            + self.query_fraction * 100_000.0
            + self.mapq * 1_000.0
            + self.alignment
        )


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []

    def commit() -> None:
        nonlocal name, chunks
        if name is None:
            return
        if name in records:
            raise ValueError(f"duplicate FASTA name: {name}")
        records[name] = "".join(chunks).upper()
        name = None
        chunks = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError("empty FASTA header")
            else:
                if name is None:
                    raise ValueError("sequence before first FASTA header")
                chunks.append(line)
    commit()
    return records


def read_gfa(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.startswith("S\t"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] == "*":
                raise ValueError(f"invalid GFA segment on line {line_number}")
            if fields[1] in sequences:
                raise ValueError(f"duplicate GFA segment name: {fields[1]}")
            sequences[fields[1]] = fields[2].upper()
    return sequences


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return sequence.translate(table)[::-1]


def normalize_identity(value: float) -> float:
    return value / 100.0 if value > 1.0 else value


def project_interval(
    query_len: int,
    query_start: int,
    query_end: int,
    target_len: int,
    target_start: int,
    target_end: int,
    reverse: bool,
) -> tuple[int, int]:
    if reverse:
        query_left = query_len - query_end
        query_right = query_start
    else:
        query_left = query_start
        query_right = query_len - query_end
    start = max(0, target_start - query_left)
    end = min(target_len, target_end + query_right)
    return start, end


def read_paf(path: Path, min_mapq: int) -> dict[tuple[str, str, bool], list[Anchor]]:
    anchors: dict[tuple[str, str, bool], list[Anchor]] = {}
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 columns")
            query, target = fields[0], fields[5]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            reverse = fields[4] == "-"
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matches = int(fields[9])
            alignment = int(fields[10])
            mapq = int(fields[11])
            if mapq < min_mapq or alignment <= 0 or query_len <= 0:
                continue
            full_start, full_end = project_interval(
                query_len,
                query_start,
                query_end,
                target_len,
                target_start,
                target_end,
                reverse,
            )
            anchor = Anchor(
                query=query,
                target=target,
                reverse=reverse,
                full_start=full_start,
                full_end=full_end,
                identity=matches / alignment,
                query_fraction=(query_end - query_start) / query_len,
                mapq=mapq,
                alignment=alignment,
            )
            anchors.setdefault((query, target, reverse), []).append(anchor)
    for key, values in anchors.items():
        values.sort(key=lambda anchor: -anchor.score)
        anchors[key] = values
    return anchors


def parse_oriented(label: str) -> tuple[str, bool]:
    if label.endswith("+"):
        return label[:-1], False
    if label.endswith("-"):
        return label[:-1], True
    raise ValueError(f"oriented endpoint lacks +/- suffix: {label}")


def choose_anchor_pair(
    source_values: list[Anchor],
    target_values: list[Anchor],
    expected_gap: int,
) -> tuple[Anchor, Anchor, int] | None:
    best: tuple[tuple[int, float], Anchor, Anchor, int] | None = None
    for source in source_values:
        for target in target_values:
            reconstructed_gap = target.full_start - source.full_end
            rank = (
                abs(reconstructed_gap - expected_gap),
                -(source.score + target.score),
            )
            if best is None or rank < best[0]:
                best = (rank, source, target, reconstructed_gap)
    if best is None:
        return None
    return best[1], best[2], best[3]


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfa", type=Path, required=True)
    parser.add_argument("--guide", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--edge-report", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--left-context", type=int, default=1024)
    parser.add_argument("--right-context", type=int, default=1024)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--max-gap-disagreement", type=int, default=30)
    parser.add_argument("--selected-only", action="store_true")
    args = parser.parse_args()

    if args.left_context <= 0 or args.right_context <= 0:
        parser.error("context sizes must be positive")
    if args.max_gap_disagreement < 0:
        parser.error("--max-gap-disagreement must be non-negative")

    segments = read_gfa(args.gfa)
    guides = read_fasta(args.guide)
    anchors = read_paf(args.paf, args.min_mapq)

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "edge_id",
        "source",
        "target",
        "guide",
        "projected_gap",
        "reconstructed_gap",
        "gap_disagreement",
        "overlap",
        "bridge_bases",
        "left_context_bases",
        "right_context_bases",
        "window_bases",
        "junction_offset",
        "source_identity",
        "target_identity",
        "source_query_fraction",
        "target_query_fraction",
        "source_mapq",
        "target_mapq",
    ]

    exported = 0
    skipped = 0
    with args.edge_report.open(newline="") as edge_handle, args.output_fasta.open(
        "w"
    ) as fasta_handle, args.manifest.open("w", newline="") as manifest_handle:
        reader = csv.DictReader(edge_handle, delimiter="\t")
        required = {
            "source",
            "target",
            "guide",
            "eligible",
            "selected",
            "projected_gap",
            "overlap",
            "guide_bases",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"edge report is missing columns: {', '.join(sorted(missing))}")
        writer = csv.DictWriter(manifest_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row_index, row in enumerate(reader, 1):
            if row["eligible"].lower() != "true":
                continue
            if args.selected_only and row["selected"].lower() != "true":
                continue
            source_name, source_reverse = parse_oriented(row["source"])
            target_name, target_reverse = parse_oriented(row["target"])
            guide_name = row["guide"]
            projected_gap = int(row["projected_gap"])
            overlap = int(row["overlap"])
            if source_name not in segments or target_name not in segments or guide_name not in guides:
                skipped += 1
                continue
            source_values = anchors.get((source_name, guide_name, source_reverse), [])
            target_values = anchors.get((target_name, guide_name, target_reverse), [])
            pair = choose_anchor_pair(source_values, target_values, projected_gap)
            if pair is None:
                skipped += 1
                continue
            source_anchor, target_anchor, reconstructed_gap = pair
            disagreement = abs(reconstructed_gap - projected_gap)
            if disagreement > args.max_gap_disagreement:
                skipped += 1
                continue

            source_sequence = segments[source_name]
            target_sequence = segments[target_name]
            if source_reverse:
                source_sequence = reverse_complement(source_sequence)
            if target_reverse:
                target_sequence = reverse_complement(target_sequence)
            if overlap < 0 or overlap > len(target_sequence):
                skipped += 1
                continue

            bridge = ""
            if projected_gap > 0:
                bridge_start = max(0, source_anchor.full_end)
                bridge_end = min(len(guides[guide_name]), target_anchor.full_start)
                if bridge_end < bridge_start:
                    skipped += 1
                    continue
                bridge = guides[guide_name][bridge_start:bridge_end]
                expected_guide_bases = int(row["guide_bases"])
                if abs(len(bridge) - expected_guide_bases) > args.max_gap_disagreement:
                    skipped += 1
                    continue

            left = source_sequence[-args.left_context :]
            right = target_sequence[overlap : overlap + args.right_context]
            window = left + bridge + right
            if not window:
                skipped += 1
                continue
            edge_id = f"junction_{row_index:06d}"
            junction_offset = len(left) + len(bridge)
            fasta_handle.write(
                f">{edge_id} source={row['source']} target={row['target']} "
                f"guide={guide_name} gap={projected_gap} overlap={overlap} "
                f"junction_offset={junction_offset}\n{wrap(window)}\n"
            )
            writer.writerow(
                {
                    "edge_id": edge_id,
                    "source": row["source"],
                    "target": row["target"],
                    "guide": guide_name,
                    "projected_gap": projected_gap,
                    "reconstructed_gap": reconstructed_gap,
                    "gap_disagreement": disagreement,
                    "overlap": overlap,
                    "bridge_bases": len(bridge),
                    "left_context_bases": len(left),
                    "right_context_bases": len(right),
                    "window_bases": len(window),
                    "junction_offset": junction_offset,
                    "source_identity": f"{source_anchor.identity:.6f}",
                    "target_identity": f"{target_anchor.identity:.6f}",
                    "source_query_fraction": f"{source_anchor.query_fraction:.6f}",
                    "target_query_fraction": f"{target_anchor.query_fraction:.6f}",
                    "source_mapq": source_anchor.mapq,
                    "target_mapq": target_anchor.mapq,
                }
            )
            exported += 1

    print(f"exported={exported}\tskipped={skipped}")


if __name__ == "__main__":
    main()
