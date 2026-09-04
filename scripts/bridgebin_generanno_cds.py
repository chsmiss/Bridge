#!/usr/bin/env python3
"""Run GENERanno prokaryotic CDS annotation and extract BridgeBin-ready ORFs.

The public ``GenerTeam/GENERanno-prokaryote-0.5b-cds-annotator`` is a two-head,
single-nucleotide token classifier: the heads represent positive and negative strands and
predict CDS vs non-coding state for each base. This adapter overlaps long contigs into
8-kb windows, averages CDS probabilities across overlaps, then calls biologically valid
open reading frames that are supported by the corresponding strand probability track.

Outputs from ``--output-prefix PREFIX``:
  PREFIX.cds.tsv          contig/gene coordinates and GENERanno support
  PREFIX.proteins.faa     translated ORFs for ESM-C
  PREFIX.protein_map.tsv  protein_id -> contig mapping for bridgebin_esmc_embed.py
  PREFIX.architecture.tsv contig-level coding architecture vector

GENERanno predicts coding state rather than protein family/function. Protein repertoire is
therefore intentionally delegated to downstream ESM-C/protein-family processing instead
of being invented from CDS labels.

Optional dependencies: torch, transformers.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class Orf:
    contig: str
    start: int
    end: int
    strand: str
    frame: int
    cds_probability: float
    nt: str
    protein: str


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument(
        "--model",
        default="GenerTeam/GENERanno-prokaryote-0.5b-cds-annotator",
    )
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--window-bp", type=int, default=8192)
    parser.add_argument("--overlap-bp", type=int, default=1024)
    parser.add_argument("--min-orf-bp", type=int, default=90)
    parser.add_argument("--max-orf-bp", type=int, default=12000)
    parser.add_argument("--min-cds-probability", type=float, default=0.60)
    parser.add_argument("--min-start-probability", type=float, default=0.45)
    parser.add_argument("--allow-edge-truncated", action="store_true")
    parser.add_argument("--bf16", action="store_true")
    return parser.parse_args(argv)


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"{path}: sequence before FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def translate(sequence: str) -> str:
    amino = []
    for offset in range(0, len(sequence) - 2, 3):
        amino.append(CODON_TABLE.get(sequence[offset : offset + 3], "X"))
    protein = "".join(amino)
    return protein[:-1] if protein.endswith("*") else protein


def load_runtime(model_name: str, requested_device: str, bf16: bool):
    try:
        import torch
        from transformers import AutoModelForTokenClassification, AutoTokenizer
    except ImportError as error:  # pragma: no cover - optional model dependency
        raise RuntimeError("install optional dependencies: pip install torch transformers") from error

    if requested_device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = requested_device
    dtype = torch.bfloat16 if bf16 and device == "cuda" else None
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    kwargs = {"trust_remote_code": True}
    if dtype is not None:
        kwargs["torch_dtype"] = dtype
    model = AutoModelForTokenClassification.from_pretrained(model_name, **kwargs)
    model.to(device)
    model.eval()

    id2label = {int(key): str(value).upper() for key, value in model.config.id2label.items()}
    cds_ids = [index for index, label in id2label.items() if label == "CDS"]
    if len(cds_ids) != 1:
        raise RuntimeError(f"expected exactly one CDS label, got {id2label}")
    num_heads = int(getattr(model, "num_prediction_heads", 1))
    if num_heads < 2:
        raise RuntimeError(
            f"expected positive/negative strand prediction heads, model reports {num_heads}"
        )
    return torch, tokenizer, model, device, cds_ids[0], num_heads


def window_ranges(length: int, window_bp: int, overlap_bp: int):
    if length <= window_bp:
        yield 0, length
        return
    stride = window_bp - overlap_bp
    start = 0
    while start < length:
        end = min(length, start + window_bp)
        yield start, end
        if end == length:
            break
        start += stride


def predict_window(sequence: str, runtime) -> Tuple[List[float], List[float]]:
    torch, tokenizer, model, device, cds_id, num_heads = runtime
    inputs = tokenizer(sequence, add_special_tokens=False, return_tensors="pt")
    if int(inputs["input_ids"].shape[1]) != len(sequence):
        raise RuntimeError(
            "GENERanno CDS tokenizer is expected to produce one token per nucleotide; "
            f"got {inputs['input_ids'].shape[1]} tokens for {len(sequence)} bases"
        )
    inputs = {key: value.to(device) for key, value in inputs.items()}
    with torch.inference_mode():
        logits = model(**inputs).logits[0]
        probs = torch.softmax(logits.float(), dim=-1)[:, cds_id].cpu()
    expected = len(sequence) * num_heads
    if probs.numel() != expected:
        raise RuntimeError(
            f"unexpected GENERanno output length {probs.numel()}, expected {expected} "
            f"({len(sequence)} bases x {num_heads} heads)"
        )
    head = probs.reshape(num_heads, len(sequence))
    return head[0].tolist(), head[1].tolist()


def predict_contig(sequence: str, runtime, window_bp: int, overlap_bp: int):
    positive_sum = [0.0] * len(sequence)
    negative_sum = [0.0] * len(sequence)
    weights = [0] * len(sequence)
    for start, end in window_ranges(len(sequence), window_bp, overlap_bp):
        positive, negative = predict_window(sequence[start:end], runtime)
        for local, position in enumerate(range(start, end)):
            positive_sum[position] += positive[local]
            negative_sum[position] += negative[local]
            weights[position] += 1
    positive = [value / max(1, weight) for value, weight in zip(positive_sum, weights)]
    negative = [value / max(1, weight) for value, weight in zip(negative_sum, weights)]
    return positive, negative


def prefix(values: Sequence[float]) -> List[float]:
    out = [0.0]
    total = 0.0
    for value in values:
        total += value
        out.append(total)
    return out


def mean_range(sums: Sequence[float], start: int, end: int) -> float:
    if end <= start:
        return 0.0
    return (sums[end] - sums[start]) / (end - start)


def scan_orfs(
    contig: str,
    oriented_sequence: str,
    probabilities: Sequence[float],
    strand: str,
    original_length: int,
    min_orf_bp: int,
    max_orf_bp: int,
    min_cds_probability: float,
    min_start_probability: float,
    allow_edge_truncated: bool,
) -> List[Orf]:
    sums = prefix(probabilities)
    candidates: List[Orf] = []
    for frame in range(3):
        active_start: Optional[int] = None
        offset = frame
        while offset + 3 <= len(oriented_sequence):
            codon = oriented_sequence[offset : offset + 3]
            if active_start is None and codon in START_CODONS:
                local_start = mean_range(sums, offset, min(len(probabilities), offset + 9))
                if local_start >= min_start_probability:
                    active_start = offset
            if active_start is not None and codon in STOP_CODONS:
                end = offset + 3
                length = end - active_start
                if min_orf_bp <= length <= max_orf_bp:
                    score = mean_range(sums, active_start, end)
                    if score >= min_cds_probability:
                        nt = oriented_sequence[active_start:end]
                        protein = translate(nt)
                        if strand == "+":
                            original_start, original_end = active_start, end
                        else:
                            original_start = original_length - end
                            original_end = original_length - active_start
                        candidates.append(
                            Orf(
                                contig=contig,
                                start=original_start,
                                end=original_end,
                                strand=strand,
                                frame=frame,
                                cds_probability=score,
                                nt=nt,
                                protein=protein,
                            )
                        )
                active_start = None
            elif active_start is not None and offset + 3 - active_start > max_orf_bp:
                active_start = None
            offset += 3

        if allow_edge_truncated and active_start is not None:
            end = len(oriented_sequence) - ((len(oriented_sequence) - active_start) % 3)
            length = end - active_start
            if min_orf_bp <= length <= max_orf_bp:
                score = mean_range(sums, active_start, end)
                if score >= min_cds_probability:
                    nt = oriented_sequence[active_start:end]
                    protein = translate(nt)
                    if strand == "+":
                        original_start, original_end = active_start, end
                    else:
                        original_start = original_length - end
                        original_end = original_length - active_start
                    candidates.append(
                        Orf(
                            contig=contig,
                            start=original_start,
                            end=original_end,
                            strand=strand,
                            frame=frame,
                            cds_probability=score,
                            nt=nt,
                            protein=protein,
                        )
                    )
    return candidates


def nonredundant(orfs: Sequence[Orf]) -> List[Orf]:
    # Multiple overlapping starts in one CDS are common. Keep the strongest supported
    # longest ORFs greedily; overlapping opposite-strand genes remain possible because
    # competition is restricted to the same strand.
    ranked = sorted(
        orfs,
        key=lambda orf: (-orf.cds_probability, -(orf.end - orf.start), orf.start, orf.strand),
    )
    kept: List[Orf] = []
    for candidate in ranked:
        conflict = False
        for existing in kept:
            if existing.strand != candidate.strand:
                continue
            overlap = max(0, min(existing.end, candidate.end) - max(existing.start, candidate.start))
            shorter = min(existing.end - existing.start, candidate.end - candidate.start)
            if shorter > 0 and overlap / shorter >= 0.60:
                conflict = True
                break
        if not conflict:
            kept.append(candidate)
    kept.sort(key=lambda orf: (orf.start, orf.end, orf.strand))
    return kept


def quantiles(values: Sequence[int]) -> Tuple[float, float, float]:
    if not values:
        return 0.0, 0.0, 0.0
    ordered = sorted(values)
    def at(fraction: float) -> float:
        position = fraction * (len(ordered) - 1)
        lo = int(math.floor(position))
        hi = int(math.ceil(position))
        if lo == hi:
            return float(ordered[lo])
        mix = position - lo
        return ordered[lo] * (1.0 - mix) + ordered[hi] * mix
    return at(0.25), at(0.50), at(0.75)


def architecture_vector(
    sequence: str, positive: Sequence[float], negative: Sequence[float], orfs: Sequence[Orf]
) -> Tuple[List[float], float]:
    length = max(1, len(sequence))
    positive_fraction = sum(value >= 0.5 for value in positive) / length
    negative_fraction = sum(value >= 0.5 for value in negative) / length
    coding_union = sum(max(p, n) >= 0.5 for p, n in zip(positive, negative)) / length
    lengths = [orf.end - orf.start for orf in orfs]
    q25, median, q75 = quantiles(lengths)
    plus = sum(orf.strand == "+" for orf in orfs)
    minus = sum(orf.strand == "-" for orf in orfs)
    strand_balance = (plus - minus) / max(1, plus + minus)
    mean_orf = sum(lengths) / max(1, len(lengths))
    mean_support = sum(orf.cds_probability for orf in orfs) / max(1, len(orfs))
    vector = [
        positive_fraction,
        negative_fraction,
        coding_union,
        len(orfs) * 1000.0 / length,
        mean_orf / 3000.0,
        q25 / 3000.0,
        median / 3000.0,
        q75 / 3000.0,
        strand_balance,
        mean_support,
    ]
    confidence = min(1.0, (1.0 - math.exp(-len(orfs) / 8.0)) * max(0.25, mean_support))
    return vector, confidence


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.window_bp < 128 or not 0 <= args.overlap_bp < args.window_bp:
        raise SystemExit("require --window-bp >= 128 and 0 <= --overlap-bp < --window-bp")
    if args.min_orf_bp < 30 or args.max_orf_bp < args.min_orf_bp:
        raise SystemExit("invalid ORF length bounds")
    if not 0.0 <= args.min_cds_probability <= 1.0:
        raise SystemExit("--min-cds-probability must be in [0,1]")

    runtime = load_runtime(args.model, args.device, args.bf16)
    prefix_path = args.output_prefix
    prefix_path.parent.mkdir(parents=True, exist_ok=True)
    cds_path = Path(str(prefix_path) + ".cds.tsv")
    protein_path = Path(str(prefix_path) + ".proteins.faa")
    map_path = Path(str(prefix_path) + ".protein_map.tsv")
    architecture_path = Path(str(prefix_path) + ".architecture.tsv")

    total_orfs = total_contigs = 0
    with cds_path.open("w", encoding="utf-8", newline="") as cds_handle, protein_path.open(
        "w", encoding="utf-8"
    ) as protein_handle, map_path.open("w", encoding="utf-8", newline="") as map_handle, architecture_path.open(
        "w", encoding="utf-8", newline=""
    ) as architecture_handle:
        cds_writer = csv.writer(cds_handle, delimiter="\t", lineterminator="\n")
        cds_writer.writerow(
            ["contig", "gene_id", "start", "end", "strand", "frame", "cds_probability", "nt_length"]
        )
        map_writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        map_writer.writerow(["protein_id", "contig"])
        architecture_writer = csv.writer(
            architecture_handle, delimiter="\t", lineterminator="\n"
        )
        architecture_writer.writerow(
            ["contig", "gene_architecture", "architecture_confidence", "model"]
        )

        for contig, sequence in read_fasta(args.contigs):
            if len(sequence) < args.min_orf_bp:
                continue
            positive, negative = predict_contig(
                sequence, runtime, args.window_bp, args.overlap_bp
            )
            plus = scan_orfs(
                contig,
                sequence,
                positive,
                "+",
                len(sequence),
                args.min_orf_bp,
                args.max_orf_bp,
                args.min_cds_probability,
                args.min_start_probability,
                args.allow_edge_truncated,
            )
            rc = reverse_complement(sequence)
            minus = scan_orfs(
                contig,
                rc,
                list(reversed(negative)),
                "-",
                len(sequence),
                args.min_orf_bp,
                args.max_orf_bp,
                args.min_cds_probability,
                args.min_start_probability,
                args.allow_edge_truncated,
            )
            orfs = nonredundant(plus + minus)
            architecture, architecture_confidence = architecture_vector(
                sequence, positive, negative, orfs
            )
            architecture_writer.writerow(
                [
                    contig,
                    ",".join(f"{value:.7g}" for value in architecture),
                    f"{architecture_confidence:.6f}",
                    args.model,
                ]
            )
            total_contigs += 1
            for serial, orf in enumerate(orfs, start=1):
                gene_id = f"{contig}|generanno_cds_{serial:05d}"
                cds_writer.writerow(
                    [
                        contig,
                        gene_id,
                        orf.start + 1,
                        orf.end,
                        orf.strand,
                        orf.frame,
                        f"{orf.cds_probability:.6f}",
                        orf.end - orf.start,
                    ]
                )
                protein_handle.write(f">{gene_id}\n")
                for position in range(0, len(orf.protein), 80):
                    protein_handle.write(orf.protein[position : position + 80] + "\n")
                map_writer.writerow([gene_id, contig])
                total_orfs += 1

    print(
        f"bridgebin-generanno: contigs={total_contigs} orfs={total_orfs} model={args.model} "
        f"cds={cds_path} proteins={protein_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
