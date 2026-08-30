#!/usr/bin/env python3
"""Use a small ESM masked language model to score protein continuity at graph junctions.

The input is the TSV produced by protein_bridge_evidence.py.  The scorer compares local
pseudo-log-likelihood with the two peptide halves joined versus scored independently.
A positive esm_delta is weak support for a biologically coherent continuation.  It is a
reranking signal only; bridgeasm-evidence-path ignores it unless nucleotide/protein
support already makes an existing GFA edge eligible.

Dependencies are intentionally optional:
    pip install torch transformers
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--window-aa", type=int, default=12)
    parser.add_argument("--flank-aa", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument(
        "--include-class",
        action="append",
        default=[],
        help="optional breakpoint_class filter; may be supplied more than once",
    )
    return parser.parse_args(argv)


def load_runtime(model_name: str, requested_device: str):
    try:
        import torch
        from transformers import AutoModelForMaskedLM, AutoTokenizer
    except ImportError as error:  # pragma: no cover - exercised only in optional workflow
        raise RuntimeError("install optional dependencies with: pip install torch transformers") from error

    if requested_device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = requested_device
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    model.eval()
    model.to(device)
    if tokenizer.mask_token_id is None:
        raise RuntimeError(f"tokenizer for {model_name!r} has no mask token")
    return torch, tokenizer, model, device


def sequence_log_probabilities(
    sequence: str,
    positions: Sequence[int],
    torch,
    tokenizer,
    model,
    device: str,
    batch_size: int,
) -> Tuple[float, int]:
    sequence = sequence.replace("*", "X").upper()
    valid_positions = sorted({position for position in positions if 0 <= position < len(sequence)})
    if not sequence or not valid_positions:
        return 0.0, 0

    encoded = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
    input_ids = encoded["input_ids"][0]
    attention = encoded.get("attention_mask")
    attention_row = attention[0] if attention is not None else None

    # ESM tokenizers add a BOS token.  Infer the residue offset rather than assuming it.
    special_before = max(0, input_ids.numel() - len(sequence) - 1)
    residue_offset = 1 if special_before >= 0 else 0
    total = 0.0
    used = 0
    with torch.no_grad():
        for start in range(0, len(valid_positions), batch_size):
            chunk = valid_positions[start : start + batch_size]
            batch_ids = input_ids.repeat(len(chunk), 1)
            batch_attention = (
                attention_row.repeat(len(chunk), 1) if attention_row is not None else None
            )
            token_positions: List[int] = []
            true_tokens: List[int] = []
            for row, residue_position in enumerate(chunk):
                token_position = residue_position + residue_offset
                if token_position >= batch_ids.shape[1]:
                    continue
                token_positions.append(token_position)
                true_tokens.append(int(batch_ids[row, token_position]))
                batch_ids[row, token_position] = tokenizer.mask_token_id
            if not token_positions:
                continue
            batch_ids = batch_ids.to(device)
            kwargs = {"input_ids": batch_ids}
            if batch_attention is not None:
                kwargs["attention_mask"] = batch_attention.to(device)
            logits = model(**kwargs).logits
            log_probs = torch.log_softmax(logits, dim=-1)
            for row, (token_position, true_token) in enumerate(zip(token_positions, true_tokens)):
                total += float(log_probs[row, token_position, true_token].cpu())
                used += 1
    return total, used


def score_candidate(
    peptide: str,
    boundary: int,
    window_aa: int,
    flank_aa: int,
    runtime,
    batch_size: int,
) -> Tuple[float, float, float, int]:
    torch, tokenizer, model, device = runtime
    start = max(0, boundary - flank_aa)
    end = min(len(peptide), boundary + flank_aa)
    local = peptide[start:end]
    local_boundary = boundary - start
    if local_boundary < 2 or len(local) - local_boundary < 2:
        return 0.0, 0.0, 0.0, 0

    joined_positions = list(
        range(max(0, local_boundary - window_aa), min(len(local), local_boundary + window_aa))
    )
    joined_total, joined_count = sequence_log_probabilities(
        local,
        joined_positions,
        torch,
        tokenizer,
        model,
        device,
        batch_size,
    )

    left = local[:local_boundary]
    right = local[local_boundary:]
    left_positions = list(range(max(0, len(left) - window_aa), len(left)))
    right_positions = list(range(0, min(window_aa, len(right))))
    left_total, left_count = sequence_log_probabilities(
        left,
        left_positions,
        torch,
        tokenizer,
        model,
        device,
        batch_size,
    )
    right_total, right_count = sequence_log_probabilities(
        right,
        right_positions,
        torch,
        tokenizer,
        model,
        device,
        batch_size,
    )
    independent_total = left_total + right_total
    independent_count = left_count + right_count
    count = min(joined_count, independent_count)
    if count == 0:
        return 0.0, 0.0, 0.0, 0
    joined_mean = joined_total / max(1, joined_count)
    independent_mean = independent_total / max(1, independent_count)
    delta = joined_mean - independent_mean
    return delta, joined_mean, independent_mean, count


def read_candidates(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source", "target", "junction_peptide", "junction_boundary_aa"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"candidate TSV is missing columns: {sorted(missing)}")
        yield from reader


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        if args.window_aa < 1 or args.flank_aa < args.window_aa:
            raise ValueError("require 1 <= --window-aa <= --flank-aa")
        if args.batch_size < 1:
            raise ValueError("--batch-size must be positive")
        runtime = load_runtime(args.model, args.device)
        included = set(args.include_class)
        written = 0
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "source",
                    "target",
                    "esm_delta",
                    "joined_pll_per_residue",
                    "independent_pll_per_residue",
                    "scored_residues",
                    "model",
                ]
            )
            for row in read_candidates(args.candidates):
                breakpoint_class = row.get("breakpoint_class", "")
                if included and breakpoint_class not in included:
                    continue
                peptide = row.get("junction_peptide", "")
                if not peptide or peptide == ".":
                    continue
                try:
                    boundary = int(row.get("junction_boundary_aa", "-1"))
                except ValueError:
                    continue
                delta, joined, independent, count = score_candidate(
                    peptide=peptide,
                    boundary=boundary,
                    window_aa=args.window_aa,
                    flank_aa=args.flank_aa,
                    runtime=runtime,
                    batch_size=args.batch_size,
                )
                writer.writerow(
                    [
                        row["source"],
                        row["target"],
                        f"{delta:.6f}",
                        f"{joined:.6f}",
                        f"{independent:.6f}",
                        count,
                        args.model,
                    ]
                )
                written += 1
        print(f"ESM breakpoint scoring: rows={written} model={args.model}", file=sys.stderr)
        return 0
    except (OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
