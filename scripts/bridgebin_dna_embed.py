#!/usr/bin/env python3
"""Extract orientation-robust DNA foundation-model embeddings for BridgeBin.

The default model is a Hugging Face DNABERT-S implementation, but any AutoModel that
returns ``last_hidden_state`` can be supplied. Long contigs are represented by overlapping
windows; each window is encoded in both forward and reverse-complement orientation and
the two pooled embeddings are averaged. The resulting per-window TSV is consumed by
``bridgebin_bio_features.py --dna-embeddings``.

Optional dependencies: torch, transformers.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", default="multimolecule/dnaberts")
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--window-bp", type=int, default=2048)
    parser.add_argument("--stride-bp", type=int, default=1024)
    parser.add_argument("--min-window-bp", type=int, default=512)
    parser.add_argument("--max-tokens", type=int, default=512)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--trust-remote-code", action="store_true")
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


def windows(sequence: str, window_bp: int, stride_bp: int, min_window_bp: int):
    if len(sequence) <= window_bp:
        if len(sequence) >= min_window_bp:
            yield 0, sequence
        return
    start = 0
    serial = 0
    while start < len(sequence):
        piece = sequence[start : start + window_bp]
        if len(piece) < min_window_bp:
            break
        yield serial, piece
        serial += 1
        if start + window_bp >= len(sequence):
            break
        start += stride_bp


def load_runtime(model_name: str, requested_device: str, trust_remote_code: bool):
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as error:  # pragma: no cover - optional dependency
        raise RuntimeError("install optional dependencies: pip install torch transformers") from error

    device = (
        "cuda"
        if requested_device == "auto" and torch.cuda.is_available()
        else "cpu"
        if requested_device == "auto"
        else requested_device
    )
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=trust_remote_code)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=trust_remote_code)
    model.to(device)
    model.eval()
    return torch, tokenizer, model, device


def mean_pool(hidden, attention_mask, torch):
    mask = attention_mask.clone()
    if mask.shape[1] > 0:
        mask[:, 0] = 0  # CLS/BOS
    lengths = attention_mask.sum(dim=1)
    for row, length in enumerate(lengths.tolist()):
        sep = int(length) - 1
        if 0 <= sep < mask.shape[1]:
            mask[row, sep] = 0
    weights = mask.unsqueeze(-1).to(hidden.dtype)
    summed = (hidden * weights).sum(dim=1)
    denom = weights.sum(dim=1).clamp_min(1.0)
    return summed / denom


def encode_sequences(sequences: Sequence[str], runtime, max_tokens: int):
    torch, tokenizer, model, device = runtime
    encoded = tokenizer(
        list(sequences),
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_tokens,
    )
    encoded = {key: value.to(device) for key, value in encoded.items()}
    with torch.no_grad():
        output = model(**encoded)
    hidden = output.last_hidden_state
    pooled = mean_pool(hidden, encoded["attention_mask"], torch)
    return pooled.detach().float().cpu()


def confidence(sequence: str) -> float:
    valid = sum(base in "ACGT" for base in sequence)
    quality = valid / max(1, len(sequence))
    amount = 1.0 - math.exp(-len(sequence) / 700.0)
    return max(0.0, min(1.0, quality * amount))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not (1 <= args.min_window_bp <= args.window_bp):
        raise SystemExit("require 1 <= --min-window-bp <= --window-bp")
    if not (1 <= args.stride_bp <= args.window_bp):
        raise SystemExit("require 1 <= --stride-bp <= --window-bp")
    if args.max_tokens < 4 or args.batch_size < 1:
        raise SystemExit("--max-tokens and --batch-size are too small")

    runtime = load_runtime(args.model, args.device, args.trust_remote_code)
    records = []
    for contig, sequence in read_fasta(args.contigs):
        for serial, piece in windows(sequence, args.window_bp, args.stride_bp, args.min_window_bp):
            records.append((contig, serial, piece))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "window_id", "embedding", "confidence", "model"])
        for start in range(0, len(records), args.batch_size):
            batch = records[start : start + args.batch_size]
            forward = [piece for _, _, piece in batch]
            reverse = [reverse_complement(piece) for piece in forward]
            forward_embeddings = encode_sequences(forward, runtime, args.max_tokens)
            reverse_embeddings = encode_sequences(reverse, runtime, args.max_tokens)
            averaged = (forward_embeddings + reverse_embeddings) / 2.0
            for (contig, serial, piece), embedding in zip(batch, averaged):
                writer.writerow(
                    [
                        contig,
                        f"w{serial}",
                        ",".join(f"{float(value):.7g}" for value in embedding.tolist()),
                        f"{confidence(piece):.6f}",
                        args.model,
                    ]
                )
                written += 1
    print(f"bridgebin-dna: windows={written} contigs={len({c for c, _, _ in records})} model={args.model}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
