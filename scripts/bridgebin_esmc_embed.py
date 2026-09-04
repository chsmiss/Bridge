#!/usr/bin/env python3
"""Extract auditable per-protein ESM-C embeddings for BridgeBin Biological Brain.

This is an optional GPU/Python adapter. The Rust binner never imports ESM-C directly.
Output is compatible with ``scripts/bridgebin_bio_features.py --protein-embeddings``.

Current EvolutionaryScale/Biohub local ESM-C API is used when available::

    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    model = ESMC.from_pretrained("esmc_600m").to("cuda")

Install the optional dependency with ``pip install esm`` and ensure model weights are
available. For reproducibility, the exact model string is written into the output.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True, help="protein FASTA")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mapping",
        type=Path,
        help="optional TSV with protein_id and contig columns; otherwise FASTA IDs must be contig|protein",
    )
    parser.add_argument("--model", default="esmc_600m")
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--min-aa", type=int, default=30)
    parser.add_argument("--max-aa", type=int, default=4096)
    parser.add_argument("--truncate", action="store_true")
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


def read_mapping(path: Path) -> Dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        fields = {name.lower(): name for name in reader.fieldnames}
        protein_col = next(
            (fields[name] for name in ("protein_id", "protein", "orf", "gene_id") if name in fields),
            None,
        )
        contig_col = next(
            (fields[name] for name in ("contig", "contig_id", "sequence") if name in fields),
            None,
        )
        if protein_col is None or contig_col is None:
            raise ValueError("mapping TSV needs protein_id and contig columns")
        return {
            row[protein_col].strip(): row[contig_col].strip()
            for row in reader
            if row.get(protein_col, "").strip() and row.get(contig_col, "").strip()
        }


def infer_contig(protein_id: str) -> str:
    if "|" not in protein_id:
        raise ValueError(
            f"cannot infer contig for protein {protein_id!r}; provide --mapping or use contig|protein FASTA IDs"
        )
    return protein_id.split("|", 1)[0]


def load_model(model_name: str, requested_device: str):
    try:
        import torch
        from esm.models.esmc import ESMC
        from esm.sdk.api import ESMProtein, LogitsConfig
    except ImportError as error:  # pragma: no cover - optional GPU dependency
        raise RuntimeError("install the optional ESM package with: pip install esm") from error

    if requested_device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = requested_device
    model = ESMC.from_pretrained(model_name).to(device)
    model.eval()
    return torch, ESMProtein, LogitsConfig, model, device


def mean_embedding(sequence: str, runtime) -> List[float]:
    torch, ESMProtein, LogitsConfig, model, _device = runtime
    protein = ESMProtein(sequence=sequence)
    with torch.no_grad():
        encoded = model.encode(protein)
        output = model.logits(encoded, LogitsConfig(sequence=True, return_embeddings=True))
    embedding = output.embeddings
    if embedding is None:
        raise RuntimeError("ESM-C returned no embeddings")
    values = embedding.detach().float().cpu()
    while values.ndim > 2 and values.shape[0] == 1:
        values = values.squeeze(0)
    if values.ndim != 2:
        raise RuntimeError(f"unexpected ESM-C embedding shape {tuple(values.shape)}")

    # ESM-C representations normally include sequence special tokens. Infer residue rows
    # from the supplied protein length rather than assuming one exact SDK version.
    if values.shape[0] >= len(sequence) + 2:
        values = values[1 : 1 + len(sequence)]
    elif values.shape[0] >= len(sequence):
        values = values[: len(sequence)]
    if values.shape[0] == 0:
        raise RuntimeError("ESM-C produced an empty residue representation")
    pooled = values.mean(dim=0)
    return [float(value) for value in pooled.tolist()]


def embedding_confidence(sequence: str) -> float:
    informative = sum(residue not in {"X", "B", "Z", "J", "U", "O", "*"} for residue in sequence)
    quality = informative / max(1, len(sequence))
    amount = 1.0 - math.exp(-len(sequence) / 120.0)
    return max(0.0, min(1.0, quality * amount))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.min_aa < 1 or args.max_aa < args.min_aa:
        raise SystemExit("require 1 <= --min-aa <= --max-aa")
    mapping = read_mapping(args.mapping) if args.mapping else {}
    runtime = load_model(args.model, args.device)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    written = skipped_short = skipped_long = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "protein_id", "embedding", "confidence", "model"])
        for protein_id, original_sequence in read_fasta(args.proteins):
            sequence = original_sequence.rstrip("*")
            if len(sequence) < args.min_aa:
                skipped_short += 1
                continue
            if len(sequence) > args.max_aa:
                if not args.truncate:
                    skipped_long += 1
                    continue
                sequence = sequence[: args.max_aa]
            contig = mapping.get(protein_id) or infer_contig(protein_id)
            embedding = mean_embedding(sequence, runtime)
            writer.writerow(
                [
                    contig,
                    protein_id,
                    ",".join(f"{value:.7g}" for value in embedding),
                    f"{embedding_confidence(sequence):.6f}",
                    args.model,
                ]
            )
            written += 1
    print(
        f"bridgebin-esmc: written={written} skipped_short={skipped_short} "
        f"skipped_long={skipped_long} model={args.model}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
