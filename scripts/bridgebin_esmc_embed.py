#!/usr/bin/env python3
"""Extract auditable per-protein ESM-C embeddings for BridgeBin Biological Brain.

The Rust binner never imports ESM-C.  This optional Python adapter accepts a sparse ORF
set (normally produced only for Biological-Brain candidate bins) and emits one row per
protein in the generic format consumed by ``bridgebin_bio_features.py``.

The loader supports the current Biohub Hugging Face API (``EsmcForMaskedLM``) and keeps a
legacy fallback for older ``esm`` releases.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True, help="protein FASTA")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mapping",
        type=Path,
        help="optional TSV with protein_id and contig columns; otherwise FASTA IDs must be contig|protein",
    )
    parser.add_argument("--model", default="biohub/ESMC-300M")
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--min-aa", type=int, default=30)
    parser.add_argument("--max-aa", type=int, default=768)
    parser.add_argument("--truncate", action="store_true")
    parser.add_argument("--batch-size", type=int, default=2)
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
    for separator in ("::", "|"):
        if separator in protein_id:
            return protein_id.split(separator, 1)[0]
    raise ValueError(
        f"cannot infer contig for protein {protein_id!r}; provide --mapping or a contig::protein ID"
    )


def embedding_confidence(sequence: str) -> float:
    informative = sum(residue not in {"X", "B", "Z", "J", "U", "O", "*"} for residue in sequence)
    quality = informative / max(1, len(sequence))
    amount = 1.0 - math.exp(-len(sequence) / 120.0)
    return max(0.0, min(1.0, quality * amount))


def format_vector(values) -> str:
    return ",".join(f"{float(value):.7g}" for value in values)


def resolve_device(torch, requested: str):
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(requested)


def current_runtime(model_name: str, requested_device: str):
    import torch
    from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

    device = resolve_device(torch, requested_device)
    try:
        model = EsmcForMaskedLM.from_pretrained(model_name, device=str(device)).eval()
    except TypeError:
        model = EsmcForMaskedLM.from_pretrained(model_name).to(device).eval()
    tokenizer = EsmcTokenizer()
    return "hf", torch, model, tokenizer, device


def legacy_runtime(model_name: str, requested_device: str):
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    device = resolve_device(torch, requested_device)
    legacy_name = model_name
    aliases = {
        "biohub/ESMC-300M": "esmc_300m",
        "biohub/ESMC-600M": "esmc_600m",
        "biohub/ESMC-6B": "esmc_6b",
    }
    legacy_name = aliases.get(legacy_name, legacy_name)
    model = ESMC.from_pretrained(legacy_name).to(device).eval()
    return "legacy", torch, model, (ESMProtein, LogitsConfig), device


def load_runtime(model_name: str, requested_device: str):
    try:
        return current_runtime(model_name, requested_device)
    except (ImportError, AttributeError):
        try:
            return legacy_runtime(model_name, requested_device)
        except ImportError as error:
            raise RuntimeError("install the optional ESM dependency with: pip install esm") from error


def extract_current_batch(runtime, sequences: List[str]):
    _kind, torch, model, tokenizer, device = runtime
    encoded = tokenizer(
        sequences,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max(len(sequence) for sequence in sequences) + 2,
    )
    encoded = {key: value.to(device) for key, value in encoded.items()}
    with torch.inference_mode():
        output = model(**encoded, output_hidden_states=True)
    hidden_states = getattr(output, "hidden_states", None)
    if not hidden_states:
        hidden = getattr(output, "last_hidden_state", None)
        if hidden is None:
            hidden = getattr(output, "embeddings", None)
    else:
        hidden = hidden_states[-1]
    if hidden is None:
        raise RuntimeError("ESM-C returned no final representations")
    hidden = hidden.float()
    mask = encoded.get("attention_mask")
    if mask is None:
        mask = torch.ones(hidden.shape[:2], dtype=torch.long, device=hidden.device)
    mask_f = mask.to(hidden.dtype).unsqueeze(-1)
    pooled = (hidden * mask_f).sum(dim=1) / mask_f.sum(dim=1).clamp_min(1.0)
    return pooled.cpu().tolist()


def extract_legacy(runtime, sequence: str):
    _kind, torch, model, api, _device = runtime
    ESMProtein, LogitsConfig = api
    with torch.inference_mode():
        encoded = model.encode(ESMProtein(sequence=sequence))
        output = model.logits(encoded, LogitsConfig(sequence=True, return_embeddings=True))
    values = output.embeddings
    if values is None:
        raise RuntimeError("legacy ESM-C returned no embeddings")
    values = values.detach().float().cpu()
    while values.ndim > 2 and values.shape[0] == 1:
        values = values.squeeze(0)
    if values.ndim != 2:
        raise RuntimeError(f"unexpected ESM-C embedding shape {tuple(values.shape)}")
    if values.shape[0] >= len(sequence) + 2:
        values = values[1 : 1 + len(sequence)]
    elif values.shape[0] >= len(sequence):
        values = values[: len(sequence)]
    return values.mean(dim=0).tolist()


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.min_aa < 1 or args.max_aa < args.min_aa or args.batch_size < 1:
        raise SystemExit("require 1 <= --min-aa <= --max-aa and positive --batch-size")
    mapping = read_mapping(args.mapping) if args.mapping else {}
    runtime = load_runtime(args.model, args.device)
    kind = runtime[0]

    records = []
    skipped_short = skipped_long = 0
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
        records.append((contig, protein_id, sequence))
    if not records:
        raise SystemExit("no proteins passed ESM-C filters")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "protein_id", "embedding", "confidence", "aa_length", "model"])
        if kind == "hf":
            for start in range(0, len(records), args.batch_size):
                batch = records[start : start + args.batch_size]
                vectors = extract_current_batch(runtime, [item[2] for item in batch])
                for (contig, protein_id, sequence), vector in zip(batch, vectors):
                    writer.writerow(
                        [contig, protein_id, format_vector(vector), f"{embedding_confidence(sequence):.6f}", len(sequence), args.model]
                    )
                    written += 1
        else:
            for contig, protein_id, sequence in records:
                vector = extract_legacy(runtime, sequence)
                writer.writerow(
                    [contig, protein_id, format_vector(vector), f"{embedding_confidence(sequence):.6f}", len(sequence), args.model]
                )
                written += 1
    print(
        f"bridgebin-esmc: written={written} skipped_short={skipped_short} skipped_long={skipped_long} "
        f"model={args.model} backend={kind} device={runtime[-1]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
