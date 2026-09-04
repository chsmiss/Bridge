#!/usr/bin/env python3
"""Extract orientation-robust DNA foundation-model embeddings for BridgeBin.

The default is the DNABERT-S authors' public checkpoint ``zhihan1996/DNABERT-S``, loaded
with Transformers custom code as documented by the model repository. MultiMolecule models
remain supported explicitly when a ``multimolecule/*`` name is requested.

Long contigs are represented by overlapping windows; each selected window is encoded in
forward and reverse-complement orientation and the pooled representations are averaged.
For expensive foundation models ``--max-windows-per-contig`` can evenly subsample the
full contig rather than spending inference on every overlapping window.  A value of zero
keeps every window.  The resulting per-window TSV is consumed by
``bridgebin_bio_features.py --dna-embeddings``.

Optional dependencies:
  default DNABERT-S:  pip install torch transformers einops
  MultiMolecule:      pip install torch multimolecule
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
    parser.add_argument(
        "--pairs",
        type=Path,
        help=(
            "optional candidate-pair TSV; only contigs appearing in left/right (or source/target) "
            "columns are embedded"
        ),
    )
    parser.add_argument("--model", default="zhihan1996/DNABERT-S")
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--window-bp", type=int, default=2048)
    parser.add_argument("--stride-bp", type=int, default=1024)
    parser.add_argument("--min-window-bp", type=int, default=512)
    parser.add_argument(
        "--max-windows-per-contig",
        type=int,
        default=0,
        help="evenly sample at most N windows per contig; 0 keeps all windows",
    )
    parser.add_argument("--max-tokens", type=int, default=512)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument(
        "--trust-remote-code",
        action="store_true",
        help="allow custom Transformers code; automatically enabled for zhihan1996/DNABERT-S",
    )
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


def read_pair_endpoints(path: Optional[Path]) -> Optional[set[str]]:
    if path is None:
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        fields = {name.lower(): name for name in reader.fieldnames}
        left_name = next((fields[name] for name in ("left", "source", "contig_a", "contig1") if name in fields), None)
        right_name = next((fields[name] for name in ("right", "target", "contig_b", "contig2") if name in fields), None)
        if left_name is None or right_name is None:
            raise ValueError(f"{path}: pair table needs left/right or source/target columns")
        endpoints: set[str] = set()
        for row in reader:
            left = (row.get(left_name) or "").strip()
            right = (row.get(right_name) or "").strip()
            if left:
                endpoints.add(left)
            if right:
                endpoints.add(right)
    return endpoints


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


def evenly_sample_windows(
    candidates: Sequence[Tuple[int, str]], limit: int
) -> List[Tuple[int, str]]:
    if limit <= 0 or len(candidates) <= limit:
        return list(candidates)
    if limit == 1:
        return [candidates[len(candidates) // 2]]
    # Include both ends and spread the remaining samples across the full contig.  The
    # integer formulation is deterministic and cannot select the same index twice when
    # limit <= len(candidates).
    indices = [round(i * (len(candidates) - 1) / (limit - 1)) for i in range(limit)]
    return [candidates[index] for index in indices]


def load_runtime(model_name: str, requested_device: str, trust_remote_code: bool):
    try:
        import torch
    except ImportError as error:  # pragma: no cover - optional dependency
        raise RuntimeError("install the optional runtime with: pip install torch") from error

    device = (
        "cuda"
        if requested_device == "auto" and torch.cuda.is_available()
        else "cpu"
        if requested_device == "auto"
        else requested_device
    )

    if model_name.startswith("multimolecule/"):
        try:
            from multimolecule import AutoModel, AutoTokenizer
        except ImportError as error:  # pragma: no cover - optional dependency
            raise RuntimeError(
                "MultiMolecule DNA models require: pip install torch multimolecule"
            ) from error
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name)
        backend = "multimolecule"
    else:
        try:
            from transformers import AutoModel, AutoTokenizer
        except ImportError as error:  # pragma: no cover - optional dependency
            raise RuntimeError(
                "generic DNA models require: pip install torch transformers"
            ) from error
        # The authors' DNABERT-S checkpoint ships its MosaicBERT/ALiBi implementation
        # as repository custom code. Make that one well-known model work out of the box;
        # other custom repositories still require an explicit opt-in.
        allow_remote = trust_remote_code or model_name.lower() == "zhihan1996/dnabert-s"
        tokenizer = AutoTokenizer.from_pretrained(
            model_name, trust_remote_code=allow_remote
        )
        model = AutoModel.from_pretrained(
            model_name, trust_remote_code=allow_remote
        )
        backend = "transformers"

    model.to(device)
    model.eval()
    return torch, tokenizer, model, device, backend


def mean_pool(hidden, attention_mask, torch):
    mask = attention_mask.clone()
    if mask.shape[1] > 0:
        mask[:, 0] = 0  # exclude CLS/BOS when present
    lengths = attention_mask.sum(dim=1)
    for row, length in enumerate(lengths.tolist()):
        terminal = int(length) - 1
        if 0 <= terminal < mask.shape[1]:
            mask[row, terminal] = 0  # exclude terminal SEP/EOS when present
    weights = mask.unsqueeze(-1).to(hidden.dtype)
    summed = (hidden * weights).sum(dim=1)
    denom = weights.sum(dim=1).clamp_min(1.0)
    return summed / denom


def encode_sequences(sequences: Sequence[str], runtime, max_tokens: int):
    torch, tokenizer, model, device, _backend = runtime
    encoded = tokenizer(
        list(sequences),
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_tokens,
    )
    if "attention_mask" not in encoded:
        encoded["attention_mask"] = torch.ones_like(encoded["input_ids"])
    encoded = {key: value.to(device) for key, value in encoded.items()}
    with torch.no_grad():
        output = model(**encoded)
    hidden = output.last_hidden_state if hasattr(output, "last_hidden_state") else output[0]
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
    if args.max_windows_per_contig < 0:
        raise SystemExit("--max-windows-per-contig must be >= 0")
    if args.max_tokens < 4 or args.batch_size < 1:
        raise SystemExit("--max-tokens and --batch-size are too small")

    endpoints = read_pair_endpoints(args.pairs)
    runtime = load_runtime(args.model, args.device, args.trust_remote_code)
    records = []
    full_window_count = 0
    selected_contigs = 0
    for contig, sequence in read_fasta(args.contigs):
        if endpoints is not None and contig not in endpoints:
            continue
        selected_contigs += 1
        candidates = list(
            windows(sequence, args.window_bp, args.stride_bp, args.min_window_bp)
        )
        full_window_count += len(candidates)
        for serial, piece in evenly_sample_windows(
            candidates, args.max_windows_per_contig
        ):
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
    backend = runtime[-1]
    embedded_contigs = len({c for c, _, _ in records})
    missing_endpoints = 0 if endpoints is None else max(0, len(endpoints) - selected_contigs)
    print(
        f"bridgebin-dna: windows={written}/{full_window_count} "
        f"contigs={embedded_contigs} selected={selected_contigs} "
        f"missing_pair_endpoints={missing_endpoints} model={args.model} backend={backend}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())