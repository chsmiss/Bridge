#!/usr/bin/env python3
"""Train or run BridgeBin's multimodal same-genome pair head.

The expensive foundation models stay outside Rust.  This script consumes contig-level
embeddings/features and learns the actual decision BridgeBin needs:

    P(contig A and contig B originate from the same genome)

Training deliberately supports asymmetric false-merge cost and reports conservative
positive/negative thresholds on a held-out split.  It has no sklearn dependency so that
the pair head can be audited and reproduced in minimal environments.

Feature TSV columns are header-based. Supported biological fields include:
  dna_embedding / dna_lm_embedding / dnabert_embedding
  gene_profile / gene_embedding
  esm_embedding / esmc_embedding / protein_embedding
  taxonomy / lineage, plus corresponding confidence columns

Pair TSVs require left/right IDs. Training additionally requires label (0/1). Optional
cheap/physical columns are coverage_similarity, composition_similarity, gc_similarity,
and physical_support. Missing modalities are represented explicitly rather than silently
being treated as zero similarity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


BIO_FEATURES = ("dna", "gene", "protein", "taxonomy")
PAIR_FEATURES = ("coverage", "composition", "gc", "physical")
FEATURE_NAMES = tuple(f"{name}_similarity" for name in BIO_FEATURES + PAIR_FEATURES) + tuple(
    f"{name}_present" for name in BIO_FEATURES + PAIR_FEATURES
)


@dataclass
class ContigFeature:
    dna: List[float]
    dna_confidence: float
    gene: List[float]
    gene_confidence: float
    protein: List[float]
    protein_confidence: float
    taxonomy: str
    taxonomy_confidence: float


@dataclass
class Example:
    left: str
    right: str
    features: List[float]
    label: Optional[int]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    train = sub.add_parser("train", help="fit a logistic pair head on labelled pairs")
    train.add_argument("--features", type=Path, required=True)
    train.add_argument("--pairs", type=Path, required=True)
    train.add_argument("--model-out", type=Path, required=True)
    train.add_argument("--validation-fraction", type=float, default=0.20)
    train.add_argument("--epochs", type=int, default=800)
    train.add_argument("--learning-rate", type=float, default=0.08)
    train.add_argument("--l2", type=float, default=1e-3)
    train.add_argument("--negative-weight", type=float, default=3.0)
    train.add_argument("--precision-target", type=float, default=0.995)
    train.add_argument("--seed", type=int, default=43)

    score = sub.add_parser("score", help="score candidate pairs with a trained head")
    score.add_argument("--features", type=Path, required=True)
    score.add_argument("--pairs", type=Path, required=True)
    score.add_argument("--model", type=Path, required=True)
    score.add_argument("--output", type=Path, required=True)
    score.add_argument("--model-name", default="bridgebin-biobrain-pair-head")
    return parser.parse_args(argv)


def rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        yield from reader


def first(row: Dict[str, str], names: Sequence[str], default: str = "") -> str:
    for name in names:
        value = row.get(name)
        if value is not None and value.strip() and value.strip() not in {".", "NA", "na"}:
            return value.strip()
    return default


def parse_vector(raw: str) -> List[float]:
    if not raw:
        return []
    values = [float(piece) for piece in raw.split(",") if piece.strip()]
    if any(not math.isfinite(value) for value in values):
        raise ValueError("embedding contains non-finite value")
    return values


def probability(raw: str, default: float = 0.0) -> float:
    if not raw:
        return default
    value = float(raw)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"probability outside [0,1]: {raw!r}")
    return value


def read_contig_features(path: Path) -> Dict[str, ContigFeature]:
    out: Dict[str, ContigFeature] = {}
    for row in rows(path):
        contig = first(row, ("contig", "contig_id", "sequence"))
        if not contig:
            continue
        dna = parse_vector(first(row, ("dna_embedding", "dna_lm_embedding", "dnabert_embedding")))
        gene = parse_vector(first(row, ("gene_profile", "gene_embedding")))
        protein = parse_vector(
            first(row, ("protein_embedding", "esm_embedding", "esmc_embedding"))
        )
        out[contig] = ContigFeature(
            dna=dna,
            dna_confidence=probability(
                first(row, ("dna_confidence", "dna_lm_confidence")), 1.0 if dna else 0.0
            ),
            gene=gene,
            gene_confidence=probability(
                first(row, ("gene_confidence", "gene_conf")), 1.0 if gene else 0.0
            ),
            protein=protein,
            protein_confidence=probability(
                first(row, ("protein_confidence", "esm_confidence", "protein_conf")),
                1.0 if protein else 0.0,
            ),
            taxonomy=first(row, ("taxonomy", "lineage", "taxon")),
            taxonomy_confidence=probability(
                first(row, ("taxonomy_confidence", "tax_confidence", "tax_conf")),
                0.0,
            ),
        )
    return out


def cosine(left: Sequence[float], right: Sequence[float]) -> Optional[float]:
    if not left or len(left) != len(right):
        return None
    dot = sum(a * b for a, b in zip(left, right))
    lnorm = math.sqrt(sum(value * value for value in left))
    rnorm = math.sqrt(sum(value * value for value in right))
    if lnorm <= 1e-12 or rnorm <= 1e-12:
        return None
    return max(-1.0, min(1.0, dot / (lnorm * rnorm)))


def lineage_parts(value: str) -> List[str]:
    return [
        part.strip().lower()
        for part in value.replace("|", ";").split(";")
        if part.strip() and part.strip() not in {".", "na"}
    ]


def taxonomy_similarity(left: str, right: str) -> Optional[float]:
    a = lineage_parts(left)
    b = lineage_parts(right)
    if not a or not b:
        return None
    shared = 0
    for x, y in zip(a, b):
        if x != y:
            break
        shared += 1
    return shared / max(1, min(len(a), len(b)))


def optional_pair_value(row: Dict[str, str], names: Sequence[str]) -> Optional[float]:
    raw = first(row, names)
    if not raw:
        return None
    value = float(raw)
    if not math.isfinite(value):
        return None
    return max(0.0, min(1.0, value))


def make_features(
    left: ContigFeature,
    right: ContigFeature,
    pair_row: Dict[str, str],
) -> List[float]:
    sims: Dict[str, Optional[float]] = {
        "dna": cosine(left.dna, right.dna),
        "gene": cosine(left.gene, right.gene),
        "protein": cosine(left.protein, right.protein),
        "taxonomy": taxonomy_similarity(left.taxonomy, right.taxonomy),
        "coverage": optional_pair_value(pair_row, ("coverage_similarity", "coverage")),
        "composition": optional_pair_value(
            pair_row, ("composition_similarity", "composition", "tnf_similarity")
        ),
        "gc": optional_pair_value(pair_row, ("gc_similarity", "gc")),
        "physical": optional_pair_value(
            pair_row, ("physical_support", "read_support", "link_support")
        ),
    }
    confidences = {
        "dna": min(left.dna_confidence, right.dna_confidence),
        "gene": min(left.gene_confidence, right.gene_confidence),
        "protein": min(left.protein_confidence, right.protein_confidence),
        "taxonomy": min(left.taxonomy_confidence, right.taxonomy_confidence),
        "coverage": 1.0,
        "composition": 1.0,
        "gc": 1.0,
        "physical": 1.0,
    }
    values = []
    present = []
    for name in BIO_FEATURES + PAIR_FEATURES:
        value = sims[name]
        if value is None:
            values.append(0.0)
            present.append(0.0)
        else:
            # Biological cosine may be negative. Map it to [0,1], then multiply by
            # modality confidence so low-quality annotation cannot masquerade as evidence.
            normalized = (value + 1.0) * 0.5 if name in {"dna", "gene", "protein"} else value
            values.append(max(0.0, min(1.0, normalized)) * confidences[name])
            present.append(confidences[name])
    return values + present


def read_examples(
    feature_path: Path,
    pair_path: Path,
    require_label: bool,
) -> List[Example]:
    features = read_contig_features(feature_path)
    examples = []
    for row in rows(pair_path):
        left = first(row, ("left", "source", "contig_a", "contig1"))
        right = first(row, ("right", "target", "contig_b", "contig2"))
        if not left or not right or left == right or left not in features or right not in features:
            continue
        label_raw = first(row, ("label", "same_genome", "truth"))
        label: Optional[int] = None
        if label_raw:
            parsed = float(label_raw)
            if parsed not in {0.0, 1.0}:
                raise ValueError(f"label must be 0/1, got {label_raw!r}")
            label = int(parsed)
        elif require_label:
            raise ValueError("training pair table requires label/same_genome column")
        examples.append(
            Example(
                left=left,
                right=right,
                features=make_features(features[left], features[right], row),
                label=label,
            )
        )
    if not examples:
        raise ValueError("no usable pair examples")
    return examples


def deterministic_validation(example: Example, fraction: float, seed: int) -> bool:
    key = "\0".join(sorted((example.left, example.right))) + f"\0{seed}"
    digest = hashlib.blake2b(key.encode("utf-8"), digest_size=8).digest()
    unit = int.from_bytes(digest, "little") / float(2**64 - 1)
    return unit < fraction


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def fit_standardizer(examples: Sequence[Example]) -> Tuple[List[float], List[float]]:
    width = len(FEATURE_NAMES)
    means = [0.0] * width
    for example in examples:
        for index, value in enumerate(example.features):
            means[index] += value
    means = [value / len(examples) for value in means]
    variances = [0.0] * width
    for example in examples:
        for index, value in enumerate(example.features):
            variances[index] += (value - means[index]) ** 2
    scales = [math.sqrt(value / len(examples)) for value in variances]
    scales = [value if value > 1e-8 else 1.0 for value in scales]
    return means, scales


def standardized(values: Sequence[float], means: Sequence[float], scales: Sequence[float]) -> List[float]:
    return [(value - mean) / scale for value, mean, scale in zip(values, means, scales)]


def train_model(args: argparse.Namespace) -> int:
    if not 0.0 < args.validation_fraction < 0.5:
        raise ValueError("--validation-fraction must be in (0,0.5)")
    if not 0.5 < args.precision_target < 1.0:
        raise ValueError("--precision-target must be in (0.5,1)")
    examples = read_examples(args.features, args.pairs, require_label=True)
    train = [
        example
        for example in examples
        if not deterministic_validation(example, args.validation_fraction, args.seed)
    ]
    valid = [
        example
        for example in examples
        if deterministic_validation(example, args.validation_fraction, args.seed)
    ]
    if not train or not valid:
        raise ValueError("deterministic split produced empty train or validation set")
    if len({example.label for example in train}) < 2:
        raise ValueError("training split needs both positive and negative labels")

    means, scales = fit_standardizer(train)
    weights = [0.0] * len(FEATURE_NAMES)
    bias = 0.0
    rng = random.Random(args.seed)
    order = list(range(len(train)))
    for epoch in range(args.epochs):
        rng.shuffle(order)
        rate = args.learning_rate / math.sqrt(1.0 + epoch / 50.0)
        grad_w = [0.0] * len(weights)
        grad_b = 0.0
        total_weight = 0.0
        for example_index in order:
            example = train[example_index]
            x = standardized(example.features, means, scales)
            label = int(example.label or 0)
            sample_weight = args.negative_weight if label == 0 else 1.0
            probability_same = sigmoid(bias + sum(w * value for w, value in zip(weights, x)))
            error = (probability_same - label) * sample_weight
            for index, value in enumerate(x):
                grad_w[index] += error * value
            grad_b += error
            total_weight += sample_weight
        inv = 1.0 / max(total_weight, 1.0)
        for index in range(len(weights)):
            grad = grad_w[index] * inv + args.l2 * weights[index]
            weights[index] -= rate * grad
        bias -= rate * grad_b * inv

    validation = []
    for example in valid:
        x = standardized(example.features, means, scales)
        probability_same = sigmoid(bias + sum(w * value for w, value in zip(weights, x)))
        validation.append((probability_same, int(example.label or 0)))

    join_threshold = precision_threshold(validation, args.precision_target, positive=True)
    split_threshold = precision_threshold(validation, args.precision_target, positive=False)
    metrics = classification_metrics(validation, 0.5)
    model = {
        "version": 1,
        "feature_names": list(FEATURE_NAMES),
        "means": means,
        "scales": scales,
        "weights": weights,
        "bias": bias,
        "training": {
            "examples": len(train),
            "validation_examples": len(valid),
            "negative_weight": args.negative_weight,
            "precision_target": args.precision_target,
            "recommended_join_min_same": join_threshold,
            "recommended_split_max_same": split_threshold,
            "validation_at_0_5": metrics,
        },
    }
    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    args.model_out.write_text(json.dumps(model, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "train": len(train),
                "validation": len(valid),
                "join_min_same": join_threshold,
                "split_max_same": split_threshold,
                **metrics,
            },
            sort_keys=True,
        )
    )
    return 0


def precision_threshold(
    predictions: Sequence[Tuple[float, int]], target: float, positive: bool
) -> float:
    if positive:
        candidates = sorted({probability for probability, _ in predictions}, reverse=True)
        best = 1.0
        for threshold in candidates:
            selected = [label for probability, label in predictions if probability >= threshold]
            if selected and sum(selected) / len(selected) >= target:
                best = threshold
        return best
    candidates = sorted({probability for probability, _ in predictions})
    best = 0.0
    for threshold in candidates:
        selected = [label for probability, label in predictions if probability <= threshold]
        if selected and sum(1 - label for label in selected) / len(selected) >= target:
            best = threshold
    return best


def classification_metrics(predictions: Sequence[Tuple[float, int]], threshold: float) -> Dict[str, float]:
    tp = fp = tn = fn = 0
    for probability, label in predictions:
        pred = int(probability >= threshold)
        if pred == 1 and label == 1:
            tp += 1
        elif pred == 1:
            fp += 1
        elif label == 0:
            tn += 1
        else:
            fn += 1
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    specificity = tn / max(1, tn + fp)
    return {"precision": precision, "recall": recall, "specificity": specificity}


def score_model(args: argparse.Namespace) -> int:
    model = json.loads(args.model.read_text(encoding="utf-8"))
    if tuple(model.get("feature_names", ())) != FEATURE_NAMES:
        raise ValueError("pair-head model feature schema does not match this script")
    means = [float(value) for value in model["means"]]
    scales = [float(value) for value in model["scales"]]
    weights = [float(value) for value in model["weights"]]
    bias = float(model["bias"])
    examples = read_examples(args.features, args.pairs, require_label=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["left", "right", "p_same", "confidence", "model"])
        for example in examples:
            x = standardized(example.features, means, scales)
            probability_same = sigmoid(bias + sum(w * value for w, value in zip(weights, x)))
            # Normalized Bernoulli entropy: 0 at p=0/1 and 1 at p=0.5.
            p = min(max(probability_same, 1e-12), 1.0 - 1e-12)
            entropy = -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p)) / math.log(2.0)
            confidence = max(0.0, min(1.0, 1.0 - entropy))
            writer.writerow(
                [
                    example.left,
                    example.right,
                    f"{probability_same:.8f}",
                    f"{confidence:.8f}",
                    args.model_name,
                ]
            )
    print(f"bridgebin-pair-head: scored={len(examples)} output={args.output}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.command == "train":
        return train_model(args)
    if args.command == "score":
        return score_model(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
