#!/usr/bin/env python3
"""Train or run BridgeBin's multimodal same-genome pair head.

The expensive foundation models stay outside Rust. This script consumes contig-level
features and learns the decision BridgeBin actually needs:

    P(contig A and contig B originate from the same genome)

Biological modalities are deliberately separated rather than concatenated blindly:
  dna           species-aware DNA LM / GENERanno-base embedding
  gene          functional gene-family profile when available
  architecture  GENERanno CDS/coding-architecture vector
  protein       pooled ESM-C protein representation
  repertoire    ESM-C protein-prototype TF-IDF profile
  taxonomy      soft lineage agreement

Optional pair-level evidence is coverage, nucleotide composition, GC, and physical
support. Every modality has an explicit presence/confidence feature, so missing model
outputs cannot silently masquerade as zero similarity.

Training supports asymmetric false-merge cost and conservative held-out calibration.
When pair rows contain genome/group columns, validation holds out whole genomes: every
pair touching a held-out genome is excluded from training. This prevents leakage from
putting contigs from the same genome into both train and validation. If no group metadata
is available the script falls back to a deterministic pair split and records that weaker
protocol in the model metadata.

For vector modalities, cosine is clipped to [0,1]. It is intentionally *not* transformed
as (cos+1)/2: orthogonal genome representations must remain zero evidence, not 0.5.
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
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


VECTOR_FEATURES = ("dna", "gene", "architecture", "protein", "repertoire")
BIO_FEATURES = VECTOR_FEATURES + ("taxonomy",)
PAIR_FEATURES = ("coverage", "composition", "gc", "physical")
ALL_MODALITIES = BIO_FEATURES + PAIR_FEATURES
FEATURE_NAMES = tuple(f"{name}_similarity" for name in ALL_MODALITIES) + tuple(
    f"{name}_present" for name in ALL_MODALITIES
)


@dataclass
class ContigFeature:
    dna: List[float]
    dna_confidence: float
    gene: List[float]
    gene_confidence: float
    architecture: List[float]
    architecture_confidence: float
    protein: List[float]
    protein_confidence: float
    repertoire: List[float]
    repertoire_confidence: float
    taxonomy: str
    taxonomy_confidence: float


@dataclass
class Example:
    left: str
    right: str
    features: List[float]
    label: Optional[int]
    left_group: str = ""
    right_group: str = ""


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    train = sub.add_parser("train", help="fit a calibrated logistic same-genome head")
    train.add_argument("--features", type=Path, required=True)
    train.add_argument("--pairs", type=Path, required=True)
    train.add_argument("--model-out", type=Path, required=True)
    train.add_argument("--validation-fraction", type=float, default=0.20)
    train.add_argument("--epochs", type=int, default=800)
    train.add_argument("--learning-rate", type=float, default=0.08)
    train.add_argument("--l2", type=float, default=1e-3)
    train.add_argument(
        "--negative-weight",
        type=float,
        default=3.0,
        help="relative training cost of a false merge / negative pair",
    )
    train.add_argument("--precision-target", type=float, default=0.995)
    train.add_argument("--seed", type=int, default=43)
    train.add_argument(
        "--require-group-holdout",
        action="store_true",
        help="fail instead of pair-splitting if genome/group columns are absent",
    )

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


def feature_confidence(
    row: Dict[str, str], names: Sequence[str], present: bool, default_if_present: float = 1.0
) -> float:
    return probability(first(row, names), default_if_present if present else 0.0)


def read_contig_features(path: Path) -> Dict[str, ContigFeature]:
    out: Dict[str, ContigFeature] = {}
    for row in rows(path):
        contig = first(row, ("contig", "contig_id", "sequence"))
        if not contig:
            continue
        dna = parse_vector(first(row, ("dna_embedding", "dna_lm_embedding", "dnabert_embedding")))
        gene = parse_vector(first(row, ("gene_profile", "gene_embedding", "gene_repertoire")))
        architecture = parse_vector(
            first(row, ("gene_architecture", "architecture_embedding", "coding_architecture"))
        )
        protein = parse_vector(
            first(row, ("protein_embedding", "esm_embedding", "esmc_embedding"))
        )
        repertoire = parse_vector(
            first(row, ("protein_repertoire", "protein_prototypes", "repertoire_embedding"))
        )
        taxonomy = first(row, ("taxonomy", "lineage", "taxon"))
        out[contig] = ContigFeature(
            dna=dna,
            dna_confidence=feature_confidence(
                row, ("dna_confidence", "dna_lm_confidence"), bool(dna)
            ),
            gene=gene,
            gene_confidence=feature_confidence(
                row, ("gene_confidence", "gene_conf", "gene_repertoire_confidence"), bool(gene)
            ),
            architecture=architecture,
            architecture_confidence=feature_confidence(
                row,
                ("architecture_confidence", "gene_architecture_confidence"),
                bool(architecture),
            ),
            protein=protein,
            protein_confidence=feature_confidence(
                row,
                ("protein_confidence", "esm_confidence", "protein_conf"),
                bool(protein),
            ),
            repertoire=repertoire,
            repertoire_confidence=feature_confidence(
                row,
                ("repertoire_confidence", "protein_repertoire_confidence"),
                bool(repertoire),
            ),
            taxonomy=taxonomy,
            taxonomy_confidence=feature_confidence(
                row,
                ("taxonomy_confidence", "tax_confidence", "tax_conf"),
                bool(taxonomy),
                default_if_present=0.5,
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


def vector_similarity(left: Sequence[float], right: Sequence[float]) -> Optional[float]:
    value = cosine(left, right)
    if value is None:
        return None
    # Negative/orthogonal representations are not positive same-genome evidence.
    return max(0.0, value)


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
    left: ContigFeature, right: ContigFeature, pair_row: Dict[str, str]
) -> List[float]:
    sims: Dict[str, Optional[float]] = {
        "dna": vector_similarity(left.dna, right.dna),
        "gene": vector_similarity(left.gene, right.gene),
        "architecture": vector_similarity(left.architecture, right.architecture),
        "protein": vector_similarity(left.protein, right.protein),
        "repertoire": vector_similarity(left.repertoire, right.repertoire),
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
        "architecture": min(left.architecture_confidence, right.architecture_confidence),
        "protein": min(left.protein_confidence, right.protein_confidence),
        "repertoire": min(left.repertoire_confidence, right.repertoire_confidence),
        "taxonomy": min(left.taxonomy_confidence, right.taxonomy_confidence),
        "coverage": 1.0,
        "composition": 1.0,
        "gc": 1.0,
        "physical": 1.0,
    }
    values: List[float] = []
    present: List[float] = []
    for name in ALL_MODALITIES:
        value = sims[name]
        confidence = confidences[name]
        if value is None or confidence <= 0.0:
            values.append(0.0)
            present.append(0.0)
        else:
            values.append(max(0.0, min(1.0, value)) * confidence)
            present.append(confidence)
    return values + present


def read_examples(feature_path: Path, pair_path: Path, require_label: bool) -> List[Example]:
    features = read_contig_features(feature_path)
    examples: List[Example] = []
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
                left_group=first(
                    row, ("left_group", "left_genome", "genome_left", "source_genome")
                ),
                right_group=first(
                    row, ("right_group", "right_genome", "genome_right", "target_genome")
                ),
            )
        )
    if not examples:
        raise ValueError("no usable pair examples")
    return examples


def hash_unit(value: str, seed: int) -> float:
    digest = hashlib.blake2b(f"{value}\0{seed}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") / float(2**64 - 1)


def pair_validation(example: Example, fraction: float, seed: int) -> bool:
    return hash_unit("\0".join(sorted((example.left, example.right))), seed) < fraction


def split_examples(
    examples: Sequence[Example], fraction: float, seed: int, require_group_holdout: bool
) -> Tuple[List[Example], List[Example], Dict[str, object]]:
    grouped = [example for example in examples if example.left_group and example.right_group]
    groups: Set[str] = {
        group
        for example in grouped
        for group in (example.left_group, example.right_group)
    }
    if len(grouped) == len(examples) and len(groups) >= 3:
        ordered_groups = sorted(groups, key=lambda group: (hash_unit(group, seed), group))
        holdout_count = max(
            1, min(len(ordered_groups) - 2, round(len(ordered_groups) * fraction))
        )
        validation_groups = set(ordered_groups[:holdout_count])
        train = [
            example
            for example in examples
            if example.left_group not in validation_groups
            and example.right_group not in validation_groups
        ]
        valid = [
            example
            for example in examples
            if example.left_group in validation_groups
            or example.right_group in validation_groups
        ]
        protocol = {
            "split_protocol": "whole_genome_holdout",
            "validation_groups": sorted(validation_groups),
            "all_groups": len(groups),
        }
    else:
        if require_group_holdout:
            raise ValueError(
                "--require-group-holdout needs group metadata on every row and at least 3 groups"
            )
        train = [example for example in examples if not pair_validation(example, fraction, seed)]
        valid = [example for example in examples if pair_validation(example, fraction, seed)]
        protocol = {
            "split_protocol": "pair_hash_fallback",
            "validation_groups": [],
            "all_groups": len(groups),
        }
    if not train or not valid:
        raise ValueError("validation split produced empty train or validation set")
    return train, valid, protocol


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
    return means, [value if value > 1e-8 else 1.0 for value in scales]


def standardized(
    values: Sequence[float], means: Sequence[float], scales: Sequence[float]
) -> List[float]:
    return [(value - mean) / scale for value, mean, scale in zip(values, means, scales)]


def train_model(args: argparse.Namespace) -> int:
    if not 0.0 < args.validation_fraction < 0.5:
        raise ValueError("--validation-fraction must be in (0,0.5)")
    if not 0.5 < args.precision_target < 1.0:
        raise ValueError("--precision-target must be in (0.5,1)")
    if args.negative_weight <= 0.0:
        raise ValueError("--negative-weight must be positive")
    examples = read_examples(args.features, args.pairs, require_label=True)
    train, valid, protocol = split_examples(
        examples, args.validation_fraction, args.seed, args.require_group_holdout
    )
    if len({example.label for example in train}) < 2:
        raise ValueError("training split needs both positive and negative labels")
    if len({example.label for example in valid}) < 2:
        raise ValueError("validation split needs both positive and negative labels")

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
            weights[index] -= rate * (grad_w[index] * inv + args.l2 * weights[index])
        bias -= rate * grad_b * inv

    validation: List[Tuple[float, int]] = []
    for example in valid:
        x = standardized(example.features, means, scales)
        probability_same = sigmoid(bias + sum(w * value for w, value in zip(weights, x)))
        validation.append((probability_same, int(example.label or 0)))

    join_threshold = precision_threshold(validation, args.precision_target, positive=True)
    split_threshold = precision_threshold(validation, args.precision_target, positive=False)
    metrics = classification_metrics(validation, 0.5)
    modality_weights = {
        name: weights[index]
        for index, name in enumerate(FEATURE_NAMES)
        if name.endswith("_similarity")
    }
    model = {
        "version": 3,
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
            "learned_similarity_weights": modality_weights,
            **protocol,
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
                **protocol,
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


def classification_metrics(
    predictions: Sequence[Tuple[float, int]], threshold: float
) -> Dict[str, float]:
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
    false_merge_rate = fp / max(1, fp + tn)
    return {
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "false_merge_rate": false_merge_rate,
    }


def score_model(args: argparse.Namespace) -> int:
    model = json.loads(args.model.read_text(encoding="utf-8"))
    if tuple(model.get("feature_names", ())) != FEATURE_NAMES:
        raise ValueError(
            "pair-head model feature schema does not match this script; retrain after modality changes"
        )
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
