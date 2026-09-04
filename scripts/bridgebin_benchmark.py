#!/usr/bin/env python3
import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

BASES = "ACGT"

MODELS = {
    "g0": {
        "A": [8, 1, 1, 6],
        "C": [5, 2, 1, 5],
        "G": [5, 1, 2, 5],
        "T": [6, 1, 1, 8],
    },
    "g1": {
        "A": [8, 1, 1, 6],
        "C": [5, 2, 1, 5],
        "G": [5, 1, 2, 5],
        "T": [6, 1, 1, 8],
    },
    "g2": {
        "A": [1, 4, 4, 1],
        "C": [4, 1, 1, 4],
        "G": [4, 1, 1, 4],
        "T": [1, 4, 4, 1],
    },
    "g3": {
        "A": [6, 1, 1, 2],
        "C": [1, 6, 5, 1],
        "G": [1, 5, 6, 1],
        "T": [2, 1, 1, 6],
    },
    "g4": {
        "A": [1, 4, 3, 1],
        "C": [1, 8, 6, 1],
        "G": [1, 6, 8, 1],
        "T": [1, 3, 4, 1],
    },
    "g5": {
        "A": [1, 3, 5, 1],
        "C": [2, 5, 8, 1],
        "G": [1, 8, 5, 2],
        "T": [1, 5, 3, 1],
    },
}

PROFILES = {
    "g0": [35.0, 5.0, 18.0, 42.0, 3.0, 20.0],
    "g1": [6.0, 32.0, 4.0, 10.0, 28.0, 7.0],
    "g2": [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
    "g3": [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
    "g4": [4.0, 10.0, 28.0, 6.0, 20.0, 32.0],
    "g5": [5.0, 11.0, 26.0, 7.0, 19.0, 30.0],
}

LENGTHS = [
    7200,
    6400,
    5600,
    4800,
    4200,
    3600,
    3200,
    2800,
    2500,
    2300,
    2100,
    1900,
    1700,
    1400,
    1200,
]
SAMPLES = [f"sample{i}" for i in range(1, 7)]


def weighted_choice(rng, weights):
    total = float(sum(weights))
    x = rng.random() * total
    acc = 0.0
    for base, weight in zip(BASES, weights):
        acc += weight
        if x <= acc:
            return base
    return BASES[-1]


def generate_sequence(rng, model, length):
    prev = rng.choice(BASES)
    seq = [prev]
    for _ in range(length - 1):
        prev = weighted_choice(rng, model[prev])
        seq.append(prev)
    return "".join(seq)


def cmd_generate(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    truth_rows = []
    coverage_rows = []

    with (out / "contigs.fa").open("w") as fasta:
        for genome in sorted(MODELS):
            for idx, length in enumerate(LENGTHS, start=1):
                contig = f"{genome}_c{idx:02d}"
                local_rng = random.Random(rng.getrandbits(64))
                seq = generate_sequence(local_rng, MODELS[genome], length)
                fasta.write(f">{contig}\n")
                for start in range(0, len(seq), 80):
                    fasta.write(seq[start : start + 80] + "\n")
                truth_rows.append((contig, genome, length, int(length >= args.min_contig)))

                contig_scale = 1.0 + rng.uniform(-0.08, 0.08)
                depths = []
                for base_depth in PROFILES[genome]:
                    sample_noise = 1.0 + rng.uniform(-0.02, 0.02)
                    depths.append(base_depth * contig_scale * sample_noise)
                coverage_rows.append((contig, depths))

    with (out / "coverage.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", *SAMPLES])
        for contig, depths in coverage_rows:
            writer.writerow([contig, *[f"{x:.6f}" for x in depths]])

    with (out / "truth.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", "genome", "length", "eligible"])
        writer.writerows(truth_rows)

    with (out / "truth_abundance.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome", "sample", "depth", "relative_abundance"])
        for sample_idx, sample in enumerate(SAMPLES):
            total = sum(PROFILES[genome][sample_idx] for genome in PROFILES)
            for genome in sorted(PROFILES):
                depth = PROFILES[genome][sample_idx]
                writer.writerow(
                    [genome, sample, f"{depth:.6f}", f"{depth / total:.8f}"]
                )

    print(
        f"generated {len(truth_rows)} contigs across {len(MODELS)} genomes in {out}"
    )


def comb2(n):
    return n * (n - 1) / 2.0


def adjusted_rand_index(true_labels, pred_labels):
    n = len(true_labels)
    if n < 2:
        return 1.0
    table = defaultdict(Counter)
    true_counts = Counter(true_labels)
    pred_counts = Counter(pred_labels)
    for truth, pred in zip(true_labels, pred_labels):
        table[truth][pred] += 1
    sum_nij = sum(comb2(v) for row in table.values() for v in row.values())
    sum_ai = sum(comb2(v) for v in true_counts.values())
    sum_bj = sum(comb2(v) for v in pred_counts.values())
    total = comb2(n)
    expected = sum_ai * sum_bj / total if total else 0.0
    maximum = 0.5 * (sum_ai + sum_bj)
    denom = maximum - expected
    if abs(denom) < 1e-15:
        return 1.0
    return (sum_nij - expected) / denom


def pearson(xs, ys):
    if len(xs) != len(ys) or not xs:
        return 0.0
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denom_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    denom_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if denom_x == 0.0 or denom_y == 0.0:
        return 0.0
    return numerator / (denom_x * denom_y)


def load_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cmd_evaluate(args):
    truth_rows = load_tsv(args.truth)
    assignment_rows = load_tsv(args.assignments)
    abundance_rows = load_tsv(args.abundance)
    truth_abundance_rows = load_tsv(args.truth_abundance)

    truth = {
        row["contig"]: {
            "genome": row["genome"],
            "length": int(row["length"]),
            "eligible": row.get("eligible", "1") == "1",
        }
        for row in truth_rows
    }
    assignments = {row["contig"]: row["bin"] for row in assignment_rows}

    eligible = [contig for contig, info in truth.items() if info["eligible"]]
    eligible_bp = sum(truth[contig]["length"] for contig in eligible)
    bin_genome_bp = defaultdict(Counter)
    genome_bin_bp = defaultdict(Counter)
    binned_bp = 0
    true_labels = []
    pred_labels = []

    for contig in eligible:
        info = truth[contig]
        pred = assignments.get(contig, "unbinned")
        if pred != "unbinned":
            bp = info["length"]
            binned_bp += bp
            bin_genome_bp[pred][info["genome"]] += bp
            genome_bin_bp[info["genome"]][pred] += bp
        true_labels.append(info["genome"])
        pred_labels.append(pred if pred != "unbinned" else f"unbinned::{contig}")

    genome_total_bp = Counter()
    for contig in eligible:
        genome_total_bp[truth[contig]["genome"]] += truth[contig]["length"]

    majority_bp = sum(max(counter.values()) for counter in bin_genome_bp.values() if counter)
    weighted_purity = majority_bp / binned_bp if binned_bp else 0.0

    best_bin = {}
    completeness_values = []
    captured_bp = 0
    hq = 0
    mq = 0
    for genome in sorted(genome_total_bp):
        counter = genome_bin_bp[genome]
        if counter:
            bin_name, overlap = max(counter.items(), key=lambda item: (item[1], item[0]))
            best_bin[genome] = bin_name
        else:
            bin_name, overlap = None, 0
            best_bin[genome] = None
        completeness = overlap / genome_total_bp[genome]
        completeness_values.append(completeness)
        captured_bp += overlap
        bin_purity = 0.0
        if bin_name is not None:
            total_bin_bp = sum(bin_genome_bp[bin_name].values())
            bin_purity = overlap / total_bin_bp if total_bin_bp else 0.0
        if completeness >= 0.90 and bin_purity >= 0.95:
            hq += 1
        if completeness >= 0.50 and bin_purity >= 0.90:
            mq += 1

    bp_recall = captured_bp / eligible_bp if eligible_bp else 0.0
    if weighted_purity + bp_recall > 0:
        f1 = 2.0 * weighted_purity * bp_recall / (weighted_purity + bp_recall)
    else:
        f1 = 0.0

    abundance = {
        (row["bin"], row["sample"]): {
            "robust_depth": float(row["robust_depth"]),
            "relative_abundance": float(row["relative_abundance"]),
        }
        for row in abundance_rows
    }
    truth_abundance = {
        (row["genome"], row["sample"]): {
            "depth": float(row["depth"]),
            "relative_abundance": float(row["relative_abundance"]),
        }
        for row in truth_abundance_rows
    }

    depth_abs_pct = []
    rel_abs_err = []
    depth_truth_values = []
    depth_est_values = []
    for genome, bin_name in best_bin.items():
        if bin_name is None:
            continue
        for sample in SAMPLES:
            truth_value = truth_abundance[(genome, sample)]
            estimate = abundance.get((bin_name, sample))
            if estimate is None:
                continue
            if truth_value["depth"] > 0:
                depth_abs_pct.append(
                    abs(estimate["robust_depth"] - truth_value["depth"])
                    / truth_value["depth"]
                )
            rel_abs_err.append(
                abs(
                    estimate["relative_abundance"]
                    - truth_value["relative_abundance"]
                )
            )
            depth_truth_values.append(truth_value["depth"])
            depth_est_values.append(estimate["robust_depth"])

    metrics = {
        "eligible_contigs": len(eligible),
        "eligible_genomes": len(genome_total_bp),
        "predicted_bins": len(bin_genome_bp),
        "bp_binned_fraction": binned_bp / eligible_bp if eligible_bp else 0.0,
        "weighted_purity": weighted_purity,
        "bp_recall": bp_recall,
        "f1": f1,
        "ari": adjusted_rand_index(true_labels, pred_labels),
        "mean_genome_completeness": sum(completeness_values)
        / len(completeness_values),
        "hq_like_genomes_90_95": hq,
        "mq_like_genomes_50_90": mq,
        "depth_mape": sum(depth_abs_pct) / len(depth_abs_pct)
        if depth_abs_pct
        else 1.0,
        "relative_abundance_mae": sum(rel_abs_err) / len(rel_abs_err)
        if rel_abs_err
        else 1.0,
        "depth_pearson": pearson(depth_truth_values, depth_est_values),
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    tsv = output.with_suffix(".tsv")
    with tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            if isinstance(value, float):
                writer.writerow([key, f"{value:.8f}"])
            else:
                writer.writerow([key, value])

    for key, value in metrics.items():
        if isinstance(value, float):
            print(f"{key}\t{value:.6f}")
        else:
            print(f"{key}\t{value}")


def main():
    parser = argparse.ArgumentParser(description="Deterministic BridgeBin v0 benchmark")
    sub = parser.add_subparsers(dest="command", required=True)

    generate = sub.add_parser("generate")
    generate.add_argument("--out", required=True)
    generate.add_argument("--seed", type=int, default=43)
    generate.add_argument("--min-contig", type=int, default=1500)
    generate.set_defaults(func=cmd_generate)

    evaluate = sub.add_parser("evaluate")
    evaluate.add_argument("--truth", required=True)
    evaluate.add_argument("--truth-abundance", required=True)
    evaluate.add_argument("--assignments", required=True)
    evaluate.add_argument("--abundance", required=True)
    evaluate.add_argument("--output", required=True)
    evaluate.set_defaults(func=cmd_evaluate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
