#!/usr/bin/env python3
"""Run the BridgeBin v2.1 Biological Brain pipeline end to end.

Production path (no truth is used):

  signed v2 core
    -> sparse candidate mining
    -> optional DNA LM / GENERanno / ESM-C feature extraction
    -> feature join
    -> calibrated same-genome pair head
    -> v2.1 split + constrained cross-bin merge + residual rescue

Heavy foundation-model inference remains optional and replaceable. Precomputed feature TSVs
can be supplied instead, which is useful on clusters where model inference is scheduled
separately. The pair model is expected to have been trained with whole-genome hold-out;
its validation-derived split/join thresholds are used automatically unless disabled.

The oracle fixture used by CI is intentionally unsupported here: this is the production
entry point and never consumes truth labels.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--coverage", type=Path)
    parser.add_argument("--markers", type=Path)
    parser.add_argument("--links", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pair-model", type=Path, required=True)

    parser.add_argument("--bridgebin-bin", type=Path, default=ROOT / "target/release/bridgebin")
    parser.add_argument(
        "--bridgebin-v21-bin", type=Path, default=ROOT / "target/release/bridgebin-v21"
    )
    parser.add_argument("--min-contig", type=int, default=1500)
    parser.add_argument("--threads", type=int, default=0, help="reserved for model adapters")

    # Precomputed/cheap biological features.
    parser.add_argument("--dna-embeddings", type=Path)
    parser.add_argument("--gene-hits", type=Path)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--gene-architecture", type=Path)
    parser.add_argument("--protein-embeddings", type=Path)
    parser.add_argument("--protein-repertoire", type=Path)

    # Optional local inference.
    parser.add_argument("--run-dna", action="store_true")
    parser.add_argument("--dna-model", default="multimolecule/dnaberts")
    parser.add_argument("--run-generanno", action="store_true")
    parser.add_argument(
        "--generanno-model", default="GenerTeam/GENERanno-prokaryote-0.5b-cds-annotator"
    )
    parser.add_argument("--run-esmc", action="store_true")
    parser.add_argument("--esmc-model", default="esmc_600m")
    parser.add_argument("--proteins", type=Path, help="protein FASTA if ESM-C is run without GENERanno")
    parser.add_argument("--protein-map", type=Path)
    parser.add_argument("--model-device", default="auto", choices=("auto", "cpu", "cuda"))

    # Repertoire codebook. A fixed codebook is preferred in production.
    parser.add_argument("--protein-codebook", type=Path)
    parser.add_argument(
        "--fit-protein-codebook",
        action="store_true",
        help="learn a sample-local prototype codebook; useful for experiments, not preferred for deployment",
    )
    parser.add_argument("--protein-codebook-clusters", type=int, default=256)

    # Candidate mining.
    parser.add_argument("--candidate-max-pairs", type=int, default=250000)
    parser.add_argument("--anchors-per-bin", type=int, default=10)

    # Calibration/partition behavior.
    parser.add_argument(
        "--ignore-model-thresholds",
        action="store_true",
        help="use bridgebin-v21 defaults instead of validation-derived pair-model thresholds",
    )
    parser.add_argument(
        "--pair-input-confidence",
        type=float,
        default=0.0,
        help="minimum independent input reliability; p_same thresholds already encode decision confidence",
    )
    parser.add_argument("--keep-work", action="store_true")
    return parser.parse_args(argv)


def run(command: Sequence[object], log: Optional[Path] = None) -> None:
    argv = [str(value) for value in command]
    print("+", " ".join(argv), flush=True)
    if log is None:
        subprocess.run(argv, check=True)
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT)


def existing(path: Optional[Path], label: str) -> Optional[Path]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path.resolve()


def append_optional(command: List[object], flag: str, path: Optional[Path]) -> None:
    if path is not None:
        command.extend([flag, path])


def read_model_thresholds(path: Path) -> Dict[str, float]:
    model = json.loads(path.read_text(encoding="utf-8"))
    training = model.get("training") or {}
    result: Dict[str, float] = {}
    for output, key in (
        ("join", "recommended_join_min_same"),
        ("split", "recommended_split_max_same"),
    ):
        value = training.get(key)
        if value is None:
            continue
        value = float(value)
        if 0.0 <= value <= 1.0:
            result[output] = value
    result["model_version"] = float(model.get("version", 0))
    return result


def build_core(args: argparse.Namespace, core_dir: Path, work: Path) -> Path:
    command: List[object] = [
        args.bridgebin_bin,
        "--contigs",
        args.contigs,
        "--out-dir",
        core_dir,
        "--algorithm",
        "v2",
        "--min-contig",
        args.min_contig,
        "--v2-core-attraction",
        "0.74",
        "--v2-component-attraction",
        "0.66",
        "--v2-rescue-attraction",
        "0.68",
        "--v2-rescue-margin",
        "0.02",
    ]
    append_optional(command, "--coverage", args.coverage)
    append_optional(command, "--markers", args.markers)
    append_optional(command, "--links", args.links)
    run(command, work / "01_core.log")
    assignments = core_dir / "assignments.tsv"
    if not assignments.exists():
        raise RuntimeError("v2 core did not produce assignments.tsv")
    return assignments


def mine_candidates(args: argparse.Namespace, assignments: Path, work: Path) -> Path:
    output = work / "candidate_pairs.tsv"
    command: List[object] = [
        sys.executable,
        SCRIPTS / "bridgebin_candidate_pairs.py",
        "--contigs",
        args.contigs,
        "--assignments",
        assignments,
        "--output",
        output,
        "--min-length",
        args.min_contig,
        "--anchors-per-bin",
        args.anchors_per_bin,
        "--max-pairs",
        args.candidate_max_pairs,
    ]
    append_optional(command, "--coverage", args.coverage)
    run(command, work / "02_candidates.log")
    return output


def prepare_features(args: argparse.Namespace, work: Path) -> Path:
    dna = existing(args.dna_embeddings, "DNA embeddings")
    architecture = existing(args.gene_architecture, "gene architecture")
    protein_embeddings = existing(args.protein_embeddings, "protein embeddings")
    repertoire = existing(args.protein_repertoire, "protein repertoire")
    gene_hits = existing(args.gene_hits, "gene hits")
    taxonomy = existing(args.taxonomy, "taxonomy")

    if args.run_dna:
        if dna is not None:
            raise ValueError("choose either --dna-embeddings or --run-dna")
        dna = work / "dna_embeddings.tsv"
        run(
            [
                sys.executable,
                SCRIPTS / "bridgebin_dna_embed.py",
                "--contigs",
                args.contigs,
                "--output",
                dna,
                "--model",
                args.dna_model,
                "--device",
                args.model_device,
            ],
            work / "03_dna.log",
        )

    generated_proteins: Optional[Path] = None
    generated_map: Optional[Path] = None
    if args.run_generanno:
        if architecture is not None:
            raise ValueError("choose either --gene-architecture or --run-generanno")
        prefix = work / "generanno"
        run(
            [
                sys.executable,
                SCRIPTS / "bridgebin_generanno_cds.py",
                "--contigs",
                args.contigs,
                "--output-prefix",
                prefix,
                "--model",
                args.generanno_model,
                "--device",
                args.model_device,
            ],
            work / "04_generanno.log",
        )
        architecture = Path(str(prefix) + ".architecture.tsv")
        generated_proteins = Path(str(prefix) + ".proteins.faa")
        generated_map = Path(str(prefix) + ".protein_map.tsv")

    if args.run_esmc:
        if protein_embeddings is not None:
            raise ValueError("choose either --protein-embeddings or --run-esmc")
        proteins = existing(args.proteins, "protein FASTA") or generated_proteins
        mapping = existing(args.protein_map, "protein mapping") or generated_map
        if proteins is None:
            raise ValueError("--run-esmc requires --proteins or --run-generanno")
        protein_embeddings = work / "esmc_embeddings.tsv"
        command = [
            sys.executable,
            SCRIPTS / "bridgebin_esmc_embed.py",
            "--proteins",
            proteins,
            "--output",
            protein_embeddings,
            "--model",
            args.esmc_model,
            "--device",
            args.model_device,
        ]
        append_optional(command, "--mapping", mapping)
        run(command, work / "05_esmc.log")

    if repertoire is None and protein_embeddings is not None and (
        args.protein_codebook is not None or args.fit_protein_codebook
    ):
        repertoire = work / "protein_repertoire.tsv"
        if args.protein_codebook is not None:
            codebook = existing(args.protein_codebook, "protein codebook")
            run(
                [
                    sys.executable,
                    SCRIPTS / "bridgebin_protein_repertoire.py",
                    "transform",
                    "--embeddings",
                    protein_embeddings,
                    "--codebook",
                    codebook,
                    "--output",
                    repertoire,
                ],
                work / "06_repertoire.log",
            )
        else:
            codebook = work / "protein_codebook.npz"
            run(
                [
                    sys.executable,
                    SCRIPTS / "bridgebin_protein_repertoire.py",
                    "fit-transform",
                    "--embeddings",
                    protein_embeddings,
                    "--codebook-out",
                    codebook,
                    "--output",
                    repertoire,
                    "--clusters",
                    args.protein_codebook_clusters,
                ],
                work / "06_repertoire.log",
            )

    inputs: List[Path] = []
    if any(path is not None for path in (dna, gene_hits, protein_embeddings, taxonomy)):
        base = work / "base_bio_features.tsv"
        command: List[object] = [
            sys.executable,
            SCRIPTS / "bridgebin_bio_features.py",
            "--output",
            base,
        ]
        append_optional(command, "--dna-embeddings", dna)
        append_optional(command, "--gene-hits", gene_hits)
        append_optional(command, "--protein-embeddings", protein_embeddings)
        append_optional(command, "--taxonomy", taxonomy)
        run(command, work / "07_base_features.log")
        inputs.append(base)
    if architecture is not None:
        inputs.append(architecture)
    if repertoire is not None:
        inputs.append(repertoire)
    if not inputs:
        raise ValueError(
            "no Biological Brain features available; provide precomputed features or enable model inference"
        )

    joined = work / "biobrain_features.tsv"
    command = [sys.executable, SCRIPTS / "bridgebin_join_features.py"]
    for path in inputs:
        command.extend(["--input", path])
    command.extend(["--output", joined])
    run(command, work / "08_join_features.log")
    return joined


def score_pairs(
    args: argparse.Namespace, features: Path, candidates: Path, work: Path
) -> Path:
    output = work / "pair_scores.tsv"
    run(
        [
            sys.executable,
            SCRIPTS / "bridgebin_pair_head.py",
            "score",
            "--features",
            features,
            "--pairs",
            candidates,
            "--model",
            args.pair_model,
            "--output",
            output,
        ],
        work / "09_pair_head.log",
    )
    return output


def run_final(
    args: argparse.Namespace,
    features: Path,
    pair_scores: Path,
    final_dir: Path,
    work: Path,
) -> Dict[str, float]:
    thresholds = {} if args.ignore_model_thresholds else read_model_thresholds(args.pair_model)
    command: List[object] = [
        args.bridgebin_v21_bin,
        "--contigs",
        args.contigs,
        "--pair-scores",
        pair_scores,
        "--bio-features",
        features,
        "--out-dir",
        final_dir,
        "--min-contig",
        args.min_contig,
        "--v21-min-confidence",
        args.pair_input_confidence,
    ]
    append_optional(command, "--coverage", args.coverage)
    append_optional(command, "--markers", args.markers)
    append_optional(command, "--links", args.links)
    if "split" in thresholds and "join" in thresholds and thresholds["split"] < thresholds["join"]:
        join = thresholds["join"]
        split = thresholds["split"]
        command.extend(
            [
                "--v21-split-max-same",
                f"{split:.8f}",
                "--v21-join-min-same",
                f"{join:.8f}",
                "--v21-rescue-min-same",
                f"{join:.8f}",
                "--v21-bin-merge-min-same",
                f"{max(join, 0.92):.8f}",
            ]
        )
    run(command, work / "10_final.log")
    return thresholds


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.contigs = existing(args.contigs, "contigs")
    args.coverage = existing(args.coverage, "coverage")
    args.markers = existing(args.markers, "markers")
    args.links = existing(args.links, "links")
    args.pair_model = existing(args.pair_model, "pair model")
    args.bridgebin_bin = existing(args.bridgebin_bin, "bridgebin binary")
    args.bridgebin_v21_bin = existing(args.bridgebin_v21_bin, "bridgebin-v21 binary")
    if args.min_contig < 1 or args.candidate_max_pairs < 1 or args.anchors_per_bin < 1:
        raise SystemExit("length/pair/anchor limits must be positive")
    if not 0.0 <= args.pair_input_confidence <= 1.0:
        raise SystemExit("--pair-input-confidence must be in [0,1]")

    out = args.out_dir.resolve()
    work = out / "work"
    core = out / "core"
    final = out / "final"
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    assignments = build_core(args, core, work)
    candidates = mine_candidates(args, assignments, work)
    features = prepare_features(args, work)
    pair_scores = score_pairs(args, features, candidates, work)
    thresholds = run_final(args, features, pair_scores, final, work)

    manifest = {
        "contigs": str(args.contigs),
        "coverage": None if args.coverage is None else str(args.coverage),
        "markers": None if args.markers is None else str(args.markers),
        "links": None if args.links is None else str(args.links),
        "pair_model": str(args.pair_model),
        "pair_model_thresholds": thresholds,
        "core_assignments": str(assignments),
        "candidate_pairs": str(candidates),
        "biobrain_features": str(features),
        "pair_scores": str(pair_scores),
        "final_dir": str(final),
        "dna_model": args.dna_model if args.run_dna else None,
        "generanno_model": args.generanno_model if args.run_generanno else None,
        "esmc_model": args.esmc_model if args.run_esmc else None,
        "production_truth_used": False,
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"bridgebin-biobrain: final={final} manifest={out / 'manifest.json'}")

    if not args.keep_work:
        # Keep the final model-facing tables for audit/reproducibility; only remove logs
        # and bulky transient model products that can be regenerated from the manifest.
        for path in work.glob("*.log"):
            path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
