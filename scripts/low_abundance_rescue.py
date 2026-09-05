#!/usr/bin/env python3
"""Taxon-agnostic low-abundance recovery on top of the promoted Stage-8 backbone.

The rescue has two independent evidence channels:

1. graph-safe rare paths: maximal non-branching paths and isolated nodes from
   the low-k/k31 graphs. We never cross an ambiguous junction. Candidates are
   required to carry substantial sequence absent from Stage 8, have coherent
   graph coverage, and are quota-balanced by abundance/GC bins so the dominant
   taxon cannot consume the rescue budget.
2. residual-read reassembly: paired reads poorly represented by Stage 8 are
   selected and reassembled at k=21. Only highly novel resulting contigs are
   retained.

Three candidates are emitted for measurement: graph_strict, residual_strict,
and hybrid_balanced. Stage 8 is never replaced; rescue sequence is only added,
deduplicated, and optionally joined through unique exact overlaps.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import resource
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def canonical(seq: str) -> str:
    seq = seq.upper()
    rev = rc(seq)
    return min(seq, rev)


def canonical_kmer(seq: str) -> str:
    return canonical(seq)


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name = None
    chunks: list[str] = []
    with path.open() as handle:
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
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_fasta(records: Iterable[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for pos in range(0, len(seq), 80):
                handle.write(seq[pos : pos + 80] + "\n")


def kmers(seq: str, k: int, stride: int = 1) -> Iterator[str]:
    if len(seq) < k:
        return
    for pos in range(0, len(seq) - k + 1, stride):
        kmer = seq[pos : pos + k]
        if "N" not in kmer:
            yield canonical_kmer(kmer)


def backbone_kmers(path: Path, k: int) -> tuple[set[str], int]:
    result: set[str] = set()
    bases = 0
    for _name, seq in fasta_records(path):
        bases += len(seq)
        result.update(kmers(seq, k))
    return result, bases


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for item in fields:
        parts = item.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    direct: int
    gapped: int
    pairs: int

    @property
    def physical(self) -> int:
        return max(self.direct, self.gapped, self.pairs)


@dataclass
class GraphData:
    name: str
    k: int
    seqs: dict[str, str]
    coverage: dict[str, float]
    out: dict[str, list[str]]
    inc: dict[str, list[str]]
    edges: dict[tuple[str, str], Edge]


def parse_gfa(path: Path) -> GraphData:
    seqs: dict[str, str] = {}
    coverage: dict[str, float] = {}
    out: dict[str, list[str]] = defaultdict(list)
    inc: dict[str, list[str]] = defaultdict(list)
    edges: dict[tuple[str, str], Edge] = {}
    overlaps: list[int] = []
    with path.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S" and len(fields) >= 3:
                uid = fields[1]
                seqs[uid] = fields[2].upper()
                tags = parse_tags(fields[3:])
                try:
                    coverage[uid] = float(tags.get("KC", "0"))
                except ValueError:
                    coverage[uid] = 0.0
            elif (
                fields[0] == "L"
                and len(fields) >= 6
                and fields[2] == "+"
                and fields[4] == "+"
            ):
                src, dst = fields[1], fields[3]
                ov = fields[5]
                if ov.endswith("M") and ov[:-1].isdigit():
                    overlaps.append(int(ov[:-1]))
                tags = parse_tags(fields[6:])
                edge = Edge(
                    src,
                    dst,
                    int(tags.get("DR", "0")),
                    int(tags.get("GR", "0")),
                    int(tags.get("PE", "0")),
                )
                edges[(src, dst)] = edge
                out[src].append(dst)
                inc[dst].append(src)
    for uid in seqs:
        out[uid] = sorted(set(out.get(uid, [])))
        inc[uid] = sorted(set(inc.get(uid, [])))
    k = statistics.mode(overlaps) if overlaps else 0
    return GraphData(path.parent.name + "/" + path.name, k, seqs, coverage, out, inc, edges)


def assemble_path(graph: GraphData, nodes: list[str]) -> str:
    if not nodes:
        return ""
    seq = graph.seqs[nodes[0]]
    overlap = graph.k
    for uid in nodes[1:]:
        seq += graph.seqs[uid][min(overlap, len(graph.seqs[uid])) :]
    return seq


def maximal_nonbranching_paths(graph: GraphData) -> list[list[str]]:
    """Return maximal paths in the graph after removing ambiguous nodes.

    A safe node has at most one predecessor and at most one successor in the
    original graph. Branch nodes themselves are excluded; edges are retained
    only when the source has one outgoing edge and the destination has one
    incoming edge. This preserves branch-adjacent terminal tips as singleton
    candidates without ever crossing the ambiguous junction.
    """

    safe_nodes = {
        uid
        for uid in graph.seqs
        if len(graph.inc[uid]) <= 1 and len(graph.out[uid]) <= 1
    }
    safe_out: dict[str, list[str]] = {uid: [] for uid in safe_nodes}
    safe_inc: dict[str, list[str]] = {uid: [] for uid in safe_nodes}
    for (src, dst), _edge in graph.edges.items():
        if src not in safe_nodes or dst not in safe_nodes:
            continue
        if len(graph.out[src]) != 1 or len(graph.inc[dst]) != 1:
            continue
        safe_out[src].append(dst)
        safe_inc[dst].append(src)

    paths: list[list[str]] = []
    visited: set[str] = set()
    starts = sorted(uid for uid in safe_nodes if not safe_inc[uid])
    for start in starts:
        if start in visited:
            continue
        path = [start]
        visited.add(start)
        cur = start
        while len(safe_out[cur]) == 1:
            nxt = safe_out[cur][0]
            if nxt in visited or len(safe_inc[nxt]) != 1:
                break
            path.append(nxt)
            visited.add(nxt)
            cur = nxt
        paths.append(path)

    for start in sorted(safe_nodes - visited):
        if start in visited:
            continue
        path = [start]
        visited.add(start)
        cur = start
        while len(safe_out[cur]) == 1:
            nxt = safe_out[cur][0]
            if nxt in visited:
                break
            path.append(nxt)
            visited.add(nxt)
            cur = nxt
        paths.append(path)
    return paths


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    vals = sorted(values)
    index = (len(vals) - 1) * q
    lo = math.floor(index)
    hi = math.ceil(index)
    if lo == hi:
        return vals[lo]
    frac = index - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


@dataclass
class GraphCandidate:
    source: str
    nodes: list[str]
    seq: str
    median_coverage: float
    min_coverage: float
    max_coverage: float
    coverage_ratio: float
    novel_kmers: int
    novel_fraction: float
    abundance_band: int = 0
    gc_bin: int = 0
    score: float = 0.0


def graph_candidates(
    graphs: list[GraphData], backbone: set[str], novel_k: int, min_length: int
) -> list[GraphCandidate]:
    candidates: list[GraphCandidate] = []
    for graph in graphs:
        local: list[GraphCandidate] = []
        for nodes in maximal_nonbranching_paths(graph):
            seq = canonical(assemble_path(graph, nodes))
            if len(seq) < min_length:
                continue
            path_kmers = set(kmers(seq, novel_k))
            if not path_kmers:
                continue
            novel = path_kmers - backbone
            coverages = [max(0.0, graph.coverage.get(uid, 0.0)) for uid in nodes]
            positive = [x for x in coverages if x > 0]
            med = statistics.median(positive) if positive else 0.0
            mn = min(positive) if positive else 0.0
            mx = max(positive) if positive else 0.0
            ratio = mx / max(mn, 1e-6) if positive else float("inf")
            local.append(
                GraphCandidate(
                    source=graph.name,
                    nodes=nodes,
                    seq=seq,
                    median_coverage=med,
                    min_coverage=mn,
                    max_coverage=mx,
                    coverage_ratio=ratio,
                    novel_kmers=len(novel),
                    novel_fraction=len(novel) / len(path_kmers),
                    gc_bin=min(9, int(gc_fraction(seq) * 10.0)),
                )
            )
        covs = [c.median_coverage for c in local if c.median_coverage > 0]
        q25, q50, q75 = quantile(covs, 0.25), quantile(covs, 0.50), quantile(covs, 0.75)
        for cand in local:
            cov = cand.median_coverage
            if cov <= q25:
                cand.abundance_band = 0
            elif cov <= q50:
                cand.abundance_band = 1
            elif cov <= q75:
                cand.abundance_band = 2
            else:
                cand.abundance_band = 3
            abundance_bonus = (1.25, 1.10, 1.00, 0.90)[cand.abundance_band]
            coherence = 1.0 / max(1.0, math.log2(max(2.0, cand.coverage_ratio)))
            cand.score = (
                cand.novel_kmers
                * cand.novel_fraction
                * abundance_bonus
                * (0.5 + 0.5 * coherence)
                * math.log2(2.0 + max(0.0, cand.median_coverage))
            )
        candidates.extend(local)
    return candidates


def select_graph_candidates(
    candidates: list[GraphCandidate],
    backbone: set[str],
    *,
    min_novel_kmers: int,
    min_novel_fraction: float,
    max_coverage_ratio: float,
    max_total_bases: int,
    per_cluster_fraction: float,
    novel_k: int,
) -> list[GraphCandidate]:
    eligible = [
        c
        for c in candidates
        if c.novel_kmers >= min_novel_kmers
        and c.novel_fraction >= min_novel_fraction
        and c.coverage_ratio <= max_coverage_ratio
        and c.median_coverage > 0
    ]
    eligible.sort(key=lambda c: (-c.score, -c.novel_kmers, -len(c.seq), c.seq))
    selected: list[GraphCandidate] = []
    selected_novel: set[str] = set()
    seen_seq: set[str] = set()
    cluster_bases: defaultdict[tuple[int, int], int] = defaultdict(int)
    total = 0
    cluster_cap = max(20_000, int(max_total_bases * per_cluster_fraction))
    for cand in eligible:
        if cand.seq in seen_seq:
            continue
        novel = set(kmers(cand.seq, novel_k)) - backbone
        fresh = novel - selected_novel
        if len(fresh) < min_novel_kmers:
            continue
        if len(fresh) < max(min_novel_kmers, int(0.5 * len(novel))):
            continue
        key = (cand.abundance_band, cand.gc_bin)
        if cluster_bases[key] + len(cand.seq) > cluster_cap:
            continue
        if total + len(cand.seq) > max_total_bases:
            continue
        selected.append(cand)
        seen_seq.add(cand.seq)
        selected_novel.update(novel)
        cluster_bases[key] += len(cand.seq)
        total += len(cand.seq)
    return selected


def read_fastq(handle) -> tuple[str, str, str, str] | None:
    h = handle.readline()
    if not h:
        return None
    s = handle.readline()
    plus = handle.readline()
    q = handle.readline()
    if not q:
        raise ValueError("truncated FASTQ")
    return h, s, plus, q


def represented_fraction(seq: str, represented: set[str], k: int, stride: int) -> float:
    vals = list(kmers(seq.strip().upper(), k, stride))
    if not vals:
        return 1.0
    return sum(q in represented for q in vals) / len(vals)


def select_residual_pairs(
    read1: Path,
    read2: Path,
    out1: Path,
    out2: Path,
    represented: set[str],
    k: int,
    stride: int,
    strong_novel: float,
    weak_novel: float,
) -> dict[str, int | float]:
    total = 0
    kept = 0
    out1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(read1, "rt") as a, gzip.open(read2, "rt") as b, gzip.open(
        out1, "wt"
    ) as oa, gzip.open(out2, "wt") as ob:
        while True:
            r1 = read_fastq(a)
            r2 = read_fastq(b)
            if r1 is None or r2 is None:
                if r1 is not None or r2 is not None:
                    raise ValueError("paired FASTQ files have different record counts")
                break
            total += 1
            f1 = represented_fraction(r1[1], represented, k, stride)
            f2 = represented_fraction(r2[1], represented, k, stride)
            keep = max(f1, f2) <= weak_novel or (
                min(f1, f2) <= strong_novel and max(f1, f2) <= 0.60
            )
            if not keep:
                continue
            kept += 1
            oa.writelines(r1)
            ob.writelines(r2)
    return {
        "pairs_total": total,
        "pairs_kept": kept,
        "kept_fraction": kept / max(1, total),
    }


@dataclass
class FastaCandidate:
    source: str
    name: str
    seq: str
    novel_kmers: int
    novel_fraction: float
    score: float


def novel_fasta_candidates(
    inputs: list[Path], backbone: set[str], novel_k: int, min_length: int
) -> list[FastaCandidate]:
    out: list[FastaCandidate] = []
    for path in inputs:
        if not path.exists():
            continue
        for name, raw_seq in fasta_records(path):
            seq = canonical(raw_seq)
            if len(seq) < min_length:
                continue
            all_kmers = set(kmers(seq, novel_k))
            if not all_kmers:
                continue
            novel = all_kmers - backbone
            frac = len(novel) / len(all_kmers)
            out.append(
                FastaCandidate(
                    source=str(path),
                    name=name,
                    seq=seq,
                    novel_kmers=len(novel),
                    novel_fraction=frac,
                    score=len(novel) * frac * math.log2(max(2, len(seq))),
                )
            )
    return out


def select_fasta_candidates(
    candidates: list[FastaCandidate],
    backbone: set[str],
    *,
    min_novel_kmers: int,
    min_novel_fraction: float,
    max_total_bases: int,
    novel_k: int,
) -> list[FastaCandidate]:
    eligible = [
        c
        for c in candidates
        if c.novel_kmers >= min_novel_kmers and c.novel_fraction >= min_novel_fraction
    ]
    eligible.sort(key=lambda c: (-c.score, -c.novel_kmers, -len(c.seq), c.seq))
    selected: list[FastaCandidate] = []
    seen: set[str] = set()
    selected_novel: set[str] = set()
    total = 0
    for cand in eligible:
        if cand.seq in seen:
            continue
        novel = set(kmers(cand.seq, novel_k)) - backbone
        fresh = novel - selected_novel
        if len(fresh) < max(min_novel_kmers, int(0.60 * len(novel))):
            continue
        if total + len(cand.seq) > max_total_bases:
            continue
        selected.append(cand)
        seen.add(cand.seq)
        selected_novel.update(novel)
        total += len(cand.seq)
    return selected


def run(cmd: list[object]) -> float:
    print("+", " ".join(map(str, cmd)), flush=True)
    started = time.monotonic()
    subprocess.run(list(map(str, cmd)), check=True)
    return time.monotonic() - started


def make_union_candidate(
    scripts: Path,
    backbone: Path,
    additions: list[Path],
    outdir: Path,
    timings: dict[str, float],
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    union = outdir / "union.fasta"
    noncontained = outdir / "noncontained.fasta"
    final = outdir / "primary_contigs.fasta"
    timings[f"merge_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "merge_fasta_unique.py",
            union,
            backbone,
            *additions,
            "--min-length",
            200,
        ]
    )
    timings[f"contain_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "filter_contained_fasta.py",
            union,
            noncontained,
            "--min-length",
            200,
            "--seed-k",
            21,
            "--window",
            12,
            "--candidate-minimizers",
            16,
            "--removed-tsv",
            outdir / "contained_removed.tsv",
            "--stats-json",
            outdir / "containment_stats.json",
        ]
    )
    timings[f"stitch_{outdir.name}"] = run(
        [
            sys.executable,
            scripts / "stitch_exact_overlaps.py",
            final,
            noncontained,
            "--min-overlap",
            31,
            "--overlap-margin",
            10,
            "--seed-length",
            31,
            "--max-seed-occurrences",
            32,
            "--min-length",
            200,
        ]
    )
    return final


def graph_metadata(selected: list[GraphCandidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "id\tsource\tnodes\tlength\tmedian_coverage\tmin_coverage\tmax_coverage"
            "\tcoverage_ratio\tnovel_kmers\tnovel_fraction\tabundance_band\tgc_bin\tscore\n"
        )
        for i, cand in enumerate(selected, 1):
            handle.write(
                f"{i}\t{cand.source}\t{len(cand.nodes)}\t{len(cand.seq)}"
                f"\t{cand.median_coverage:.4f}\t{cand.min_coverage:.4f}"
                f"\t{cand.max_coverage:.4f}\t{cand.coverage_ratio:.4f}"
                f"\t{cand.novel_kmers}\t{cand.novel_fraction:.6f}"
                f"\t{cand.abundance_band}\t{cand.gc_bin}\t{cand.score:.6f}\n"
            )


def fasta_metadata(selected: list[FastaCandidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("id\tsource\tname\tlength\tnovel_kmers\tnovel_fraction\tscore\n")
        for i, cand in enumerate(selected, 1):
            handle.write(
                f"{i}\t{cand.source}\t{cand.name}\t{len(cand.seq)}\t{cand.novel_kmers}"
                f"\t{cand.novel_fraction:.6f}\t{cand.score:.6f}\n"
            )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridgeasm", type=Path, required=True)
    ap.add_argument("--pipeline-dir", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--novel-k", type=int, default=31)
    ap.add_argument("--pair-stride", type=int, default=3)
    ap.add_argument("--pair-strong-novel", type=float, default=0.10)
    ap.add_argument("--pair-weak-novel", type=float, default=0.35)
    ap.add_argument("--graph-strict-fraction", type=float, default=0.08)
    ap.add_argument("--graph-balanced-fraction", type=float, default=0.15)
    ap.add_argument("--residual-strict-fraction", type=float, default=0.10)
    ap.add_argument("--residual-balanced-fraction", type=float, default=0.15)
    args = ap.parse_args()

    started = time.monotonic()
    timings: dict[str, float] = {}
    scripts = Path(__file__).resolve().parent
    out = args.pipeline_dir
    backbone = out / "bridge_backbone.fasta"
    base = out / "current_pipeline" / "iterative"
    k21_gfa = base / "k21_recall" / "assembly.gfa"
    k31_gfa = base / "k31_resolve" / "assembly.gfa"
    required = [backbone, k21_gfa, k31_gfa, args.read1, args.read2]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("missing low-abundance rescue inputs: " + ", ".join(missing))

    represented, backbone_bases = backbone_kmers(backbone, args.novel_k)
    rescue_dir = out / "low_abundance_rescue"
    rescue_dir.mkdir(parents=True, exist_ok=True)

    graphs = [parse_gfa(k21_gfa), parse_gfa(k31_gfa)]
    all_graph = graph_candidates(graphs, represented, args.novel_k, 200)
    graph_strict = select_graph_candidates(
        all_graph,
        represented,
        min_novel_kmers=64,
        min_novel_fraction=0.65,
        max_coverage_ratio=4.0,
        max_total_bases=max(20_000, int(backbone_bases * args.graph_strict_fraction)),
        per_cluster_fraction=0.30,
        novel_k=args.novel_k,
    )
    graph_balanced = select_graph_candidates(
        all_graph,
        represented,
        min_novel_kmers=40,
        min_novel_fraction=0.45,
        max_coverage_ratio=6.0,
        max_total_bases=max(40_000, int(backbone_bases * args.graph_balanced_fraction)),
        per_cluster_fraction=0.35,
        novel_k=args.novel_k,
    )
    graph_strict_fasta = rescue_dir / "graph_strict_additions.fasta"
    graph_balanced_fasta = rescue_dir / "graph_balanced_additions.fasta"
    write_fasta(
        ((f"graph_strict_{i:06d}", c.seq) for i, c in enumerate(graph_strict, 1)),
        graph_strict_fasta,
    )
    write_fasta(
        ((f"graph_balanced_{i:06d}", c.seq) for i, c in enumerate(graph_balanced, 1)),
        graph_balanced_fasta,
    )
    graph_metadata(graph_strict, rescue_dir / "graph_strict.tsv")
    graph_metadata(graph_balanced, rescue_dir / "graph_balanced.tsv")

    rare_r1 = rescue_dir / "rare_R1.fastq.gz"
    rare_r2 = rescue_dir / "rare_R2.fastq.gz"
    pair_stats = select_residual_pairs(
        args.read1,
        args.read2,
        rare_r1,
        rare_r2,
        represented,
        args.novel_k,
        args.pair_stride,
        args.pair_strong_novel,
        args.pair_weak_novel,
    )

    residual_asm = rescue_dir / "residual_k21"
    if pair_stats["pairs_kept"] >= 100:
        timings["residual_k21_assembly"] = run(
            [
                args.bridgeasm,
                "assemble",
                "-1",
                rare_r1,
                "-2",
                rare_r2,
                "-o",
                residual_asm,
                "-k",
                21,
                "--min-count",
                2,
                "--mercy-max-kmers",
                48,
                "--mercy-min-support",
                1,
                "--mercy-min-quality",
                30,
                "--min-read-support",
                2,
                "--min-pair-support",
                2,
                "--min-primary-support",
                3,
                "--primary-dominance",
                0.80,
                "--threaded-path-cover",
                "--major-path-cover",
                "--path-cover-secondary-dominance",
                0.20,
                "--min-contig-length",
                200,
                "--threads",
                args.threads,
            ]
        )
        residual_inputs = [
            residual_asm / "primary_contigs.fasta",
            residual_asm / "haplotigs.fasta",
        ]
    else:
        residual_asm.mkdir(parents=True, exist_ok=True)
        residual_inputs = []

    all_residual = novel_fasta_candidates(
        residual_inputs, represented, args.novel_k, 200
    )
    residual_strict = select_fasta_candidates(
        all_residual,
        represented,
        min_novel_kmers=64,
        min_novel_fraction=0.70,
        max_total_bases=max(20_000, int(backbone_bases * args.residual_strict_fraction)),
        novel_k=args.novel_k,
    )
    residual_balanced = select_fasta_candidates(
        all_residual,
        represented,
        min_novel_kmers=40,
        min_novel_fraction=0.50,
        max_total_bases=max(40_000, int(backbone_bases * args.residual_balanced_fraction)),
        novel_k=args.novel_k,
    )
    residual_strict_fasta = rescue_dir / "residual_strict_additions.fasta"
    residual_balanced_fasta = rescue_dir / "residual_balanced_additions.fasta"
    write_fasta(
        ((f"residual_strict_{i:06d}", c.seq) for i, c in enumerate(residual_strict, 1)),
        residual_strict_fasta,
    )
    write_fasta(
        ((f"residual_balanced_{i:06d}", c.seq) for i, c in enumerate(residual_balanced, 1)),
        residual_balanced_fasta,
    )
    fasta_metadata(residual_strict, rescue_dir / "residual_strict.tsv")
    fasta_metadata(residual_balanced, rescue_dir / "residual_balanced.tsv")

    graph_final = make_union_candidate(
        scripts,
        backbone,
        [graph_strict_fasta],
        rescue_dir / "candidate_graph_strict",
        timings,
    )
    residual_final = make_union_candidate(
        scripts,
        backbone,
        [residual_strict_fasta],
        rescue_dir / "candidate_residual_strict",
        timings,
    )
    hybrid_final = make_union_candidate(
        scripts,
        backbone,
        [graph_balanced_fasta, residual_balanced_fasta],
        rescue_dir / "candidate_hybrid_balanced",
        timings,
    )

    outputs = {
        "stage8_backbone": str(backbone),
        "graph_strict": str(graph_final),
        "residual_strict": str(residual_final),
        "hybrid_balanced": str(hybrid_final),
    }
    stats = {
        "pipeline": "bridge-low-abundance-rescue-v1",
        "policy": {
            "ambiguous_graph_junctions_crossed": False,
            "production_backbone_replaced": False,
            "continuity_operation": "unique_exact_overlap_only",
            "abundance_balancing": "coverage_quartile_x_gc_decile_quota",
        },
        "backbone_bases": backbone_bases,
        "backbone_novel_kmers": len(represented),
        "graph": {
            "raw_candidates": len(all_graph),
            "strict_selected": len(graph_strict),
            "strict_bases": sum(len(c.seq) for c in graph_strict),
            "balanced_selected": len(graph_balanced),
            "balanced_bases": sum(len(c.seq) for c in graph_balanced),
        },
        "residual_pairs": pair_stats,
        "residual_assembly": {
            "raw_candidates": len(all_residual),
            "strict_selected": len(residual_strict),
            "strict_bases": sum(len(c.seq) for c in residual_strict),
            "balanced_selected": len(residual_balanced),
            "balanced_bases": sum(len(c.seq) for c in residual_balanced),
        },
        "timings_seconds": timings,
        "wall_seconds": time.monotonic() - started,
        "peak_child_rss_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        / 1024.0,
        "outputs": outputs,
    }
    (rescue_dir / "low_abundance_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    manifest_path = out / "pipeline_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["low_abundance_rescue"] = stats
    manifest.setdefault("outputs", {}).update(
        {f"low_abundance_{k}": v for k, v in outputs.items() if k != "stage8_backbone"}
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
