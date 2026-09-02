#!/usr/bin/env python3
"""Context-aware graph path phasing for BridgeAsm GFA.

Four cumulative stages are emitted:
  1. full-read path phasing on the target graph
  2. low-k path projection into the target graph
  3. high-k unitig projection used only across target-graph ambiguities
  4. second-pass read threading using preliminary paths to disambiguate repeats

The optimizer never invents graph edges: every emitted adjacency must already
exist in the target GFA. Low-k/high-k evidence only ranks ambiguous exits.
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

BASE = {65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def canonical(seq: str) -> str:
    seq = seq.upper()
    rev = rc(seq)
    return rev if rev < seq else seq


def rolling_keys(seq: str, k: int) -> Iterator[tuple[int, int]]:
    mask = (1 << (2 * k)) - 1
    key = 0
    valid = 0
    for i, byte in enumerate(seq.encode("ascii", "ignore")):
        value = BASE.get(byte)
        if value is None:
            key = 0
            valid = 0
            continue
        key = ((key << 2) | value) & mask
        valid += 1
        if valid >= k:
            yield i - k + 1, key


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name = None
    chunks: list[str] = []
    with open_text(path) as handle:
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


def read_fastq(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        while True:
            h = handle.readline()
            if not h:
                return
            seq = handle.readline().strip().upper()
            handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ: {path}")
            yield h.strip().split()[0], seq


@dataclass(frozen=True)
class EdgeEvidence:
    direct: int = 0
    gapped: int = 0
    pairs: int = 0


class Graph:
    def __init__(self) -> None:
        self.k = 31
        self.seqs: dict[str, str] = {}
        self.coverage: dict[str, float] = {}
        self.out: dict[str, list[str]] = defaultdict(list)
        self.inc: dict[str, list[str]] = defaultdict(list)
        self.edge: dict[tuple[str, str], EdgeEvidence] = {}
        self.rev: dict[str, str] = {}

    @classmethod
    def from_gfa(cls, path: Path) -> "Graph":
        g = cls()
        overlaps: list[int] = []
        with path.open() as handle:
            for raw in handle:
                fields = raw.rstrip("\n").split("\t")
                if not fields:
                    continue
                if fields[0] == "S" and len(fields) >= 3:
                    uid, seq = fields[1], fields[2].upper()
                    g.seqs[uid] = seq
                    cov = 0.0
                    for tag in fields[3:]:
                        if tag.startswith("KC:f:"):
                            cov = float(tag[5:])
                    g.coverage[uid] = cov
                elif fields[0] == "L" and len(fields) >= 6:
                    src, dst = fields[1], fields[3]
                    ov = fields[5]
                    if ov.endswith("M") and ov[:-1].isdigit():
                        overlaps.append(int(ov[:-1]))
                    tags = {
                        part[:2]: part[5:]
                        for part in fields[6:]
                        if len(part) > 5 and part[2:5] == ":i:"
                    }
                    ev = EdgeEvidence(
                        int(tags.get("DR", 0)),
                        int(tags.get("GR", 0)),
                        int(tags.get("PE", 0)),
                    )
                    g.edge[(src, dst)] = ev
                    g.out[src].append(dst)
                    g.inc[dst].append(src)
        if overlaps:
            g.k = Counter(overlaps).most_common(1)[0][0]
        for uid in g.seqs:
            g.out[uid] = sorted(set(g.out.get(uid, [])))
            g.inc[uid] = sorted(set(g.inc.get(uid, [])))
        seq_index: dict[str, str] = {}
        for uid, seq in g.seqs.items():
            seq_index.setdefault(seq, uid)
        for uid, seq in g.seqs.items():
            g.rev[uid] = seq_index.get(rc(seq), uid)
        return g

    def ambiguous(self, uid: str) -> bool:
        return len(self.out[uid]) > 1 or len(self.inc[uid]) > 1


class KmerIndex:
    """Exact-orientation k-mer -> unique unitig, plus bounded ambiguous candidates."""

    def __init__(self, graph: Graph, k: int, max_ambiguous: int = 8) -> None:
        self.k = k
        self.unique: dict[int, str | None] = {}
        self.ambig: dict[int, list[str]] = {}
        for uid, seq in graph.seqs.items():
            seen: set[int] = set()
            for _, key in rolling_keys(seq, k):
                if key in seen:
                    continue
                seen.add(key)
                if key not in self.unique:
                    self.unique[key] = uid
                    continue
                old = self.unique[key]
                if old == uid:
                    continue
                if old is not None:
                    vals = [old]
                    if uid not in vals:
                        vals.append(uid)
                    self.ambig[key] = vals
                    self.unique[key] = None
                else:
                    vals = self.ambig.get(key)
                    if vals is not None and uid not in vals and len(vals) < max_ambiguous:
                        vals.append(uid)
        self.ambig = {
            key: vals
            for key, vals in self.ambig.items()
            if 1 < len(vals) <= max_ambiguous
        }


def unique_short_bridge(
    graph: Graph, source: str, target: str, max_edges: int = 3
) -> list[str] | None:
    if source == target:
        return [source]
    frontier: list[tuple[str, list[str]]] = [(source, [source])]
    solutions: list[list[str]] = []
    for _ in range(max_edges):
        nxt: list[tuple[str, list[str]]] = []
        for node, path in frontier:
            for child in graph.out[node]:
                if child in path:
                    continue
                path2 = path + [child]
                if child == target:
                    solutions.append(path2)
                    if len(solutions) > 1:
                        return None
                else:
                    nxt.append((child, path2))
        frontier = nxt
        if solutions:
            break
    return solutions[0] if len(solutions) == 1 else None


def preliminary_membership(
    paths: list[list[str]],
) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for pid, path in enumerate(paths):
        for pos, uid in enumerate(path):
            result[uid].append((pid, pos))
    return result


def choose_ambiguous(
    candidates: list[str],
    last_uid: str | None,
    graph: Graph,
    membership: dict[str, list[tuple[int, int]]] | None,
) -> str | None:
    if last_uid is None:
        return None
    viable = [uid for uid in candidates if uid == last_uid or uid in graph.out[last_uid]]
    if len(viable) == 1:
        return viable[0]
    if membership is None:
        return None
    last_places = membership.get(last_uid, [])
    if not last_places:
        return None
    contextual: list[str] = []
    for uid in candidates:
        places = membership.get(uid, [])
        ok = any(
            pid == qid and 0 <= qpos - pos <= 2
            for pid, pos in last_places
            for qid, qpos in places
        )
        if ok:
            contextual.append(uid)
    contextual = sorted(set(contextual))
    return contextual[0] if len(contextual) == 1 else None


def thread_sequence(
    seq: str,
    graph: Graph,
    index: KmerIndex,
    membership: dict[str, list[tuple[int, int]]] | None = None,
    max_bridge_edges: int = 3,
) -> list[list[str]]:
    segments: list[list[str]] = []
    current: list[str] = []
    last_anchor_pos: int | None = None
    for pos, key in rolling_keys(seq, index.k):
        uid = index.unique.get(key)
        if uid is None and membership is not None:
            cands = index.ambig.get(key)
            if cands:
                uid = choose_ambiguous(
                    cands, current[-1] if current else None, graph, membership
                )
        if uid is None:
            continue
        if current and uid == current[-1]:
            last_anchor_pos = pos
            continue
        if not current:
            current = [uid]
            last_anchor_pos = pos
            continue
        last = current[-1]
        if uid in graph.out[last]:
            current.append(uid)
            last_anchor_pos = pos
            continue
        bridge = unique_short_bridge(graph, last, uid, max_bridge_edges)
        if bridge is not None and len(bridge) > 1:
            span = 0 if last_anchor_pos is None else pos - last_anchor_pos
            added_bp = sum(
                max(1, len(graph.seqs[node]) - graph.k) for node in bridge[1:]
            )
            if span <= len(seq) and added_bp <= len(seq) + graph.k:
                current.extend(bridge[1:])
                last_anchor_pos = pos
                continue
        if len(current) >= 2:
            segments.append(current)
        current = [uid]
        last_anchor_pos = pos
    if len(current) >= 2:
        segments.append(current)
    return segments


def add_context(
    counter: Counter[tuple[str, ...]],
    path: list[str],
    max_context: int,
    weight: int = 1,
) -> None:
    seen: set[tuple[str, ...]] = set()
    for n in range(2, min(max_context, len(path)) + 1):
        for i in range(len(path) - n + 1):
            key = tuple(path[i : i + n])
            if key not in seen:
                counter[key] += weight
                seen.add(key)


def collect_read_contexts(
    graph: Graph,
    index: KmerIndex,
    r1: Path,
    r2: Path | None,
    membership: dict[str, list[tuple[int, int]]] | None,
    max_context: int,
) -> tuple[Counter[tuple[str, ...]], dict[str, int]]:
    counter: Counter[tuple[str, ...]] = Counter()
    stats = {"reads": 0, "threaded_reads": 0, "segments": 0, "contexts": 0}
    files = [r1] + ([r2] if r2 else [])
    for path in files:
        assert path is not None
        for _, seq in read_fastq(path):
            stats["reads"] += 1
            segments = thread_sequence(seq, graph, index, membership)
            if segments:
                stats["threaded_reads"] += 1
            stats["segments"] += len(segments)
            for segment in segments:
                add_context(counter, segment, max_context)
    stats["contexts"] = len(counter)
    return counter, stats


def collect_fasta_contexts(
    graph: Graph,
    index: KmerIndex,
    files: list[Path],
    max_context: int,
    only_ambiguous: bool = False,
) -> tuple[Counter[tuple[str, ...]], dict[str, int]]:
    counter: Counter[tuple[str, ...]] = Counter()
    records = segments_n = accepted = 0
    for path in files:
        if not path.exists() or path.stat().st_size == 0:
            continue
        for _, seq in read_fasta(path):
            records += 1
            for segment in thread_sequence(seq, graph, index):
                segments_n += 1
                if only_ambiguous:
                    crosses = any(graph.ambiguous(uid) for uid in segment[1:-1]) or any(
                        len(graph.out[a]) > 1 or len(graph.inc[b]) > 1
                        for a, b in zip(segment, segment[1:])
                    )
                    if not crosses:
                        continue
                accepted += 1
                add_context(counter, segment, max_context)
    return counter, {
        "records": records,
        "segments": segments_n,
        "accepted_segments": accepted,
        "contexts": len(counter),
    }


def context_strength(
    counter: Counter[tuple[str, ...]],
    history: list[str],
    candidate: str,
    forward: bool,
) -> tuple[int, int]:
    best_support = 0
    best_len = 0
    max_hist = min(5, len(history))
    for n in range(1, max_hist + 1):
        if forward:
            key = tuple(history[-n:] + [candidate])
        else:
            key = tuple([candidate] + history[:n])
        support = counter.get(key, 0)
        if support > 0 and (
            n + 1 > best_len or (n + 1 == best_len and support > best_support)
        ):
            best_support = support
            best_len = n + 1
    return best_support, best_len


@dataclass
class Choice:
    uid: str
    score: float
    raw_support: int
    raw_len: int
    proj_support: int
    proj_len: int
    high_support: int
    high_len: int
    second_support: int
    second_len: int
    direct: int
    pairs: int


def candidate_choice(
    graph: Graph,
    history: list[str],
    candidate: str,
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    second_ctx: Counter[tuple[str, ...]],
) -> Choice:
    u, v = (history[-1], candidate) if forward else (candidate, history[0])
    ev = graph.edge.get((u, v), EdgeEvidence())
    rs, rl = context_strength(raw_ctx, history, candidate, forward)
    ps, pl = context_strength(proj_ctx, history, candidate, forward)
    hs, hl = context_strength(high_ctx, history, candidate, forward)
    ss, sl = context_strength(second_ctx, history, candidate, forward)
    score = (
        ev.direct * 2.0
        + ev.gapped * 1.0
        + ev.pairs * 3.0
        + rs * max(2, rl) * 7.0
        + ps * max(2, pl) * 6.0
        + hs * max(2, hl) * 9.0
        + ss * max(2, sl) * 10.0
    )
    return Choice(
        candidate,
        score,
        rs,
        rl,
        ps,
        pl,
        hs,
        hl,
        ss,
        sl,
        ev.direct,
        ev.pairs,
    )


def strong_context(choice: Choice) -> bool:
    return (
        (choice.second_len >= 3 and choice.second_support >= 2)
        or (choice.raw_len >= 3 and choice.raw_support >= 2)
        or (choice.high_len >= 3 and choice.high_support >= 1)
        or (choice.proj_len >= 3 and choice.proj_support >= 1)
    )


def choose_extension(
    graph: Graph,
    history: list[str],
    candidates: list[str],
    used: set[str],
    forward: bool,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    second_ctx: Counter[tuple[str, ...]],
    dominance: float,
    min_direct: int,
) -> Choice | None:
    available = [
        uid
        for uid in candidates
        if uid not in used and graph.rev.get(uid, uid) not in used
    ]
    if not available:
        return None
    if len(available) == 1:
        choice = candidate_choice(
            graph,
            history,
            available[0],
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            second_ctx,
        )
        return choice if choice.score > 0 or len(candidates) == 1 else None
    choices = [
        candidate_choice(
            graph,
            history,
            uid,
            forward,
            raw_ctx,
            proj_ctx,
            high_ctx,
            second_ctx,
        )
        for uid in available
    ]
    choices.sort(
        key=lambda item: (
            -item.score,
            -item.second_len,
            -item.raw_len,
            -item.high_len,
            -item.proj_len,
            item.uid,
        )
    )
    best = choices[0]
    total = sum(max(0.0, choice.score) for choice in choices)
    dominance_fraction = best.score / total if total > 0 else 0.0
    if strong_context(best):
        second = choices[1].score if len(choices) > 1 else 0.0
        if best.score > second * 1.10:
            return best
    if best.direct >= min_direct and dominance_fraction >= dominance:
        return best
    return None


def claim(path: list[str], graph: Graph, used: set[str]) -> None:
    for uid in path:
        used.add(uid)
        used.add(graph.rev.get(uid, uid))


def path_sequence(path: list[str], graph: Graph) -> str:
    if not path:
        return ""
    seq = graph.seqs[path[0]]
    for uid in path[1:]:
        seq += graph.seqs[uid][graph.k :]
    return seq


def resolve_paths(
    graph: Graph,
    raw_ctx: Counter[tuple[str, ...]],
    proj_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    second_ctx: Counter[tuple[str, ...]],
    dominance: float,
    min_direct: int,
    min_length: int,
) -> tuple[list[list[str]], dict[str, int]]:
    used: set[str] = set()
    paths: list[list[str]] = []
    phased_extensions = 0
    branch_stops = 0
    seeds = sorted(
        graph.seqs,
        key=lambda uid: (
            -(graph.coverage.get(uid, 0.0) * len(graph.seqs[uid])),
            -len(graph.seqs[uid]),
            uid,
        ),
    )
    for seed in seeds:
        if seed in used or graph.rev.get(seed, seed) in used:
            continue
        path = [seed]
        local_seen = {seed}
        while True:
            current = path[0]
            choice = choose_extension(
                graph,
                path,
                graph.inc[current],
                used | local_seen,
                False,
                raw_ctx,
                proj_ctx,
                high_ctx,
                second_ctx,
                dominance,
                min_direct,
            )
            if choice is None:
                if len(graph.inc[current]) > 1:
                    branch_stops += 1
                break
            path.insert(0, choice.uid)
            local_seen.add(choice.uid)
            local_seen.add(graph.rev.get(choice.uid, choice.uid))
            if strong_context(choice):
                phased_extensions += 1
        while True:
            current = path[-1]
            choice = choose_extension(
                graph,
                path,
                graph.out[current],
                used | local_seen,
                True,
                raw_ctx,
                proj_ctx,
                high_ctx,
                second_ctx,
                dominance,
                min_direct,
            )
            if choice is None:
                if len(graph.out[current]) > 1:
                    branch_stops += 1
                break
            path.append(choice.uid)
            local_seen.add(choice.uid)
            local_seen.add(graph.rev.get(choice.uid, choice.uid))
            if strong_context(choice):
                phased_extensions += 1
        claim(path, graph, used)
        if len(path_sequence(path, graph)) >= min_length:
            paths.append(path)
    return paths, {
        "paths": len(paths),
        "phased_extensions": phased_extensions,
        "branch_stops": branch_stops,
        "claimed_orientations": len(used),
    }


def write_paths(
    paths: list[list[str]],
    graph: Graph,
    fasta: Path,
    path_tsv: Path,
    min_length: int,
) -> dict[str, int]:
    records: list[tuple[str, list[str]]] = []
    seen: set[str] = set()
    for path in paths:
        seq = path_sequence(path, graph)
        if len(seq) < min_length:
            continue
        can = canonical(seq)
        if can in seen:
            continue
        seen.add(can)
        records.append((can, path))
    records.sort(key=lambda item: (-len(item[0]), item[0]))
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("w") as out:
        for i, (seq, path) in enumerate(records, 1):
            out.write(f">graph_path_{i:07d} len={len(seq)} nodes={len(path)}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j : j + 80] + "\n")
    with path_tsv.open("w") as out:
        out.write("path_id\tlength\tnodes\tunitig_path\n")
        for i, (seq, path) in enumerate(records, 1):
            out.write(
                f"graph_path_{i:07d}\t{len(seq)}\t{len(path)}\t{','.join(path)}\n"
            )
    lengths = sorted((len(seq) for seq, _ in records), reverse=True)
    total = sum(lengths)
    halfway = total / 2
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running >= halfway:
            n50 = length
            break
    return {
        "records": len(records),
        "bases": total,
        "n50": n50,
        "largest": lengths[0] if lengths else 0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", type=Path, required=True)
    ap.add_argument("-1", "--read1", type=Path, required=True)
    ap.add_argument("-2", "--read2", type=Path)
    ap.add_argument("--projection", type=Path, action="append", default=[])
    ap.add_argument("--highk-gfa", type=Path, action="append", default=[])
    ap.add_argument("-o", "--output-dir", type=Path, required=True)
    ap.add_argument("--anchor-k", type=int, default=31)
    ap.add_argument("--max-context", type=int, default=6)
    ap.add_argument("--dominance", type=float, default=0.72)
    ap.add_argument("--min-direct", type=int, default=4)
    ap.add_argument("--min-length", type=int, default=200)
    args = ap.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    graph = Graph.from_gfa(args.gfa)
    if args.anchor_k > graph.k:
        raise SystemExit("anchor-k must be <= target graph k")
    index = KmerIndex(graph, args.anchor_k)
    stats: dict[str, object] = {
        "target_gfa": str(args.gfa),
        "target_k": graph.k,
        "nodes": len(graph.seqs),
        "edges": len(graph.edge),
        "ambiguous_nodes": sum(graph.ambiguous(uid) for uid in graph.seqs),
        "unique_anchor_kmers": sum(
            1 for uid in index.unique.values() if uid is not None
        ),
        "ambiguous_anchor_kmers": len(index.ambig),
    }

    raw_ctx, raw_stats = collect_read_contexts(
        graph, index, args.read1, args.read2, None, args.max_context
    )
    stats["full_read_threading"] = raw_stats
    empty: Counter[tuple[str, ...]] = Counter()

    stage1_paths, stage1_resolve = resolve_paths(
        graph,
        raw_ctx,
        empty,
        empty,
        empty,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    stats["stage1_full_read"] = {
        **stage1_resolve,
        **write_paths(
            stage1_paths,
            graph,
            out / "stage1_full_read.fasta",
            out / "stage1_full_read.paths.tsv",
            args.min_length,
        ),
    }

    proj_ctx, proj_stats = collect_fasta_contexts(
        graph, index, args.projection, args.max_context
    )
    stats["projection_threading"] = proj_stats
    stage2_paths, stage2_resolve = resolve_paths(
        graph,
        raw_ctx,
        proj_ctx,
        empty,
        empty,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    stats["stage2_iterative_projection"] = {
        **stage2_resolve,
        **write_paths(
            stage2_paths,
            graph,
            out / "stage2_iterative_projection.fasta",
            out / "stage2_iterative_projection.paths.tsv",
            args.min_length,
        ),
    }

    high_ctx: Counter[tuple[str, ...]] = Counter()
    high_stats = {"records": 0, "segments": 0, "accepted_segments": 0, "contexts": 0}
    for high_gfa in args.highk_gfa:
        high_graph = Graph.from_gfa(high_gfa)
        tmp = out / (high_gfa.parent.name + ".highk_unitigs.fasta")
        with tmp.open("w") as handle:
            for uid, seq in high_graph.seqs.items():
                handle.write(
                    f">{uid} cov={high_graph.coverage.get(uid, 0):.3f}\n{seq}\n"
                )
        ctx, stage_stats = collect_fasta_contexts(
            graph, index, [tmp], args.max_context, only_ambiguous=True
        )
        high_ctx.update(ctx)
        for key in high_stats:
            high_stats[key] += stage_stats.get(key, 0)
    high_stats["contexts"] = len(high_ctx)
    stats["local_highk_threading"] = high_stats
    stage3_paths, stage3_resolve = resolve_paths(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        empty,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    stats["stage3_local_highk"] = {
        **stage3_resolve,
        **write_paths(
            stage3_paths,
            graph,
            out / "stage3_local_highk.fasta",
            out / "stage3_local_highk.paths.tsv",
            args.min_length,
        ),
    }

    membership = preliminary_membership(stage3_paths)
    second_ctx, second_stats = collect_read_contexts(
        graph,
        index,
        args.read1,
        args.read2,
        membership,
        args.max_context,
    )
    for key in list(second_ctx):
        baseline = raw_ctx.get(key, 0)
        if second_ctx[key] <= baseline:
            del second_ctx[key]
        else:
            second_ctx[key] -= baseline
    second_stats["novel_contexts"] = len(second_ctx)
    second_stats["novel_support"] = sum(second_ctx.values())
    stats["second_pass_threading"] = second_stats
    stage4_paths, stage4_resolve = resolve_paths(
        graph,
        raw_ctx,
        proj_ctx,
        high_ctx,
        second_ctx,
        args.dominance,
        args.min_direct,
        args.min_length,
    )
    stats["stage4_second_pass"] = {
        **stage4_resolve,
        **write_paths(
            stage4_paths,
            graph,
            out / "stage4_second_pass.fasta",
            out / "stage4_second_pass.paths.tsv",
            args.min_length,
        ),
    }

    (out / "graph_optimizer_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
