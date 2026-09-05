#!/usr/bin/env python3
"""Stage33: fragment-aware local microassembly from immutable Stage24 seed ends.

Stage31 whole-contig bridges and Stage32 exact GFA endpoint walks were both safe
but too sparse.  The remaining first-principles failure mode is sequence removed
from the global k31 graph by the min_count>=2 gate.  This stage therefore works
locally around each immutable seed end:

  * terminal k21 signatures are used for *mapping evidence only*; assembly still
    uses directed k31/k32 sequence transitions;
  * support is counted once per original physical read pair;
  * paired mates are recruited together and recruitment is iterated from newly
    assembled terminal sequence;
  * ordinary k31 graph edges need >=2 physical fragments;
  * high-quality singleton edges are eligible only when their real read path has
    Stage25-like singleton_density<=0.30 and >=40 solid k31 nodes;
  * production singleton paths additionally require singleton edges from >=2
    distinct physical fragments, retaining the Stage23 safeguard;
  * a bridge is emitted only when the observed local sequence reaches >=40 bp of
    another immutable seed end exactly.  No N gap or unobserved insert sequence
    is ever synthesized.

The output variants intentionally separate solid-only, conservative two-fragment
singleton rescue, and one-fragment diagnostic rescue so accuracy/recall can be
measured rather than hidden in one threshold.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import stage23_guided_singleton_rescue as s23
import stage31_multik_seed_bridge as s31

BASE = {65: 0, 67: 1, 71: 2, 84: 3}
BITS = b"ACGT"
K = 31
MASK31 = (1 << (2 * K)) - 1

VARIANTS = {
    "solid2": {"allow_mercy": False, "min_mercy_fragments": 0},
    "oracle2": {"allow_mercy": True, "min_mercy_fragments": 2},
    "oracle1": {"allow_mercy": True, "min_mercy_fragments": 1},
}


@dataclass(frozen=True)
class Fragment:
    fid: int
    r1: bytes
    q1: bytes
    r2: bytes
    q2: bytes


@dataclass(frozen=True)
class LocalResult:
    endpoint: tuple[int, str]
    sequence: bytes
    added_bp: int
    bucket_fragments: int
    min_edge_support: int
    mercy_edges: int
    mercy_fragments: tuple[int, ...]
    branch_decisions: int
    stop_reason: str

    @property
    def evidence_rank(self) -> tuple[int, int, int, int, int, int]:
        return (
            int(self.mercy_edges == 0),
            self.min_edge_support,
            len(self.mercy_fragments),
            self.bucket_fragments,
            self.added_bp,
            -self.branch_decisions,
        )


@dataclass(frozen=True)
class BridgeEdge:
    proposal: s31.P
    local: LocalResult

    @property
    def rank(self) -> tuple[int, ...]:
        return self.local.evidence_rank + (self.proposal.lo + self.proposal.ro, -len(self.proposal.mid))


def rc(seq: bytes) -> bytes:
    return s31.rc(seq)


def rc_qual(qual: bytes) -> bytes:
    return qual[::-1]


def endpoint_id(sid: int, end: str) -> int:
    return sid * 2 + (1 if end == "R" else 0)


def endpoint_from_id(eid: int) -> tuple[int, str]:
    return eid // 2, "R" if eid & 1 else "L"


def oriented_seed(seedrows: list[tuple[str, bytes]], eid: int) -> bytes:
    sid, end = endpoint_from_id(eid)
    seq = seedrows[sid][1]
    return seq if end == "R" else rc(seq)


def directed_key(seq: bytes) -> int | None:
    if len(seq) != K:
        return None
    z = 0
    for ch in seq:
        value = BASE.get(ch)
        if value is None:
            return None
        z = (z << 2) | value
    return z


def directed_nodes(seq: bytes) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    z = valid = 0
    for i, ch in enumerate(seq):
        value = BASE.get(ch)
        if value is None:
            z = valid = 0
            continue
        z = ((z << 2) | value) & MASK31
        valid += 1
        if valid >= K:
            out.append((i - K + 1, z))
    return out


def canonical_21_keys(seq: bytes) -> list[tuple[int, int]]:
    try:
        text = seq.decode("ascii")
    except UnicodeDecodeError:
        return []
    return list(s23.rolling_keys(text, 21))


def terminal_signature_index(
    seedrows: list[tuple[str, bytes]],
    window: int,
) -> tuple[dict[int, int], dict[int, set[int]], dict[str, int]]:
    owners: dict[int, set[int]] = defaultdict(set)
    endpoint_keys: dict[int, set[int]] = {}
    for sid in range(len(seedrows)):
        for end in ("L", "R"):
            eid = endpoint_id(sid, end)
            seq = oriented_seed(seedrows, eid)[-window:]
            keys = {key for _pos, key in canonical_21_keys(seq)}
            endpoint_keys[eid] = keys
            for key in keys:
                owners[key].add(eid)
    unique = {key: next(iter(eids)) for key, eids in owners.items() if len(eids) == 1}
    unique_by_endpoint: dict[int, set[int]] = defaultdict(set)
    for key, eid in unique.items():
        unique_by_endpoint[eid].add(key)
    stats = {
        "physical_endpoints": len(seedrows) * 2,
        "terminal_k21_total": len(owners),
        "terminal_k21_unique": len(unique),
        "endpoints_with_unique_k21": sum(bool(unique_by_endpoint[eid]) for eid in range(len(seedrows) * 2)),
    }
    return unique, unique_by_endpoint, stats


def active_tail_index(
    local: dict[int, LocalResult],
    tail_bp: int,
    min_added: int,
) -> tuple[dict[int, int], dict[str, int]]:
    owners: dict[int, set[int]] = defaultdict(set)
    active_endpoints = 0
    for eid, result in local.items():
        if result.added_bp < min_added:
            continue
        active_endpoints += 1
        tail = result.sequence[-tail_bp:]
        for _pos, key in canonical_21_keys(tail):
            owners[key].add(eid)
    unique = {key: next(iter(eids)) for key, eids in owners.items() if len(eids) == 1}
    return unique, {
        "active_endpoints": active_endpoints,
        "tail_k21_total": len(owners),
        "tail_k21_unique": len(unique),
    }


def read_signature_hits(
    seq: bytes,
    key_owner: dict[int, int],
    min_hits: int,
    min_span: int,
) -> set[int]:
    by_ep: dict[int, dict[int, int]] = defaultdict(dict)
    for pos, key in canonical_21_keys(seq):
        eid = key_owner.get(key)
        if eid is not None:
            by_ep[eid].setdefault(key, pos)
    out: set[int] = set()
    for eid, hits in by_ep.items():
        if len(hits) < min_hits:
            continue
        positions = list(hits.values())
        if max(positions) - min(positions) >= min_span:
            out.add(eid)
    return out


def fastq_pairs(read1: Path, read2: Path) -> Iterator[Fragment]:
    for fid, (left, right) in enumerate(zip(s23.fastq_records(read1), s23.fastq_records(read2)), 1):
        s1, q1 = left
        s2, q2 = right
        yield Fragment(fid, s1.encode(), q1.encode(), s2.encode(), q2.encode())


def recruit_round(
    read1: Path,
    read2: Path,
    key_owner: dict[int, int],
    buckets: dict[int, set[int]],
    fragments: dict[int, Fragment],
    min_hits: int,
    min_span: int,
) -> dict[str, int]:
    touched_fragments = new_fragments = assignments = new_assignments = 0
    for frag in fastq_pairs(read1, read2):
        hits = read_signature_hits(frag.r1, key_owner, min_hits, min_span)
        hits.update(read_signature_hits(frag.r2, key_owner, min_hits, min_span))
        if not hits:
            continue
        touched_fragments += 1
        if frag.fid not in fragments:
            fragments[frag.fid] = frag
            new_fragments += 1
        for eid in hits:
            assignments += 1
            bucket = buckets.setdefault(eid, set())
            if frag.fid not in bucket:
                bucket.add(frag.fid)
                new_assignments += 1
    return {
        "signature_keys": len(key_owner),
        "touched_fragments": touched_fragments,
        "new_fragments": new_fragments,
        "assignments": assignments,
        "new_assignments": new_assignments,
        "endpoints_with_buckets": sum(bool(v) for v in buckets.values()),
        "total_bucket_assignments": sum(len(v) for v in buckets.values()),
    }


def oriented_reads(frag: Fragment) -> tuple[tuple[bytes, bytes], ...]:
    return (
        (frag.r1, frag.q1),
        (rc(frag.r1), rc_qual(frag.q1)),
        (frag.r2, frag.q2),
        (rc(frag.r2), rc_qual(frag.q2)),
    )


def mean_phred(qual: bytes) -> float:
    return sum(max(0, q - 33) for q in qual) / max(1, len(qual))


def build_local_support(
    fids: set[int],
    fragments: dict[int, Fragment],
    singleton_density: float,
    min_solid_nodes: int,
    singleton_quality: float,
) -> tuple[Counter[int], Counter[tuple[int, int]], dict[tuple[int, int], int], dict[str, int]]:
    node_support: Counter[int] = Counter()
    edge_support: Counter[tuple[int, int]] = Counter()

    # First pass: one vote per physical fragment for each directed node/edge.
    for fid in fids:
        frag = fragments[fid]
        frag_nodes: set[int] = set()
        frag_edges: set[tuple[int, int]] = set()
        for seq, _qual in oriented_reads(frag):
            nodes = directed_nodes(seq)
            frag_nodes.update(key for _pos, key in nodes)
            for (p0, k0), (p1, k1) in zip(nodes, nodes[1:]):
                if p1 != p0 + 1:
                    continue
                base = BASE.get(seq[p1 + K - 1])
                if base is not None:
                    frag_edges.add((k0, base))
        node_support.update(frag_nodes)
        edge_support.update(frag_edges)

    # Second pass: identify Stage25-like high-quality singleton paths and record
    # the one physical fragment that can contribute mercy to each edge.
    mercy_owner: dict[tuple[int, int], int] = {}
    conflicts: set[tuple[int, int]] = set()
    eligible_reads = 0
    for fid in fids:
        frag = fragments[fid]
        for seq, qual in oriented_reads(frag):
            nodes = directed_nodes(seq)
            if not nodes or mean_phred(qual) < singleton_quality:
                continue
            solid = sum(node_support[key] >= 2 for _pos, key in nodes)
            singleton = sum(node_support[key] == 1 for _pos, key in nodes)
            density = singleton / max(1, len(nodes))
            if solid < min_solid_nodes or density > singleton_density:
                continue
            eligible_reads += 1
            for (p0, k0), (p1, _k1) in zip(nodes, nodes[1:]):
                if p1 != p0 + 1:
                    continue
                base = BASE.get(seq[p1 + K - 1])
                if base is None:
                    continue
                edge = (k0, base)
                if edge_support.get(edge, 0) != 1:
                    continue
                old = mercy_owner.get(edge)
                if old is None:
                    mercy_owner[edge] = fid
                elif old != fid:
                    conflicts.add(edge)
    for edge in conflicts:
        mercy_owner.pop(edge, None)
    return node_support, edge_support, mercy_owner, {
        "fragments": len(fids),
        "nodes": len(node_support),
        "edges": len(edge_support),
        "solid_nodes": sum(v >= 2 for v in node_support.values()),
        "solid_edges": sum(v >= 2 for v in edge_support.values()),
        "eligible_singleton_reads": eligible_reads,
        "eligible_singleton_edges": len(mercy_owner),
    }


def extend_local(
    seedrows: list[tuple[str, bytes]],
    eid: int,
    fids: set[int],
    fragments: dict[int, Fragment],
    *,
    seed_overlap: int,
    max_extension: int,
    dominance: float,
    allow_mercy: bool,
    min_mercy_fragments: int,
    singleton_density: float,
    min_solid_nodes: int,
    singleton_quality: float,
) -> tuple[LocalResult, dict[str, int]]:
    source = oriented_seed(seedrows, eid)
    seed_tail = source[-seed_overlap:]
    start = directed_key(source[-K:])
    if start is None or not fids:
        result = LocalResult(endpoint_from_id(eid), seed_tail, 0, len(fids), 0, 0, (), 0, "no_bucket")
        return result, {"fragments": len(fids), "nodes": 0, "edges": 0, "solid_nodes": 0, "solid_edges": 0, "eligible_singleton_reads": 0, "eligible_singleton_edges": 0}

    _nodes, edges, mercy_owner, stats = build_local_support(
        fids,
        fragments,
        singleton_density,
        min_solid_nodes,
        singleton_quality,
    )
    node = start
    added = bytearray()
    visited = {node}
    min_support = 10**9
    mercy_fids: set[int] = set()
    mercy_edges = branch_decisions = 0
    stop = "max_extension"

    for _ in range(max_extension):
        solid: list[tuple[int, int]] = []
        mercy: list[tuple[int, int, int]] = []
        for base in range(4):
            edge = (node, base)
            support = edges.get(edge, 0)
            if support >= 2:
                solid.append((support, base))
            elif allow_mercy and support == 1 and edge in mercy_owner:
                mercy.append((support, base, mercy_owner[edge]))
        chosen_base: int | None = None
        chosen_support = 0
        chosen_mercy: int | None = None
        if solid:
            solid.sort(reverse=True)
            if len(solid) > 1:
                if solid[0][0] == solid[1][0] or solid[0][0] < math.ceil(dominance * solid[1][0]):
                    stop = "ambiguous_solid_branch"
                    break
                branch_decisions += 1
            chosen_support, chosen_base = solid[0]
        elif mercy:
            # Singleton rescue cannot resolve a branch: a single physical
            # fragment is permitted to repair a missing count, not choose among
            # alternative topologies.
            unique_bases = {base for _support, base, _fid in mercy}
            if len(unique_bases) != 1:
                stop = "ambiguous_singleton_branch"
                break
            _support, chosen_base, chosen_mercy = mercy[0]
            chosen_support = 1
            mercy_edges += 1
            mercy_fids.add(chosen_mercy)
        else:
            stop = "no_supported_edge"
            break
        assert chosen_base is not None
        next_node = ((node << 2) | chosen_base) & MASK31
        if next_node in visited:
            stop = "cycle"
            break
        visited.add(next_node)
        added.append(BITS[chosen_base])
        min_support = min(min_support, chosen_support)
        node = next_node

    if mercy_edges and len(mercy_fids) < min_mercy_fragments:
        # Conservative production gate: discard the entire mercy-derived tail
        # rather than pretending a single fragment can create a bridge.
        return LocalResult(endpoint_from_id(eid), seed_tail, 0, len(fids), 0, mercy_edges, tuple(sorted(mercy_fids)), branch_decisions, "insufficient_mercy_fragments"), stats

    sequence = seed_tail + bytes(added)
    return LocalResult(
        endpoint_from_id(eid),
        sequence,
        len(added),
        len(fids),
        0 if min_support == 10**9 else min_support,
        mercy_edges,
        tuple(sorted(mercy_fids)),
        branch_decisions,
        stop,
    ), stats


def assemble_all_endpoints(
    seedrows: list[tuple[str, bytes]],
    buckets: dict[int, set[int]],
    fragments: dict[int, Fragment],
    cfg: dict[str, object],
    args: argparse.Namespace,
) -> tuple[dict[int, LocalResult], dict[str, object]]:
    results: dict[int, LocalResult] = {}
    support_summary: Counter[str] = Counter()
    added = []
    for eid, fids in buckets.items():
        if not fids:
            continue
        result, stats = extend_local(
            seedrows,
            eid,
            fids,
            fragments,
            seed_overlap=args.seed_overlap,
            max_extension=args.max_extension,
            dominance=args.edge_dominance,
            allow_mercy=bool(cfg["allow_mercy"]),
            min_mercy_fragments=int(cfg["min_mercy_fragments"]),
            singleton_density=args.max_singleton_density,
            min_solid_nodes=args.min_solid_nodes,
            singleton_quality=args.singleton_quality,
        )
        results[eid] = result
        added.append(result.added_bp)
        for key in ("fragments", "nodes", "edges", "solid_nodes", "solid_edges", "eligible_singleton_reads", "eligible_singleton_edges"):
            support_summary[key] += int(stats[key])
    stats_out: dict[str, object] = dict(support_summary)
    stats_out.update(
        {
            "endpoints_assembled": len(results),
            "endpoints_extended": sum(r.added_bp > 0 for r in results.values()),
            "endpoints_extended_ge_50": sum(r.added_bp >= 50 for r in results.values()),
            "endpoints_extended_ge_100": sum(r.added_bp >= 100 for r in results.values()),
            "endpoints_extended_ge_250": sum(r.added_bp >= 250 for r in results.values()),
            "total_added_endpoint_bp": sum(added),
            "max_endpoint_extension": max(added, default=0),
            "mercy_endpoints": sum(r.mercy_edges > 0 for r in results.values()),
            "stop_reasons": dict(Counter(r.stop_reason for r in results.values())),
        }
    )
    return results, stats_out


def parse_source_eid(name: str) -> int | None:
    marker = "eid="
    pos = name.find(marker)
    if pos < 0:
        return None
    pos += len(marker)
    end = pos
    while end < len(name) and name[end].isdigit():
        end += 1
    return int(name[pos:end]) if end > pos else None


def discover_bridges(
    seedrows: list[tuple[str, bytes]],
    local: dict[int, LocalResult],
    min_overlap: int,
    overlap_margin: int,
) -> tuple[list[BridgeEdge], dict[str, int]]:
    candidates = [
        (f"stage33_eid={eid} added={r.added_bp} mercy={r.mercy_edges}", r.sequence)
        for eid, r in local.items()
        if r.added_bp > 0
    ]
    props, raw_stats = s31.discover(
        seedrows,
        candidates,
        K,
        minov=min_overlap,
        margin=overlap_margin,
    )
    accepted: list[BridgeEdge] = []
    wrong_source = 0
    for p in props:
        eid = parse_source_eid(p.name)
        if eid is None:
            wrong_source += 1
            continue
        result = local.get(eid)
        if result is None or p.le != result.endpoint:
            wrong_source += 1
            continue
        accepted.append(BridgeEdge(p, result))
    return accepted, {
        **raw_stats,
        "source_consistent_proposals": len(accepted),
        "wrong_source_rejected": wrong_source,
    }


def select_bridges(
    edges: list[BridgeEdge],
    nseed: int,
    max_component_seeds: int,
) -> tuple[list[BridgeEdge], dict[str, int]]:
    # First collapse exact duplicate sequence hypotheses for one physical pair.
    by_pair_mid: dict[tuple[tuple[tuple[int, str], tuple[int, str]], bytes], list[BridgeEdge]] = defaultdict(list)
    for edge in edges:
        by_pair_mid[(edge.proposal.pair, edge.proposal.midkey)].append(edge)
    unique_sequences: list[BridgeEdge] = []
    for vals in by_pair_mid.values():
        unique_sequences.append(max(vals, key=lambda e: e.rank))

    # Competing sequences for the same endpoint pair must have a unique evidence
    # winner. Equal-rank alternatives are rejected rather than guessed.
    by_pair: dict[tuple[tuple[int, str], tuple[int, str]], list[BridgeEdge]] = defaultdict(list)
    for edge in unique_sequences:
        by_pair[edge.proposal.pair].append(edge)
    pair_best: list[BridgeEdge] = []
    sequence_ties = 0
    for vals in by_pair.values():
        vals.sort(key=lambda e: (e.rank, e.proposal.midkey), reverse=True)
        if len(vals) > 1 and vals[0].rank == vals[1].rank:
            sequence_ties += 1
            continue
        pair_best.append(vals[0])

    incident: dict[tuple[int, str], list[BridgeEdge]] = defaultdict(list)
    for edge in pair_best:
        incident[edge.proposal.le].append(edge)
        incident[edge.proposal.re].append(edge)
    winner: dict[tuple[int, str], BridgeEdge | None] = {}
    endpoint_ties = 0
    for ep, vals in incident.items():
        vals.sort(key=lambda e: (e.rank, e.proposal.pair), reverse=True)
        if len(vals) > 1 and vals[0].rank == vals[1].rank:
            winner[ep] = None
            endpoint_ties += 1
        else:
            winner[ep] = vals[0]
    reciprocal = [
        edge for edge in pair_best
        if winner.get(edge.proposal.le) is edge and winner.get(edge.proposal.re) is edge
    ]

    parent = list(range(nseed))
    size = [1] * nseed

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    chosen: list[BridgeEdge] = []
    cycles = component_rejected = 0
    reciprocal.sort(key=lambda e: (e.rank, e.proposal.pair), reverse=True)
    for edge in reciprocal:
        a, b = edge.proposal.le[0], edge.proposal.re[0]
        ra, rb = find(a), find(b)
        if ra == rb:
            cycles += 1
            continue
        if size[ra] + size[rb] > max_component_seeds:
            component_rejected += 1
            continue
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]
        chosen.append(edge)
    return chosen, {
        "raw_edges": len(edges),
        "physical_endpoint_pairs": len(by_pair),
        "sequence_ties": sequence_ties,
        "pair_best": len(pair_best),
        "endpoint_ties": endpoint_ties,
        "reciprocal_edges": len(reciprocal),
        "cycle_rejected": cycles,
        "component_size_rejected": component_rejected,
        "selected_bridges": len(chosen),
    }


def write_local_fasta(path: Path, local: dict[int, LocalResult]) -> None:
    rows = []
    for eid, result in sorted(local.items()):
        sid, end = result.endpoint
        rows.append((
            f"stage33_eid={eid} seed={sid+1} end={end} added={result.added_bp} mercy={result.mercy_edges} fragments={result.bucket_fragments}",
            result.sequence,
        ))
    s31.write_fa(path, rows)


def run_seed_lock(scripts: Path, output: Path, seed: Path, strict: Path, stats: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(scripts / "seed_locked_extensions.py"),
            str(output),
            str(seed),
            str(strict),
            "--min-overlap", "500",
            "--overlap-margin", "30",
            "--seed-length", "31",
            "--min-extension", "20",
            "--max-seed-occurrences", "64",
            "--stats-json", str(stats),
        ],
        check=True,
    )


def max_chain_seeds(rows: list[tuple[str, bytes]]) -> int:
    return max(
        (name.count(",") + 1 for name, _seq in rows if name.startswith("stage31_chain")),
        default=1,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--seed", type=Path, required=True)
    ap.add_argument("--strict-contigs", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--terminal-window", type=int, default=240)
    ap.add_argument("--signature-hits", type=int, default=2)
    ap.add_argument("--signature-span", type=int, default=28)
    ap.add_argument("--active-tail", type=int, default=100)
    ap.add_argument("--active-min-added", type=int, default=20)
    ap.add_argument("--seed-overlap", type=int, default=120)
    ap.add_argument("--max-extension", type=int, default=1500)
    ap.add_argument("--edge-dominance", type=float, default=1.75)
    ap.add_argument("--max-singleton-density", type=float, default=0.30)
    ap.add_argument("--min-solid-nodes", type=int, default=40)
    ap.add_argument("--singleton-quality", type=float, default=30.0)
    ap.add_argument("--bridge-min-overlap", type=int, default=40)
    ap.add_argument("--bridge-overlap-margin", type=int, default=10)
    ap.add_argument("--max-component-seeds", type=int, default=4)
    args = ap.parse_args()

    if args.seed_overlap < K or args.max_extension < 1:
        raise SystemExit("invalid local assembly geometry")
    args.output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    seedrows = list(s31.fasta(args.seed))
    if not seedrows or any(set(seq) - set(b"ACGT") for _name, seq in seedrows):
        raise SystemExit("Stage33 requires a nonempty ACGT-only immutable seed")

    seed_key_owner, _seed_unique_by_ep, signature_stats = terminal_signature_index(
        seedrows, args.terminal_window
    )
    buckets: dict[int, set[int]] = {}
    fragments: dict[int, Fragment] = {}
    recruitment: list[dict[str, object]] = []

    first = recruit_round(
        args.read1, args.read2, seed_key_owner, buckets, fragments,
        args.signature_hits, args.signature_span,
    )
    recruitment.append({"round": 0, "source": "immutable_seed_terminal_k21", **first})

    # Use the one-fragment oracle mode only as a *recruitment probe*. Every final
    # variant is rebuilt independently from the resulting real fragment buckets
    # and enforces its own solid/singleton policy.
    probe_cfg = VARIANTS["oracle1"]
    probe_local: dict[int, LocalResult] = {}
    probe_stats: dict[str, object] = {}
    for round_index in range(1, args.rounds + 1):
        probe_local, probe_stats = assemble_all_endpoints(
            seedrows, buckets, fragments, probe_cfg, args
        )
        if round_index == args.rounds:
            break
        active, active_stats = active_tail_index(
            probe_local, args.active_tail, args.active_min_added
        )
        if not active:
            recruitment.append({"round": round_index, "source": "assembled_tail_k21", **active_stats, "new_assignments": 0})
            break
        rec = recruit_round(
            args.read1, args.read2, active, buckets, fragments,
            args.signature_hits, args.signature_span,
        )
        recruitment.append({"round": round_index, "source": "assembled_tail_k21", **active_stats, **rec})
        if rec["new_assignments"] == 0:
            break

    variants: dict[str, object] = {}
    for name, cfg in VARIANTS.items():
        local, local_stats = assemble_all_endpoints(seedrows, buckets, fragments, cfg, args)
        proposals, discover_stats = discover_bridges(
            seedrows, local, args.bridge_min_overlap, args.bridge_overlap_margin
        )
        chosen, select_stats = select_bridges(
            proposals, len(seedrows), args.max_component_seeds
        )
        outdir = args.output / name
        outdir.mkdir(parents=True, exist_ok=True)
        write_local_fasta(outdir / "endpoint_extensions.fasta", local)
        selected31 = [(edge.proposal, (31,), 1) for edge in chosen]
        bridged = s31.assemble(seedrows, selected31)
        bridged_path = outdir / "bridged_seed.fasta"
        s31.write_fa(bridged_path, bridged)
        rows = [
            "left_seed\tleft_end\tright_seed\tright_end\tinternal_bp\tadded_bp\tmin_edge_support\tmercy_edges\tmercy_fragments\tbucket_fragments\trank"
        ]
        for edge in chosen:
            p, r = edge.proposal, edge.local
            rows.append(
                f"{p.le[0]+1}\t{p.le[1]}\t{p.re[0]+1}\t{p.re[1]}\t{len(p.mid)}\t{r.added_bp}\t"
                f"{r.min_edge_support}\t{r.mercy_edges}\t{len(r.mercy_fragments)}\t{r.bucket_fragments}\t"
                f"{','.join(map(str,edge.rank))}"
            )
        (outdir / "bridges.tsv").write_text("\n".join(rows) + "\n")
        extension_stats = outdir / "extension_stats.json"
        final = outdir / "primary_contigs.fasta"
        run_seed_lock(scripts, final, bridged_path, args.strict_contigs, extension_stats)
        bridge_stats = {
            **discover_stats,
            **select_stats,
            "seed_records": len(seedrows),
            "bridged_records": len(bridged),
            "seed_n50": s31.n50(seq for _name, seq in seedrows),
            "bridge_n50": s31.n50(seq for _name, seq in bridged),
            "max_chain_seeds": max_chain_seeds(bridged),
        }
        (outdir / "local_stats.json").write_text(json.dumps(local_stats, indent=2, sort_keys=True) + "\n")
        (outdir / "bridge_stats.json").write_text(json.dumps(bridge_stats, indent=2, sort_keys=True) + "\n")
        variants[name] = {
            "local": local_stats,
            "bridge": bridge_stats,
            "extension": json.loads(extension_stats.read_text()),
            "final": str(final),
        }

    bucket_sizes = sorted((len(v) for v in buckets.values()), reverse=True)
    manifest = {
        "pipeline": "stage33-fragment-aware-seed-end-microassembly-v1",
        "signature": signature_stats,
        "recruitment": recruitment,
        "probe": probe_stats,
        "bucket_summary": {
            "fragments_recruited": len(fragments),
            "endpoints_with_buckets": len(buckets),
            "total_assignments": sum(bucket_sizes),
            "median_bucket": bucket_sizes[len(bucket_sizes)//2] if bucket_sizes else 0,
            "p90_bucket": bucket_sizes[max(0, len(bucket_sizes)//10 - 1)] if bucket_sizes else 0,
            "max_bucket": max(bucket_sizes, default=0),
        },
        "parameters": {
            "rounds": args.rounds,
            "terminal_window": args.terminal_window,
            "signature_hits": args.signature_hits,
            "signature_span": args.signature_span,
            "active_tail": args.active_tail,
            "seed_overlap": args.seed_overlap,
            "max_extension": args.max_extension,
            "edge_dominance": args.edge_dominance,
            "max_singleton_density": args.max_singleton_density,
            "min_solid_nodes": args.min_solid_nodes,
            "singleton_quality": args.singleton_quality,
            "bridge_min_overlap": args.bridge_min_overlap,
            "max_component_seeds": args.max_component_seeds,
        },
        "variants": variants,
    }
    (args.output / "stage33_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
