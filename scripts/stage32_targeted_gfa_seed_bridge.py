#!/usr/bin/env python3
"""Stage32: targeted evidence-guided seed-end bridging on the k31 GFA.

Stage31 showed that exact whole-contig seed bridges are safe but extremely rare.
This stage changes only bridge discovery: immutable Stage24 seed ends are mapped
back to unique exact anchors in the k31 unitig graph, then each end is walked
outward locally. Unique topology is followed directly; ambiguous exits must be
resolved by the proven Stage8 raw-read/high-k/pair-context bounded lookahead.
Only paths that reach another exact seed end are materialized, so no graph edge
or N-gap is invented. Selection is reciprocal at physical seed ends, cycle-free,
and component size is bounded before the conservative k31 seed-lock extension.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import graph_path_phaser as gp
import repeat_graph_optimizer as rg
import stage31_multik_seed_bridge as s31
import stage789_optimizer as s78

PASS_ORDER = ("rawhigh", "paired", "simplified")
VARIANTS = {
    "consensus3": {"min_passes": 3, "require_pair": True},
    "consensus2": {"min_passes": 2, "require_pair": True},
    "paired1": {"min_passes": 1, "require_pair": True},
    "any_graph": {"min_passes": 1, "require_pair": False},
}


@dataclass(frozen=True)
class Anchor:
    state: int
    source: bool
    node: str
    seed_pos: int
    unitig_pos: int
    overlap: int

    @property
    def sid(self) -> int:
        return self.state // 2

    @property
    def endpoint(self) -> tuple[int, str]:
        return s31.endpoint(self.state, self.source)


@dataclass(frozen=True)
class GraphProposal:
    p: s31.P
    evidence_pass: str
    path_nodes: tuple[str, ...]
    branch_decisions: int
    strong_branches: int
    min_edge_support: int
    sum_edge_support: int


@dataclass
class CandidateEdge:
    p: s31.P
    passes: tuple[str, ...]
    evidence_rank: tuple[int, int, int, int, int, int]
    proposals: list[GraphProposal]


def oriented_seed(seedrows: list[tuple[str, bytes]], state: int) -> bytes:
    seq = seedrows[state // 2][1]
    return s31.rc(seq) if state & 1 else seq


def key31(seq: bytes) -> int | None:
    if len(seq) != 31:
        return None
    try:
        text = seq.decode("ascii")
    except UnicodeDecodeError:
        return None
    vals = list(gp.rolling_keys(text, 31))
    return vals[0][1] if len(vals) == 1 else None


def anchor_candidates(
    seedrows: list[tuple[str, bytes]],
    graph: gp.Graph,
    index: gp.KmerIndex,
    min_overlap: int,
    max_overlap: int,
) -> dict[tuple[tuple[int, str], bool], list[Anchor]]:
    """Return deepest exact unique graph anchors for every physical seed end.

    A physical end has one source orientation and one target orientation.  We
    retain multiple candidate anchors until cross-seed anchor collisions are
    known, then choose the deepest collision-free candidate.
    """
    out: dict[tuple[tuple[int, str], bool], list[Anchor]] = defaultdict(list)
    for sid in range(len(seedrows)):
        for rev in (0, 1):
            state = sid * 2 + rev
            seq = oriented_seed(seedrows, state)
            upper = min(max_overlap, len(seq))
            if upper < min_overlap:
                continue
            for source in (True, False):
                ep = s31.endpoint(state, source)
                seen: set[tuple[str, int]] = set()
                for overlap in range(upper, min_overlap - 1, -1):
                    seed_pos = len(seq) - overlap if source else overlap - 31
                    if seed_pos < 0 or seed_pos + 31 > len(seq):
                        continue
                    kmer = seq[seed_pos : seed_pos + 31]
                    key = key31(kmer)
                    if key is None:
                        continue
                    node = index.unique.get(key)
                    if node is None:
                        continue
                    useq = graph.seqs[node].encode()
                    unitig_pos = useq.find(kmer)
                    if unitig_pos < 0 or useq.rfind(kmer) != unitig_pos:
                        continue
                    loc = (node, unitig_pos)
                    if loc in seen:
                        continue
                    seen.add(loc)
                    # Verify as much same-seed sequence as is present inside the
                    # anchor unitig. This rejects unique 31-mers that happen to
                    # occur in a wrong local context.
                    if source:
                        take = min(len(useq) - unitig_pos, len(seq) - seed_pos)
                        if take < 31 or useq[unitig_pos : unitig_pos + take] != seq[seed_pos : seed_pos + take]:
                            continue
                    else:
                        left = min(unitig_pos, seed_pos)
                        if useq[unitig_pos - left : unitig_pos + 31] != seq[seed_pos - left : seed_pos + 31]:
                            continue
                    out[(ep, source)].append(
                        Anchor(state, source, node, seed_pos, unitig_pos, overlap)
                    )
    return out


def choose_anchors(
    candidates: dict[tuple[tuple[int, str], bool], list[Anchor]],
) -> tuple[list[Anchor], dict[str, int]]:
    # An anchor location is usable only when it identifies one physical endpoint
    # for this direction. A repetitive terminal k-mer can therefore fall back to
    # a slightly deeper/shallower unique candidate instead of deleting the seed.
    owners: dict[tuple[bool, str, int], set[tuple[int, str]]] = defaultdict(set)
    for (ep, source), vals in candidates.items():
        for a in vals:
            owners[(source, a.node, a.unitig_pos)].add(ep)
    chosen: list[Anchor] = []
    missing = collision_rejected = 0
    for (ep, source), vals in sorted(candidates.items()):
        pick = None
        for a in vals:
            if len(owners[(source, a.node, a.unitig_pos)]) == 1:
                pick = a
                break
            collision_rejected += 1
        if pick is None:
            missing += 1
        else:
            chosen.append(pick)
    expected = len({key for key in candidates})
    return chosen, {
        "candidate_endpoint_directions": expected,
        "chosen_anchors": len(chosen),
        "missing_after_collision_filter": missing,
        "collision_candidates_rejected": collision_rejected,
    }


def edge_physical(graph: gp.Graph, u: str, v: str) -> int:
    ev = graph.edge.get((u, v), gp.EdgeEvidence())
    return max(ev.direct, ev.gapped, ev.pairs)


def path_core(
    graph: gp.Graph, path: list[str], source: Anchor, target: Anchor
) -> bytes | None:
    if not path or path[0] != source.node or path[-1] != target.node:
        return None
    if len(path) == 1:
        if target.unitig_pos + 31 <= source.unitig_pos:
            return None
        return graph.seqs[path[0]][source.unitig_pos : target.unitig_pos + 31].encode()
    first = graph.seqs[path[0]]
    last = graph.seqs[path[-1]]
    target_end = target.unitig_pos + 31
    # Bases before graph.k in a downstream unitig are already represented by
    # the preceding k-overlap. If the target anchor ends entirely inside that
    # overlap, the path has already passed the physical endpoint ambiguously.
    if target_end <= graph.k:
        return None
    chunks = [first[source.unitig_pos :]]
    for uid in path[1:-1]:
        chunks.append(graph.seqs[uid][graph.k :])
    chunks.append(last[graph.k : target_end])
    return "".join(chunks).encode()


def proposal_from_path(
    seedrows: list[tuple[str, bytes]],
    graph: gp.Graph,
    path: list[str],
    source: Anchor,
    target: Anchor,
    evidence_pass: str,
    branch_decisions: int,
    strong_branches: int,
) -> GraphProposal | None:
    if source.sid == target.sid or source.endpoint == target.endpoint:
        return None
    core = path_core(graph, path, source, target)
    if core is None:
        return None
    left = oriented_seed(seedrows, source.state)
    right = oriented_seed(seedrows, target.state)
    left_suffix = left[source.seed_pos :]
    right_prefix = right[: target.seed_pos + 31]
    lo, ro = len(left_suffix), len(right_prefix)
    if len(core) <= lo + ro:
        return None
    if not core.startswith(left_suffix) or not core.endswith(right_prefix):
        return None
    mid = core[lo : len(core) - ro]
    if not mid or set(core) - set(b"ACGT"):
        return None
    supports = [edge_physical(graph, u, v) for u, v in zip(path, path[1:])]
    p = s31.P(
        31,
        source.state,
        target.state,
        lo,
        ro,
        mid,
        core,
        f"stage32:{evidence_pass}:{source.node}->{target.node}",
    )
    return GraphProposal(
        p=p,
        evidence_pass=evidence_pass,
        path_nodes=tuple(path),
        branch_decisions=branch_decisions,
        strong_branches=strong_branches,
        min_edge_support=min(supports) if supports else 0,
        sum_edge_support=sum(supports),
    )


def discover_pass(
    seedrows: list[tuple[str, bytes]],
    graph: gp.Graph,
    anchors: list[Anchor],
    raw_ctx: Counter[tuple[str, ...]],
    high_ctx: Counter[tuple[str, ...]],
    repeat_ctx: Counter[tuple[str, ...]],
    evidence_pass: str,
    max_nodes: int,
    max_bridge_bp: int,
    dominance: float = 0.70,
    min_direct: int = 4,
    lookahead_depth: int = 3,
    lookahead_max_branch: int = 4,
    lookahead_discount: float = 0.70,
    lookahead_dominance: float = 0.60,
    lookahead_margin: float = 1.15,
) -> tuple[list[GraphProposal], dict[str, int]]:
    targets: dict[str, list[Anchor]] = defaultdict(list)
    sources = []
    for a in anchors:
        (sources if a.source else targets[a.node]).append(a)
    proposals: list[GraphProposal] = []
    stats = {
        "sources": len(sources),
        "reached_seed_end": 0,
        "ambiguous_target_stop": 0,
        "branch_decisions": 0,
        "branch_stops": 0,
        "span_stops": 0,
        "node_stops": 0,
        "proposals": 0,
    }
    empty: Counter[tuple[str, ...]] = Counter()
    seen: set[tuple[tuple[tuple[int, str], tuple[int, str]], bytes]] = set()

    for source in sources:
        if source.node not in graph.seqs:
            continue
        path = [source.node]
        local_seen = {source.node, graph.rev.get(source.node, source.node)}
        branch_decisions = strong_branches = 0
        walked_bp = max(0, len(graph.seqs[source.node]) - source.unitig_pos)
        for _step in range(max_nodes):
            current = path[-1]
            hits: list[tuple[Anchor, GraphProposal]] = []
            for target in targets.get(current, []):
                q = proposal_from_path(
                    seedrows,
                    graph,
                    path,
                    source,
                    target,
                    evidence_pass,
                    branch_decisions,
                    strong_branches,
                )
                if q is not None:
                    hits.append((target, q))
            if hits:
                by_ep = {target.endpoint for target, _ in hits}
                if len(by_ep) != 1:
                    stats["ambiguous_target_stop"] += 1
                    break
                q = sorted(
                    (q for _, q in hits),
                    key=lambda x: (-x.p.lo - x.p.ro, len(x.p.mid), x.p.pair),
                )[0]
                key = (q.p.pair, q.p.midkey)
                if key not in seen:
                    seen.add(key)
                    proposals.append(q)
                stats["reached_seed_end"] += 1
                break

            children = [
                uid
                for uid in graph.out.get(current, [])
                if uid not in local_seen and graph.rev.get(uid, uid) not in local_seen
            ]
            if not children:
                break
            if len(children) == 1:
                nxt = children[0]
            else:
                choice, _rescued = s78.choose_extension_lookahead(
                    graph,
                    path,
                    children,
                    local_seen,
                    True,
                    raw_ctx,
                    empty,
                    high_ctx,
                    repeat_ctx,
                    dominance,
                    min_direct,
                    lookahead_depth,
                    lookahead_max_branch,
                    lookahead_discount,
                    lookahead_dominance,
                    lookahead_margin,
                )
                if choice is None:
                    stats["branch_stops"] += 1
                    break
                nxt = choice.uid
                branch_decisions += 1
                strong_branches += int(gp.strong_context(choice))
                stats["branch_decisions"] += 1
            walked_bp += max(1, len(graph.seqs[nxt]) - graph.k)
            if walked_bp > max_bridge_bp:
                stats["span_stops"] += 1
                break
            path.append(nxt)
            local_seen.add(nxt)
            local_seen.add(graph.rev.get(nxt, nxt))
        else:
            stats["node_stops"] += 1
    stats["proposals"] = len(proposals)
    return proposals, stats


def aggregate_edges(
    proposals: Iterable[GraphProposal],
    min_passes: int,
    require_pair: bool,
) -> tuple[list[CandidateEdge], dict[str, int]]:
    groups: dict[
        tuple[tuple[tuple[int, str], tuple[int, str]], bytes], list[GraphProposal]
    ] = defaultdict(list)
    for q in proposals:
        groups[(q.p.pair, q.p.midkey)].append(q)
    edges: list[CandidateEdge] = []
    rejected_passes = rejected_pair = 0
    for _key, vals in groups.items():
        passes = tuple(p for p in PASS_ORDER if any(v.evidence_pass == p for v in vals))
        if len(passes) < min_passes:
            rejected_passes += 1
            continue
        pair_supported = "paired" in passes or "simplified" in passes
        if require_pair and not pair_supported:
            rejected_pair += 1
            continue
        rep = max(
            vals,
            key=lambda v: (
                v.strong_branches,
                v.min_edge_support,
                v.sum_edge_support,
                -len(v.p.mid),
            ),
        )
        rank = (
            len(passes),
            int("simplified" in passes),
            int("paired" in passes),
            max(v.strong_branches for v in vals),
            max(v.min_edge_support for v in vals),
            max(v.sum_edge_support for v in vals),
        )
        edges.append(CandidateEdge(rep.p, passes, rank, vals))
    return edges, {
        "exact_bridge_hypotheses": len(groups),
        "candidate_edges": len(edges),
        "rejected_min_passes": rejected_passes,
        "rejected_pair_requirement": rejected_pair,
    }


def select_edges(
    edges: list[CandidateEdge], nseed: int, max_component_seeds: int
) -> tuple[list[CandidateEdge], dict[str, int]]:
    by_pair: dict[tuple[tuple[int, str], tuple[int, str]], list[CandidateEdge]] = defaultdict(list)
    for e in edges:
        by_pair[e.p.pair].append(e)
    sequence_unresolved = 0
    pair_best: list[CandidateEdge] = []
    for _pair, vals in by_pair.items():
        vals.sort(key=lambda e: (e.evidence_rank, -len(e.p.mid)), reverse=True)
        if len(vals) > 1 and vals[0].evidence_rank == vals[1].evidence_rank:
            sequence_unresolved += 1
            continue
        pair_best.append(vals[0])

    incident: dict[tuple[int, str], list[CandidateEdge]] = defaultdict(list)
    for e in pair_best:
        incident[e.p.le].append(e)
        incident[e.p.re].append(e)
    winners: dict[tuple[int, str], CandidateEdge | None] = {}
    endpoint_ties = 0
    for ep, vals in incident.items():
        vals.sort(key=lambda e: (e.evidence_rank, -len(e.p.mid)), reverse=True)
        if len(vals) > 1 and vals[0].evidence_rank == vals[1].evidence_rank:
            winners[ep] = None
            endpoint_ties += 1
        else:
            winners[ep] = vals[0]
    reciprocal = [
        e
        for e in pair_best
        if winners.get(e.p.le) is e and winners.get(e.p.re) is e
    ]

    parent = list(range(nseed))
    size = [1] * nseed

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    chosen: list[CandidateEdge] = []
    cycle_rejected = component_rejected = 0
    reciprocal.sort(key=lambda e: (e.evidence_rank, -len(e.p.mid)), reverse=True)
    for e in reciprocal:
        a, b = e.p.le[0], e.p.re[0]
        ra, rb = find(a), find(b)
        if ra == rb:
            cycle_rejected += 1
            continue
        if size[ra] + size[rb] > max_component_seeds:
            component_rejected += 1
            continue
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]
        chosen.append(e)
    return chosen, {
        "physical_endpoint_pairs": len(by_pair),
        "sequence_unresolved": sequence_unresolved,
        "pair_best": len(pair_best),
        "endpoint_ties": endpoint_ties,
        "reciprocal_edges": len(reciprocal),
        "cycle_rejected": cycle_rejected,
        "component_size_rejected": component_rejected,
        "selected_bridges": len(chosen),
    }


def subtract_counter(base: Counter, second: Counter) -> Counter:
    out = Counter(second)
    for key in list(out):
        if out[key] <= base.get(key, 0):
            del out[key]
        else:
            out[key] -= base.get(key, 0)
    return out


def build_contexts(
    graph: gp.Graph,
    read1: Path,
    read2: Path,
    highk_gfas: list[Path],
    work: Path,
) -> tuple[
    Counter[tuple[str, ...]],
    Counter[tuple[str, ...]],
    Counter[tuple[str, ...]],
    gp.Graph,
    dict[str, object],
]:
    index = gp.KmerIndex(graph, 31)
    raw_ctx, raw_stats = gp.collect_read_contexts(graph, index, read1, read2, None, 8)
    _proj, high_ctx, high_stats = rg.collect_projection_contexts(
        graph, index, [], highk_gfas, work / "highk_projection", 8
    )
    prelim, prelim_stats = s78.resolve_lookahead_seeded_paths(
        graph,
        raw_ctx,
        Counter(),
        high_ctx,
        Counter(),
        0.70,
        4,
        200,
        3,
        4,
        0.70,
        0.60,
        1.15,
    )
    membership = gp.preliminary_membership(prelim)
    second_all, second_stats = gp.collect_read_contexts(
        graph, index, read1, read2, membership, 8
    )
    second_ctx = subtract_counter(raw_ctx, second_all)
    pair_ctx, pair_stats = rg.collect_pair_contexts(
        graph, index, read1, read2, membership, 8, 6, 320
    )
    repeat_ctx = rg.combined_contexts(second_ctx, pair_ctx)
    simplified, simplify_stats = rg.simplify_graph(
        graph, rg.combined_contexts(raw_ctx, high_ctx, repeat_ctx)
    )
    stats: dict[str, object] = {
        "raw": raw_stats,
        "highk": high_stats,
        "preliminary_paths": prelim_stats,
        "second_pass": second_stats,
        "second_pass_novel_contexts": len(second_ctx),
        "pairs": pair_stats,
        "simplification": simplify_stats,
    }
    return raw_ctx, high_ctx, repeat_ctx, simplified, stats


def max_chain_seeds(rows: list[tuple[str, bytes]]) -> int:
    return max(
        (name.count(",") + 1 for name, _ in rows if name.startswith("stage31_chain")),
        default=1,
    )


def run_seed_lock(scripts: Path, output: Path, seed: Path, candidates: Path, stats: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(scripts / "seed_locked_extensions.py"),
            str(output),
            str(seed),
            str(candidates),
            "--min-overlap",
            "500",
            "--overlap-margin",
            "30",
            "--seed-length",
            "31",
            "--min-extension",
            "20",
            "--max-seed-occurrences",
            "64",
            "--stats-json",
            str(stats),
        ],
        check=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--read1", type=Path, required=True)
    ap.add_argument("--read2", type=Path, required=True)
    ap.add_argument("--seed", type=Path, required=True)
    ap.add_argument("--k31-gfa", type=Path, required=True)
    ap.add_argument("--strict-contigs", type=Path, required=True)
    ap.add_argument("--highk-gfa", type=Path, action="append", default=[])
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--anchor-min-overlap", type=int, default=80)
    ap.add_argument("--anchor-max-overlap", type=int, default=240)
    ap.add_argument("--max-bridge-bp", type=int, default=8000)
    ap.add_argument("--max-bridge-nodes", type=int, default=40)
    ap.add_argument("--max-component-seeds", type=int, default=4)
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent
    seedrows = list(s31.fasta(args.seed))
    if not seedrows or any(set(seq) - set(b"ACGT") for _, seq in seedrows):
        raise SystemExit("Stage32 requires a nonempty ACGT-only immutable seed")
    graph = gp.Graph.from_gfa(args.k31_gfa)
    if graph.k != 31:
        raise SystemExit(f"expected k31 GFA, observed overlap k={graph.k}")
    index = gp.KmerIndex(graph, 31)
    raw_candidates = anchor_candidates(
        seedrows, graph, index, args.anchor_min_overlap, args.anchor_max_overlap
    )
    anchors, anchor_stats = choose_anchors(raw_candidates)

    raw_ctx, high_ctx, repeat_ctx, simplified, context_stats = build_contexts(
        graph, args.read1, args.read2, args.highk_gfa, args.output
    )
    pass_cfg = {
        "rawhigh": (graph, Counter()),
        "paired": (graph, repeat_ctx),
        "simplified": (simplified, repeat_ctx),
    }
    proposals: list[GraphProposal] = []
    discovery_stats = {}
    for name in PASS_ORDER:
        use_graph, use_repeat = pass_cfg[name]
        vals, stats = discover_pass(
            seedrows,
            use_graph,
            anchors,
            raw_ctx,
            high_ctx,
            use_repeat,
            name,
            args.max_bridge_nodes,
            args.max_bridge_bp,
        )
        proposals.extend(vals)
        discovery_stats[name] = stats

    variants = {}
    for name, cfg in VARIANTS.items():
        candidates, agg_stats = aggregate_edges(
            proposals, cfg["min_passes"], cfg["require_pair"]
        )
        chosen, select_stats = select_edges(
            candidates, len(seedrows), args.max_component_seeds
        )
        bridge_dir = args.output / name
        bridge_dir.mkdir(parents=True, exist_ok=True)
        chosen31 = [
            (e.p, tuple(31 for _ in e.passes), e.evidence_rank[0]) for e in chosen
        ]
        bridged = s31.assemble(seedrows, chosen31)
        bridged_path = bridge_dir / "bridged_seed.fasta"
        s31.write_fa(bridged_path, bridged)
        rows = [
            "left_seed\tleft_end\tright_seed\tright_end\tpasses\tevidence_rank\tinternal_bp\tpath_nodes"
        ]
        for e in chosen:
            rep = max(
                e.proposals,
                key=lambda q: (
                    q.strong_branches,
                    q.min_edge_support,
                    q.sum_edge_support,
                    -len(q.p.mid),
                ),
            )
            rows.append(
                f"{e.p.le[0]+1}\t{e.p.le[1]}\t{e.p.re[0]+1}\t{e.p.re[1]}\t"
                f"{','.join(e.passes)}\t{','.join(map(str,e.evidence_rank))}\t"
                f"{len(e.p.mid)}\t{','.join(rep.path_nodes)}"
            )
        (bridge_dir / "bridges.tsv").write_text("\n".join(rows) + "\n")
        final = bridge_dir / "primary_contigs.fasta"
        extension_stats = bridge_dir / "extension_stats.json"
        run_seed_lock(scripts, final, bridged_path, args.strict_contigs, extension_stats)
        bridge_stats = {
            **agg_stats,
            **select_stats,
            "seed_records": len(seedrows),
            "bridged_records": len(bridged),
            "seed_n50": s31.n50(seq for _, seq in seedrows),
            "bridge_n50": s31.n50(seq for _, seq in bridged),
            "max_chain_seeds": max_chain_seeds(bridged),
        }
        (bridge_dir / "bridge_stats.json").write_text(
            json.dumps(bridge_stats, indent=2, sort_keys=True) + "\n"
        )
        variants[name] = {
            "bridge": bridge_stats,
            "extension": json.loads(extension_stats.read_text()),
            "final": str(final),
        }

    manifest = {
        "pipeline": "stage32-targeted-gfa-seed-bridge-v1",
        "anchor": anchor_stats,
        "contexts": context_stats,
        "discovery": discovery_stats,
        "total_graph_proposals": len(proposals),
        "parameters": {
            "anchor_min_overlap": args.anchor_min_overlap,
            "anchor_max_overlap": args.anchor_max_overlap,
            "max_bridge_bp": args.max_bridge_bp,
            "max_bridge_nodes": args.max_bridge_nodes,
            "max_component_seeds": args.max_component_seeds,
        },
        "variants": variants,
    }
    (args.output / "stage32_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
