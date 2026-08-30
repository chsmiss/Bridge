from pathlib import Path
import re

path = Path("src/assembler.rs")
text = path.read_text()

old = '''#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub read_pairs: u32,
}
'''
new = '''#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub read_pairs: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct TransitionCandidate {
    node: u32,
    direct_reads: u32,
    read_pairs: u32,
    solid_topology: bool,
}

#[derive(Clone, Copy, Debug)]
struct PathSelectionConfig {
    min_read_support: u32,
    min_pair_support: u32,
    min_primary_support: u32,
    min_count: u32,
    dominance: f32,
}
'''
if old not in text:
    raise SystemExit("TransitionEvidence block not found")
text = text.replace(old, new, 1)

old = '''        config.min_read_support,
        config.min_primary_support,
        config.primary_dominance,
    );
'''
new = '''        PathSelectionConfig {
            min_read_support: config.min_read_support,
            min_pair_support: config.min_pair_support,
            min_primary_support: config.min_primary_support,
            min_count: config.min_count,
            dominance: config.primary_dominance,
        },
    );
'''
if old not in text:
    raise SystemExit("primary_paths call not found")
text = text.replace(old, new, 1)

old = '''    min_read_support: u32,
    min_primary_support: u32,
    dominance: f32,
) -> (Vec<Vec<u32>>, usize) {
    let unitig_count = unitigs.unitigs.len();
    let excluded = non_primary_bubble_alleles(bubble_alleles);
    let mut outgoing_candidates: Vec<Vec<(u32, u32)>> = vec![Vec::new(); unitig_count];
    let mut incoming_candidates: Vec<Vec<(u32, u32)>> = vec![Vec::new(); unitig_count];
'''
new = '''    selection: PathSelectionConfig,
) -> (Vec<Vec<u32>>, usize) {
    let unitig_count = unitigs.unitigs.len();
    let excluded = non_primary_bubble_alleles(bubble_alleles);
    let mut outgoing_candidates: Vec<Vec<TransitionCandidate>> = vec![Vec::new(); unitig_count];
    let mut incoming_candidates: Vec<Vec<TransitionCandidate>> = vec![Vec::new(); unitig_count];
'''
if old not in text:
    raise SystemExit("primary_paths signature block not found")
text = text.replace(old, new, 1)

old = '''            let support = transitions
                .get(&(source, target))
                .map_or(0, |evidence| evidence.direct_reads);
            outgoing_candidates[source as usize].push((target, support));
            incoming_candidates[target as usize].push((source, support));
'''
new = '''            let evidence = transitions.get(&(source, target)).copied().unwrap_or_default();
            let solid_topology = unitigs.unitigs[source as usize].min_coverage >= selection.min_count
                && unitigs.unitigs[target as usize].min_coverage >= selection.min_count;
            outgoing_candidates[source as usize].push(TransitionCandidate {
                node: target,
                direct_reads: evidence.direct_reads,
                read_pairs: evidence.read_pairs,
                solid_topology,
            });
            incoming_candidates[target as usize].push(TransitionCandidate {
                node: source,
                direct_reads: evidence.direct_reads,
                read_pairs: evidence.read_pairs,
                solid_topology,
            });
'''
if old not in text:
    raise SystemExit("candidate construction block not found")
text = text.replace(old, new, 1)

old = '''            choose_transition(candidates, min_read_support, min_primary_support, dominance)
'''
new = '''            choose_transition(
                candidates,
                selection.min_read_support,
                selection.min_pair_support,
                selection.min_primary_support,
                selection.dominance,
            )
'''
if text.count(old) != 2:
    raise SystemExit(f"expected two choose_transition calls, found {text.count(old)}")
text = text.replace(old, new)

pattern = re.compile(
    r"fn choose_transition\(\n.*?\n\}\n\nfn extend_selected_path",
    re.S,
)
replacement = '''fn choose_transition(
    candidates: &[TransitionCandidate],
    min_read_support: u32,
    min_pair_support: u32,
    min_primary_support: u32,
    dominance: f32,
) -> Option<u32> {
    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        let candidate = candidates[0];
        // A unique retained edge has already passed k-mer/edge filtering. It
        // may be extended without a separately materialized read transition
        // when both incident unitigs are solid. Mercy-only links still need
        // direct-read or paired-fragment support.
        return (candidate.direct_reads >= min_read_support
            || candidate.read_pairs >= min_pair_support
            || candidate.solid_topology)
            .then_some(candidate.node);
    }

    let max_direct = candidates
        .iter()
        .map(|candidate| candidate.direct_reads)
        .max()
        .unwrap_or(0);
    let use_direct = max_direct >= min_read_support;
    let mut ranked = candidates.to_vec();
    if use_direct {
        ranked.sort_unstable_by(|left, right| {
            right
                .direct_reads
                .cmp(&left.direct_reads)
                .then_with(|| right.read_pairs.cmp(&left.read_pairs))
                .then_with(|| right.solid_topology.cmp(&left.solid_topology))
                .then_with(|| left.node.cmp(&right.node))
        });
    } else {
        ranked.sort_unstable_by(|left, right| {
            right
                .read_pairs
                .cmp(&left.read_pairs)
                .then_with(|| right.direct_reads.cmp(&left.direct_reads))
                .then_with(|| right.solid_topology.cmp(&left.solid_topology))
                .then_with(|| left.node.cmp(&right.node))
        });
    }

    let total: u64 = ranked
        .iter()
        .map(|candidate| {
            u64::from(if use_direct {
                candidate.direct_reads
            } else {
                candidate.read_pairs
            })
        })
        .sum();
    let best = ranked[0];
    let best_support = if use_direct {
        best.direct_reads
    } else {
        best.read_pairs
    };
    let minimum = if use_direct {
        min_primary_support
    } else {
        min_pair_support
    };
    let fraction = best_support as f32 / total.max(1) as f32;
    (best_support >= minimum && fraction >= dominance).then_some(best.node)
}

fn extend_selected_path'''
text, count = pattern.subn(replacement, text, count=1)
if count != 1:
    raise SystemExit(f"choose_transition replacements: {count}")

text += '''

#[cfg(test)]
mod transition_selection_tests {
    use super::*;

    fn candidate(
        node: u32,
        direct_reads: u32,
        read_pairs: u32,
        solid_topology: bool,
    ) -> TransitionCandidate {
        TransitionCandidate {
            node,
            direct_reads,
            read_pairs,
            solid_topology,
        }
    }

    #[test]
    fn unique_solid_topology_is_safe_without_materialized_transition() {
        let candidates = [candidate(7, 0, 0, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(7));
    }

    #[test]
    fn unique_mercy_like_edge_requires_physical_support() {
        let candidates = [candidate(7, 0, 0, false)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), None);
    }

    #[test]
    fn pair_support_can_resolve_an_adjacent_edge_when_reads_do_not() {
        let candidates = [candidate(3, 0, 7, false), candidate(4, 0, 1, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(3));
    }

    #[test]
    fn direct_reads_take_precedence_when_they_clear_the_read_gate() {
        let candidates = [candidate(3, 6, 0, true), candidate(4, 1, 20, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(3));
    }

    #[test]
    fn ambiguous_physical_evidence_remains_unresolved() {
        let candidates = [candidate(3, 0, 5, true), candidate(4, 0, 4, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), None);
    }
}
'''

path.write_text(text)
