use crate::multik_hybrid::{HybridRunStats, LocalKSummary};

/// Apply a conservative post-validation gate to Stage34 hybrid local refinements.
///
/// The local builder already requires zero multi-in/multi-out ambiguity and exact
/// projection through the backbone before setting `resolved`. This additional
/// guard rejects a degenerate high-k collapse where a few retained k-mers remain
/// but no retained transition survives. Such an edgeless graph cannot establish
/// that the original neighborhood was resolved.
#[inline]
pub fn candidate_has_connectivity(candidate: &LocalKSummary) -> bool {
    candidate.retained_kmers > 1 && candidate.directed_edges > 0
}

pub fn enforce_hybrid_safe_gate(run: &mut HybridRunStats) {
    for neighborhood in &mut run.summary.local {
        for candidate in &mut neighborhood.candidates {
            candidate.resolved = candidate.resolved && candidate_has_connectivity(candidate);
        }

        neighborhood.resolved = neighborhood.candidates.iter().any(|candidate| candidate.resolved);
        if let Some(candidate) = neighborhood
            .candidates
            .iter()
            .find(|candidate| candidate.resolved)
        {
            neighborhood.selected_k = Some(candidate.k);
        } else if let Some(candidate) = neighborhood.candidates.last() {
            // Preserve the last attempted k for diagnostics, but keep the
            // neighborhood explicitly unresolved.
            neighborhood.selected_k = Some(candidate.k);
        } else {
            neighborhood.selected_k = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multik_hybrid::{
        HybridNeighborhoodSummary, HybridRunStats, HybridSummary, LocalKSummary,
    };

    fn candidate(k: usize, retained_kmers: usize, directed_edges: usize, resolved: bool) -> LocalKSummary {
        LocalKSummary {
            k,
            retained_kmers,
            directed_edges,
            ambiguous_states: 0,
            dead_end_states: 0,
            projection_missing_nodes: 0,
            projection_missing_edges: 0,
            resolved,
        }
    }

    fn run_with(candidates: Vec<LocalKSummary>) -> HybridRunStats {
        HybridRunStats {
            summary: HybridSummary {
                version: "test".to_string(),
                read_pairs: 1,
                threads: 1,
                backbone_k: 31,
                backbone_retained_kmers: 1,
                backbone_directed_edges: 1,
                backbone_ambiguous_states: 1,
                neighborhoods: 1,
                routed_fragment_copies: 1,
                local_start_k: 55,
                max_local_k: 255,
                backbone_seconds: 0.0,
                routing_seconds: 0.0,
                local_refinement_seconds: 0.0,
                total_seconds: 0.0,
                local: vec![HybridNeighborhoodSummary {
                    id: 1,
                    seed_state: 0,
                    states: 10,
                    ambiguous_states: 1,
                    routed_fragments: 4,
                    supported_max_k: 101,
                    selected_k: Some(101),
                    resolved: true,
                    candidates,
                }],
            },
        }
    }

    #[test]
    fn edgeless_high_k_cannot_resolve_neighborhood() {
        let mut run = run_with(vec![candidate(101, 6, 0, true)]);
        enforce_hybrid_safe_gate(&mut run);
        assert!(!run.summary.local[0].resolved);
        assert!(!run.summary.local[0].candidates[0].resolved);
    }

    #[test]
    fn connected_exact_candidate_remains_resolved() {
        let mut run = run_with(vec![candidate(79, 120, 220, true)]);
        enforce_hybrid_safe_gate(&mut run);
        assert!(run.summary.local[0].resolved);
        assert_eq!(run.summary.local[0].selected_k, Some(79));
    }
}
