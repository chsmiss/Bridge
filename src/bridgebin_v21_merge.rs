use crate::bridgebin::{Assignment, BinSummary, BinningResult, Contig};
use crate::bridgebin_reconcile::MarkerTable;
use crate::bridgebin_v21::{BridgeBinV21Config, PairScoreTable};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug)]
pub struct BridgeBinV21MergeConfig {
    pub min_same: f64,
    pub min_support: usize,
    pub top_support: usize,
}

impl Default for BridgeBinV21MergeConfig {
    fn default() -> Self {
        Self {
            min_same: 0.92,
            min_support: 3,
            top_support: 12,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct BridgeBinV21MergeStats {
    pub input_bins: usize,
    pub output_bins: usize,
    pub candidate_bin_pairs: usize,
    pub accepted_merges: usize,
    pub hard_blocked_bin_pairs: usize,
    pub marker_blocked_bin_pairs: usize,
}

#[derive(Clone, Debug)]
struct Node {
    members: Vec<usize>,
    bp: usize,
    markers: HashSet<String>,
}

#[derive(Clone, Debug)]
struct Candidate {
    left: usize,
    right: usize,
    score: f64,
    support: usize,
}

pub fn merge_bins_v21(
    contigs: &[Contig],
    markers: Option<&MarkerTable>,
    input: BinningResult,
    pair_scores: &PairScoreTable,
    pair_cfg: &BridgeBinV21Config,
    cfg: &BridgeBinV21MergeConfig,
) -> (BinningResult, BridgeBinV21MergeStats) {
    let id_to_contig: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| (contig.id.as_str(), index))
        .collect();
    let mut contig_to_bin: HashMap<usize, usize> = HashMap::new();
    for assignment in &input.assignments {
        let (Some(bin), Some(&contig_index)) = (
            assignment.bin_index,
            id_to_contig.get(assignment.contig_id.as_str()),
        ) else {
            continue;
        };
        contig_to_bin.insert(contig_index, bin);
    }

    let mut members = vec![Vec::new(); input.bins.len()];
    for (&contig_index, &bin) in &contig_to_bin {
        if let Some(target) = members.get_mut(bin) {
            target.push(contig_index);
        }
    }
    let marker_sets: Vec<HashSet<String>> = contigs
        .iter()
        .map(|contig| {
            markers
                .and_then(|table| table.values.get(&contig.id))
                .cloned()
                .unwrap_or_default()
        })
        .collect();
    let nodes: Vec<Node> = members
        .into_iter()
        .map(|members| make_node(members, contigs, &marker_sets))
        .collect();

    let mut stats = BridgeBinV21MergeStats {
        input_bins: nodes.len(),
        ..Default::default()
    };
    if nodes.len() <= 1 || pair_scores.values.is_empty() {
        stats.output_bins = nodes.len();
        return (input, stats);
    }

    let id_to_bin: HashMap<&str, usize> = input
        .assignments
        .iter()
        .filter_map(|assignment| {
            assignment
                .bin_index
                .map(|bin| (assignment.contig_id.as_str(), bin))
        })
        .collect();

    let mut hard: HashSet<(usize, usize)> = HashSet::new();
    for left in 0..nodes.len() {
        for right in (left + 1)..nodes.len() {
            if !nodes[left].markers.is_disjoint(&nodes[right].markers) {
                hard.insert((left, right));
                stats.marker_blocked_bin_pairs += 1;
            }
        }
    }

    let mut support: HashMap<(usize, usize), Vec<(f64, f64)>> = HashMap::new();
    for ((left_id, right_id), pair) in &pair_scores.values {
        if pair.confidence < pair_cfg.min_pair_confidence {
            continue;
        }
        let (Some(&left), Some(&right)) = (
            id_to_bin.get(left_id.as_str()),
            id_to_bin.get(right_id.as_str()),
        ) else {
            continue;
        };
        if left == right {
            continue;
        }
        let key = ordered_index_pair(left, right);
        if pair.same_probability <= pair_cfg.split_max_same {
            if hard.insert(key) {
                stats.hard_blocked_bin_pairs += 1;
            }
            continue;
        }
        if pair.same_probability >= 0.5 {
            support
                .entry(key)
                .or_default()
                .push((pair.same_probability, pair.confidence));
        }
    }

    let mut candidates = Vec::new();
    for ((left, right), values) in support {
        if hard.contains(&(left, right)) {
            continue;
        }
        if let Some((score, count)) = aggregate_support(values, cfg.top_support) {
            if count >= cfg.min_support && score >= cfg.min_same {
                candidates.push(Candidate {
                    left,
                    right,
                    score,
                    support: count,
                });
            }
        }
    }
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| b.support.cmp(&a.support))
            .then_with(|| a.left.cmp(&b.left))
            .then_with(|| a.right.cmp(&b.right))
    });
    stats.candidate_bin_pairs = candidates.len();

    let mut parent: Vec<usize> = (0..nodes.len()).collect();
    let mut component_bins: Vec<HashSet<usize>> = (0..nodes.len())
        .map(|index| HashSet::from([index]))
        .collect();

    for candidate in candidates {
        let left_root = find(&mut parent, candidate.left);
        let right_root = find(&mut parent, candidate.right);
        if left_root == right_root {
            continue;
        }
        if component_conflict(
            &component_bins[left_root],
            &component_bins[right_root],
            &hard,
        ) {
            continue;
        }
        let (keep, remove) = if component_bins[left_root].len() >= component_bins[right_root].len()
        {
            (left_root, right_root)
        } else {
            (right_root, left_root)
        };
        parent[remove] = keep;
        let moved: Vec<usize> = component_bins[remove].drain().collect();
        component_bins[keep].extend(moved);
        stats.accepted_merges += 1;
    }

    let mut root_members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (bin, node) in nodes.iter().enumerate() {
        let root = find(&mut parent, bin);
        root_members
            .entry(root)
            .or_default()
            .extend(node.members.iter().copied());
    }
    let mut merged_nodes: Vec<Node> = root_members
        .into_values()
        .map(|members| make_node(members, contigs, &marker_sets))
        .collect();
    merged_nodes.sort_by_key(|node| std::cmp::Reverse(node.bp));
    stats.output_bins = merged_nodes.len();

    let mut final_bin = vec![None; contigs.len()];
    for (bin, node) in merged_nodes.iter().enumerate() {
        for &contig_index in &node.members {
            final_bin[contig_index] = Some(bin);
        }
    }
    let previous_score: HashMap<&str, f64> = input
        .assignments
        .iter()
        .map(|assignment| (assignment.contig_id.as_str(), assignment.score))
        .collect();
    let assignments = contigs
        .iter()
        .enumerate()
        .map(|(contig_index, contig)| Assignment {
            contig_id: contig.id.clone(),
            bin_index: final_bin[contig_index],
            score: previous_score
                .get(contig.id.as_str())
                .copied()
                .unwrap_or(0.0),
            length: contig.seq.len(),
        })
        .collect();
    let bins = merged_nodes
        .iter()
        .enumerate()
        .map(|(bin_index, node)| BinSummary {
            bin_index,
            contig_count: node.members.len(),
            total_bp: node.bp,
            mean_gc: mean_gc(node, contigs),
        })
        .collect();
    (BinningResult { assignments, bins }, stats)
}

fn aggregate_support(mut values: Vec<(f64, f64)>, top_support: usize) -> Option<(f64, usize)> {
    values.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
    values.truncate(top_support.max(1));
    if values.is_empty() {
        return None;
    }
    let total_weight = values.iter().map(|(_, confidence)| confidence).sum::<f64>();
    if total_weight <= 0.0 {
        return None;
    }
    let score = values
        .iter()
        .map(|(probability, confidence)| probability * confidence)
        .sum::<f64>()
        / total_weight;
    Some((score.clamp(0.0, 1.0), values.len()))
}

fn component_conflict(
    left: &HashSet<usize>,
    right: &HashSet<usize>,
    hard: &HashSet<(usize, usize)>,
) -> bool {
    left.iter().any(|a| {
        right
            .iter()
            .any(|b| hard.contains(&ordered_index_pair(*a, *b)))
    })
}

fn find(parent: &mut [usize], mut index: usize) -> usize {
    let mut root = index;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[index] != index {
        let next = parent[index];
        parent[index] = root;
        index = next;
    }
    root
}

fn ordered_index_pair(left: usize, right: usize) -> (usize, usize) {
    (left.min(right), left.max(right))
}

fn make_node(
    members: Vec<usize>,
    contigs: &[Contig],
    marker_sets: &[HashSet<String>],
) -> Node {
    let bp = members.iter().map(|index| contigs[*index].seq.len()).sum();
    let mut markers = HashSet::new();
    for &member in &members {
        markers.extend(marker_sets[member].iter().cloned());
    }
    Node {
        members,
        bp,
        markers,
    }
}

fn mean_gc(node: &Node, contigs: &[Contig]) -> f64 {
    let mut gc = 0usize;
    let mut valid = 0usize;
    for &member in &node.members {
        for &base in &contigs[member].seq {
            match base.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc += 1;
                    valid += 1;
                }
                b'A' | b'T' => valid += 1,
                _ => {}
            }
        }
    }
    if valid == 0 {
        0.0
    } else {
        gc as f64 / valid as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridgebin_v21::PairScore;

    fn contig(id: &str) -> Contig {
        Contig {
            id: id.to_string(),
            seq: b"ACGT".repeat(100),
        }
    }

    fn two_bins(contigs: &[Contig]) -> BinningResult {
        BinningResult {
            assignments: contigs
                .iter()
                .enumerate()
                .map(|(index, contig)| Assignment {
                    contig_id: contig.id.clone(),
                    bin_index: Some(index / 2),
                    score: 1.0,
                    length: contig.seq.len(),
                })
                .collect(),
            bins: vec![
                BinSummary {
                    bin_index: 0,
                    contig_count: 2,
                    total_bp: 800,
                    mean_gc: 0.5,
                },
                BinSummary {
                    bin_index: 1,
                    contig_count: 2,
                    total_bp: 800,
                    mean_gc: 0.5,
                },
            ],
        }
    }

    fn pair(left: &str, right: &str, same: f64) -> ((String, String), PairScore) {
        let key = if left <= right {
            (left.to_string(), right.to_string())
        } else {
            (right.to_string(), left.to_string())
        };
        (
            key,
            PairScore {
                same_probability: same,
                confidence: 0.99,
                source: "test".to_string(),
            },
        )
    }

    #[test]
    fn multiple_biological_links_merge_fragmented_bins() {
        let contigs = vec![contig("a1"), contig("a2"), contig("a3"), contig("a4")];
        let scores = PairScoreTable {
            values: HashMap::from([
                pair("a1", "a3", 0.99),
                pair("a1", "a4", 0.98),
                pair("a2", "a3", 0.97),
            ]),
        };
        let cfg = BridgeBinV21MergeConfig {
            min_support: 3,
            ..Default::default()
        };
        let (result, stats) = merge_bins_v21(
            &contigs,
            None,
            two_bins(&contigs),
            &scores,
            &BridgeBinV21Config::default(),
            &cfg,
        );
        assert_eq!(stats.accepted_merges, 1);
        assert_eq!(result.bins.len(), 1);
    }

    #[test]
    fn one_hard_negative_blocks_transitive_bin_merge() {
        let contigs = vec![contig("a1"), contig("a2"), contig("b1"), contig("b2")];
        let scores = PairScoreTable {
            values: HashMap::from([
                pair("a1", "b1", 0.99),
                pair("a1", "b2", 0.99),
                pair("a2", "b1", 0.99),
                pair("a2", "b2", 0.01),
            ]),
        };
        let cfg = BridgeBinV21MergeConfig {
            min_support: 3,
            ..Default::default()
        };
        let (result, stats) = merge_bins_v21(
            &contigs,
            None,
            two_bins(&contigs),
            &scores,
            &BridgeBinV21Config::default(),
            &cfg,
        );
        assert_eq!(stats.accepted_merges, 0);
        assert_eq!(result.bins.len(), 2);
    }
}
