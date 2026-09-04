// BridgeBin v2.1 pair refinement wrapper.
//
// The original implementation is retained verbatim in bridgebin_v21_legacy.rs.  Before
// delegating to it, this wrapper interprets calibrated hard negatives as *cut evidence*
// over the current v2 bin instead of forcing a singleton rebuild that needs dense positive
// pair scores.  This preserves recall when the Biological Brain is sparse but precise.
mod legacy {
    include!("bridgebin_v21_legacy.rs");
}

pub use legacy::{
    read_pair_score_table, BridgeBinV21Config, BridgeBinV21Stats, PairScore, PairScoreTable,
};

use crate::bridgebin::{BinningResult, Contig};
use crate::bridgebin_reconcile::MarkerTable;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, Default)]
struct PrepartitionStats {
    conflicted_input_bins: usize,
    split_bins: usize,
    hard_negative_pairs: usize,
    marker_negative_pairs: usize,
    positive_pairs: usize,
}

#[derive(Clone, Copy, Debug)]
struct WeightedEdge {
    left: usize,
    right: usize,
    weight: f64,
}

#[derive(Clone, Debug)]
struct CutCandidate {
    left: Vec<usize>,
    right: Vec<usize>,
    objective: f64,
    cut_model_hard: usize,
    cut_marker_hard: usize,
    total_model_hard: usize,
    min_side_bp: usize,
}

pub fn refine_bins_v21(
    contigs: &[Contig],
    markers: Option<&MarkerTable>,
    initial: BinningResult,
    pair_scores: &PairScoreTable,
    cfg: &BridgeBinV21Config,
) -> (BinningResult, BridgeBinV21Stats) {
    let original_input_bins = initial
        .assignments
        .iter()
        .filter_map(|assignment| assignment.bin_index)
        .collect::<HashSet<_>>()
        .len();

    let (prepartitioned, filtered_scores, prep) =
        prepartition_initial(contigs, markers, initial, pair_scores, cfg);
    let (result, mut stats) =
        legacy::refine_bins_v21(contigs, markers, prepartitioned, &filtered_scores, cfg);

    // The legacy kernel only sees already separated sub-bins.  Restore statistics to the
    // user's original v2 input so logs remain interpretable and comparable across versions.
    stats.input_bins = original_input_bins;
    stats.conflicted_input_bins += prep.conflicted_input_bins;
    stats.split_bins += prep.split_bins;
    stats.hard_negative_pairs += prep.hard_negative_pairs;
    stats.marker_negative_pairs += prep.marker_negative_pairs;
    stats.positive_pairs += prep.positive_pairs;
    (result, stats)
}

fn prepartition_initial(
    contigs: &[Contig],
    markers: Option<&MarkerTable>,
    mut initial: BinningResult,
    pair_scores: &PairScoreTable,
    cfg: &BridgeBinV21Config,
) -> (BinningResult, PairScoreTable, PrepartitionStats) {
    let by_id: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| (contig.id.as_str(), index))
        .collect();
    let marker_sets: Vec<HashSet<String>> = contigs
        .iter()
        .map(|contig| {
            markers
                .and_then(|table| table.values.get(&contig.id))
                .cloned()
                .unwrap_or_default()
        })
        .collect();

    let mut indexed_scores: HashMap<(usize, usize), PairScore> = HashMap::new();
    for ((left, right), score) in &pair_scores.values {
        let (Some(&a), Some(&b)) = (by_id.get(left.as_str()), by_id.get(right.as_str())) else {
            continue;
        };
        indexed_scores.insert((a.min(b), a.max(b)), score.clone());
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for (contig_index, assignment) in initial.assignments.iter().enumerate() {
        if let Some(bin_index) = assignment.bin_index {
            groups.entry(bin_index).or_default().push(contig_index);
        }
    }
    let mut group_ids: Vec<usize> = groups.keys().copied().collect();
    group_ids.sort_unstable();

    let mut prep = PrepartitionStats::default();
    let mut new_bin_of = vec![None; contigs.len()];
    let mut next_bin = 0usize;

    for group_id in group_ids {
        let members = groups.remove(&group_id).unwrap_or_default();
        let member_set: HashSet<usize> = members.iter().copied().collect();
        let mut model_hard: HashMap<(usize, usize), f64> = HashMap::new();
        let mut positive: HashMap<(usize, usize), f64> = HashMap::new();

        for (&pair, score) in &indexed_scores {
            if !member_set.contains(&pair.0) || !member_set.contains(&pair.1) {
                continue;
            }
            if score.confidence >= cfg.min_pair_confidence
                && score.same_probability <= cfg.split_max_same
            {
                model_hard.insert(pair, hard_weight(score, cfg));
            }
            if score.confidence >= cfg.min_pair_confidence
                && score.same_probability >= cfg.join_min_same
            {
                positive.insert(pair, positive_weight(score, cfg));
            }
        }
        prep.hard_negative_pairs += model_hard.len();
        prep.positive_pairs += positive.len();

        let mut marker_hard: HashSet<(usize, usize)> = HashSet::new();
        let mut marker_to_members: HashMap<&str, Vec<usize>> = HashMap::new();
        for &member in &members {
            for marker in &marker_sets[member] {
                marker_to_members
                    .entry(marker.as_str())
                    .or_default()
                    .push(member);
            }
        }
        for marker_members in marker_to_members.values() {
            for left_pos in 0..marker_members.len() {
                for right_pos in (left_pos + 1)..marker_members.len() {
                    let left = marker_members[left_pos];
                    let right = marker_members[right_pos];
                    let pair = (left.min(right), left.max(right));
                    if marker_hard.insert(pair) && !model_hard.contains_key(&pair) {
                        prep.marker_negative_pairs += 1;
                    }
                }
            }
        }

        if model_hard.is_empty() && marker_hard.is_empty() {
            for member in members {
                new_bin_of[member] = Some(next_bin);
            }
            next_bin += 1;
            continue;
        }
        prep.conflicted_input_bins += 1;

        let mut pieces = recursive_signed_partition(
            &members,
            &model_hard,
            &marker_hard,
            &positive,
            contigs,
            cfg,
            0,
        );
        pieces.sort_by(|a, b| {
            piece_bp(b, contigs)
                .cmp(&piece_bp(a, contigs))
                .then_with(|| a.first().cmp(&b.first()))
        });

        let actually_split = pieces.len() > 1;
        let mut retained = 0usize;
        for piece in pieces {
            let bp = piece_bp(&piece, contigs);
            if actually_split && bp < cfg.min_subbin_bp {
                continue;
            }
            for member in piece {
                new_bin_of[member] = Some(next_bin);
            }
            next_bin += 1;
            retained += 1;
        }
        if retained > 1 {
            prep.split_bins += 1;
        }
    }

    for (contig_index, assignment) in initial.assignments.iter_mut().enumerate() {
        assignment.bin_index = new_bin_of[contig_index];
    }

    // Once a globally coherent cut has been selected, a small number of model hard
    // negatives that remain inside one retained side are treated as calibration noise.
    // Removing only these internal hard edges prevents the legacy singleton splitter from
    // undoing the robust partition.  Cross-side hard negatives remain available to rescue
    // vetoes and the later cross-bin merge stage.
    let mut filtered_scores = pair_scores.clone();
    filtered_scores.values.retain(|(left, right), score| {
        if score.confidence < cfg.min_pair_confidence || score.same_probability > cfg.split_max_same
        {
            return true;
        }
        let (Some(&a), Some(&b)) = (by_id.get(left.as_str()), by_id.get(right.as_str())) else {
            return true;
        };
        match (new_bin_of[a], new_bin_of[b]) {
            (Some(left_bin), Some(right_bin)) if left_bin == right_bin => false,
            _ => true,
        }
    });

    (initial, filtered_scores, prep)
}

fn recursive_signed_partition(
    members: &[usize],
    model_hard: &HashMap<(usize, usize), f64>,
    marker_hard: &HashSet<(usize, usize)>,
    positive: &HashMap<(usize, usize), f64>,
    contigs: &[Contig],
    cfg: &BridgeBinV21Config,
    depth: usize,
) -> Vec<Vec<usize>> {
    if members.len() < 2 || depth >= 16 {
        return vec![members.to_vec()];
    }
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let local_model = model_hard
        .iter()
        .filter(|(pair, _)| member_set.contains(&pair.0) && member_set.contains(&pair.1))
        .count();
    let local_marker = marker_hard
        .iter()
        .filter(|pair| member_set.contains(&pair.0) && member_set.contains(&pair.1))
        .count();
    if local_model + local_marker == 0 {
        return vec![members.to_vec()];
    }

    let Some(candidate) = best_binary_cut(members, model_hard, marker_hard, positive, contigs)
    else {
        return vec![members.to_vec()];
    };
    if !accept_cut(&candidate, cfg) {
        return vec![members.to_vec()];
    }

    let mut pieces = recursive_signed_partition(
        &candidate.left,
        model_hard,
        marker_hard,
        positive,
        contigs,
        cfg,
        depth + 1,
    );
    pieces.extend(recursive_signed_partition(
        &candidate.right,
        model_hard,
        marker_hard,
        positive,
        contigs,
        cfg,
        depth + 1,
    ));
    pieces
}

fn accept_cut(candidate: &CutCandidate, cfg: &BridgeBinV21Config) -> bool {
    if candidate.left.is_empty() || candidate.right.is_empty() {
        return false;
    }
    if candidate.cut_marker_hard > 0 {
        // Single-copy marker duplication is a strict biological veto.  Recursive cuts will
        // keep separating a marker-conflict component until no duplicated marker remains.
        return true;
    }
    if candidate.total_model_hard == 0 {
        return false;
    }
    let explained = candidate.cut_model_hard as f64 / candidate.total_model_hard as f64;
    if explained < 0.70 || candidate.cut_model_hard < cfg.min_pair_support.max(1) {
        return false;
    }
    // Do not eject a tiny component on only one or two noisy model negatives.  Small
    // contamination can still be removed when at least three independent hard negatives
    // support it; two large genome components only need the normal min_pair_support.
    if candidate.min_side_bp < cfg.min_subbin_bp
        && candidate.cut_model_hard < cfg.min_pair_support.max(3)
    {
        return false;
    }
    true
}

fn best_binary_cut(
    members: &[usize],
    model_hard: &HashMap<(usize, usize), f64>,
    marker_hard: &HashSet<(usize, usize)>,
    positive: &HashMap<(usize, usize), f64>,
    contigs: &[Contig],
) -> Option<CutCandidate> {
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let mut diff_weights: HashMap<(usize, usize), f64> = HashMap::new();
    for (&pair, &weight) in model_hard {
        if member_set.contains(&pair.0) && member_set.contains(&pair.1) {
            *diff_weights.entry(pair).or_insert(0.0) += weight;
        }
    }
    for &pair in marker_hard {
        if member_set.contains(&pair.0) && member_set.contains(&pair.1) {
            *diff_weights.entry(pair).or_insert(0.0) += 4.0;
        }
    }
    if diff_weights.is_empty() {
        return None;
    }
    let same_edges: Vec<WeightedEdge> = positive
        .iter()
        .filter_map(|(&pair, &weight)| {
            (member_set.contains(&pair.0) && member_set.contains(&pair.1)).then_some(WeightedEdge {
                left: pair.0,
                right: pair.1,
                weight,
            })
        })
        .collect();
    let mut diff_edges: Vec<WeightedEdge> = diff_weights
        .iter()
        .map(|(&(left, right), &weight)| WeightedEdge {
            left,
            right,
            weight,
        })
        .collect();
    diff_edges.sort_by(|a, b| {
        b.weight
            .partial_cmp(&a.weight)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.left.cmp(&b.left))
            .then_with(|| a.right.cmp(&b.right))
    });

    let mut best: Option<CutCandidate> = None;
    for seed in diff_edges.iter().take(8) {
        let labels = greedy_signed_cut(members, seed, &diff_edges, &same_edges, contigs);
        let candidate = evaluate_cut(
            members,
            &labels,
            model_hard,
            marker_hard,
            &diff_edges,
            &same_edges,
            contigs,
        );
        let replace = match &best {
            None => true,
            Some(existing) => {
                candidate.objective > existing.objective + 1e-9
                    || ((candidate.objective - existing.objective).abs() <= 1e-9
                        && (
                            candidate.cut_marker_hard,
                            candidate.cut_model_hard,
                            candidate.min_side_bp,
                        ) > (
                            existing.cut_marker_hard,
                            existing.cut_model_hard,
                            existing.min_side_bp,
                        ))
            }
        };
        if replace {
            best = Some(candidate);
        }
    }
    best
}

fn greedy_signed_cut(
    members: &[usize],
    seed: &WeightedEdge,
    diff_edges: &[WeightedEdge],
    same_edges: &[WeightedEdge],
    contigs: &[Contig],
) -> HashMap<usize, u8> {
    let mut labels = HashMap::from([(seed.left, 0u8), (seed.right, 1u8)]);
    let mut unassigned: HashSet<usize> = members.iter().copied().collect();
    unassigned.remove(&seed.left);
    unassigned.remove(&seed.right);

    while !unassigned.is_empty() {
        let mut best_member: Option<usize> = None;
        let mut best_evidence = -1.0f64;
        let mut best_margin = -1.0f64;
        let mut best_scores = [0.0f64; 2];
        for &member in &unassigned {
            let scores = assignment_scores(member, &labels, diff_edges, same_edges);
            let evidence = scores[0] + scores[1];
            let margin = (scores[0] - scores[1]).abs();
            let replace = evidence > best_evidence + 1e-9
                || ((evidence - best_evidence).abs() <= 1e-9
                    && (margin > best_margin + 1e-9
                        || ((margin - best_margin).abs() <= 1e-9
                            && best_member.is_none_or(|old| {
                                contigs[member].seq.len() > contigs[old].seq.len()
                                    || (contigs[member].seq.len() == contigs[old].seq.len()
                                        && member < old)
                            }))));
            if replace {
                best_member = Some(member);
                best_evidence = evidence;
                best_margin = margin;
                best_scores = scores;
            }
        }
        let member = best_member.expect("unassigned member");
        let side = if (best_scores[0] - best_scores[1]).abs() > 1e-9 {
            usize::from(best_scores[1] > best_scores[0])
        } else {
            let bp0 = labeled_bp(&labels, 0, contigs);
            let bp1 = labeled_bp(&labels, 1, contigs);
            usize::from(bp1 > bp0)
        } as u8;
        labels.insert(member, side);
        unassigned.remove(&member);
    }

    // Deterministic 1-opt cleanup.  It corrects seed orientation mistakes while keeping
    // both sides non-empty and converges quickly on the sparse candidate graph.
    for _ in 0..8 {
        let mut changed = false;
        let mut order = members.to_vec();
        order.sort_by(|a, b| {
            contigs[*b]
                .seq
                .len()
                .cmp(&contigs[*a].seq.len())
                .then_with(|| a.cmp(b))
        });
        for member in order {
            let current = labels[&member];
            let current_count = labels.values().filter(|&&side| side == current).count();
            if current_count <= 1 {
                continue;
            }
            let before = cut_objective(&labels, diff_edges, same_edges);
            labels.insert(member, 1 - current);
            let after = cut_objective(&labels, diff_edges, same_edges);
            if after > before + 1e-9 {
                changed = true;
            } else {
                labels.insert(member, current);
            }
        }
        if !changed {
            break;
        }
    }
    labels
}

fn assignment_scores(
    member: usize,
    labels: &HashMap<usize, u8>,
    diff_edges: &[WeightedEdge],
    same_edges: &[WeightedEdge],
) -> [f64; 2] {
    let mut scores = [0.0f64; 2];
    for edge in diff_edges {
        let other = if edge.left == member {
            edge.right
        } else if edge.right == member {
            edge.left
        } else {
            continue;
        };
        if let Some(&label) = labels.get(&other) {
            scores[usize::from(1 - label)] += edge.weight;
        }
    }
    for edge in same_edges {
        let other = if edge.left == member {
            edge.right
        } else if edge.right == member {
            edge.left
        } else {
            continue;
        };
        if let Some(&label) = labels.get(&other) {
            scores[usize::from(label)] += edge.weight;
        }
    }
    scores
}

fn evaluate_cut(
    members: &[usize],
    labels: &HashMap<usize, u8>,
    model_hard: &HashMap<(usize, usize), f64>,
    marker_hard: &HashSet<(usize, usize)>,
    diff_edges: &[WeightedEdge],
    same_edges: &[WeightedEdge],
    contigs: &[Contig],
) -> CutCandidate {
    let mut left = Vec::new();
    let mut right = Vec::new();
    for &member in members {
        if labels.get(&member).copied().unwrap_or(0) == 0 {
            left.push(member);
        } else {
            right.push(member);
        }
    }
    let cut_model_hard = model_hard
        .keys()
        .filter(|&&(a, b)| labels.get(&a) != labels.get(&b))
        .count();
    let cut_marker_hard = marker_hard
        .iter()
        .filter(|&&(a, b)| labels.get(&a) != labels.get(&b))
        .count();
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let total_model_hard = model_hard
        .keys()
        .filter(|&&(a, b)| member_set.contains(&a) && member_set.contains(&b))
        .count();
    CutCandidate {
        min_side_bp: piece_bp(&left, contigs).min(piece_bp(&right, contigs)),
        objective: cut_objective(labels, diff_edges, same_edges),
        left,
        right,
        cut_model_hard,
        cut_marker_hard,
        total_model_hard,
    }
}

fn cut_objective(
    labels: &HashMap<usize, u8>,
    diff_edges: &[WeightedEdge],
    same_edges: &[WeightedEdge],
) -> f64 {
    let diff_score = diff_edges
        .iter()
        .filter(|edge| labels.get(&edge.left) != labels.get(&edge.right))
        .map(|edge| edge.weight)
        .sum::<f64>();
    let same_score = same_edges
        .iter()
        .filter(|edge| labels.get(&edge.left) == labels.get(&edge.right))
        .map(|edge| edge.weight)
        .sum::<f64>();
    diff_score + same_score
}

fn labeled_bp(labels: &HashMap<usize, u8>, side: u8, contigs: &[Contig]) -> usize {
    labels
        .iter()
        .filter(|(_, label)| **label == side)
        .map(|(&member, _)| contigs[member].seq.len())
        .sum()
}

fn piece_bp(members: &[usize], contigs: &[Contig]) -> usize {
    members
        .iter()
        .map(|&member| contigs[member].seq.len())
        .sum()
}

fn hard_weight(score: &PairScore, cfg: &BridgeBinV21Config) -> f64 {
    let denominator = cfg.split_max_same.max(1e-9);
    let margin = ((cfg.split_max_same - score.same_probability) / denominator).clamp(0.0, 1.0);
    (1.0 + margin) * (0.5 + 0.5 * score.confidence)
}

fn positive_weight(score: &PairScore, cfg: &BridgeBinV21Config) -> f64 {
    let denominator = (1.0 - cfg.join_min_same).max(1e-9);
    let margin = ((score.same_probability - cfg.join_min_same) / denominator).clamp(0.0, 1.0);
    (0.75 + margin) * (0.5 + 0.5 * score.confidence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridgebin::{Assignment, BinSummary};

    fn contig(id: &str, bp: usize) -> Contig {
        Contig {
            id: id.to_string(),
            seq: vec![b'A'; bp],
        }
    }

    fn one_bin(contigs: &[Contig]) -> BinningResult {
        BinningResult {
            assignments: contigs
                .iter()
                .map(|contig| Assignment {
                    contig_id: contig.id.clone(),
                    bin_index: Some(0),
                    score: 1.0,
                    length: contig.seq.len(),
                })
                .collect(),
            bins: vec![BinSummary {
                bin_index: 0,
                contig_count: contigs.len(),
                total_bp: contigs.iter().map(|contig| contig.seq.len()).sum(),
                mean_gc: 0.0,
            }],
        }
    }

    fn insert_score(table: &mut PairScoreTable, left: &str, right: &str, same: f64) {
        let key = if left <= right {
            (left.to_string(), right.to_string())
        } else {
            (right.to_string(), left.to_string())
        };
        table.values.insert(
            key,
            PairScore {
                same_probability: same,
                confidence: 1.0,
                source: "test".to_string(),
            },
        );
    }

    #[test]
    fn sparse_positive_edges_do_not_destroy_conflict_split_recall() {
        let contigs = vec![
            contig("a1", 10_000),
            contig("a2", 10_000),
            contig("a3", 10_000),
            contig("b1", 10_000),
            contig("b2", 10_000),
            contig("b3", 10_000),
            contig("b4", 10_000),
        ];
        let mut scores = PairScoreTable::default();
        for a in ["a1", "a2", "a3"] {
            for b in ["b1", "b2", "b3", "b4"] {
                insert_score(&mut scores, a, b, 0.02);
            }
        }
        // Two realistic false hard negatives inside one genome should be outvoted by the
        // coherent cross-genome cut rather than causing singleton fragmentation.
        insert_score(&mut scores, "b1", "b3", 0.025);
        insert_score(&mut scores, "b2", "b3", 0.030);
        // Only two positive edges: deliberately far too sparse for the legacy rebuild.
        insert_score(&mut scores, "a1", "a2", 0.60);
        insert_score(&mut scores, "b1", "b2", 0.60);

        let cfg = BridgeBinV21Config {
            min_pair_confidence: 0.0,
            split_max_same: 0.075,
            join_min_same: 0.37,
            rescue_min_same: 0.37,
            min_pair_support: 2,
            min_subbin_bp: 20_000,
            ..Default::default()
        };
        let (result, stats) = refine_bins_v21(&contigs, None, one_bin(&contigs), &scores, &cfg);
        let bins: HashMap<&str, Option<usize>> = result
            .assignments
            .iter()
            .map(|assignment| (assignment.contig_id.as_str(), assignment.bin_index))
            .collect();

        assert_eq!(stats.split_bins, 1);
        assert!(bins.values().all(Option::is_some));
        assert_eq!(bins["a1"], bins["a2"]);
        assert_eq!(bins["a2"], bins["a3"]);
        assert_eq!(bins["b1"], bins["b2"]);
        assert_eq!(bins["b2"], bins["b3"]);
        assert_eq!(bins["b3"], bins["b4"]);
        assert_ne!(bins["a1"], bins["b1"]);
    }
}
