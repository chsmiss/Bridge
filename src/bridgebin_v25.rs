use crate::bridgebin::{BinSummary, BinningResult, Contig, CoverageTable};
use crate::bridgebin_reconcile::MarkerTable;
use crate::bridgebin_v21::{refine_bins_v21, BridgeBinV21Config, BridgeBinV21Stats, PairScore, PairScoreTable};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, Default)]
pub struct BridgeBinV25Stats {
    pub input_bins: usize,
    pub biological_conflict_bins: usize,
    pub accepted_seed_splits: usize,
    pub biological_anchor_contigs: usize,
    pub propagated_contigs: usize,
    pub ambiguous_propagations: usize,
    pub v21: BridgeBinV21Stats,
}

#[derive(Clone, Debug)]
struct CheapFeature {
    gc: f64,
    kmer: [f64; 1024],
    coverage: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
struct SignedEdge {
    left: usize,
    right: usize,
    weight: f64,
    hard: bool,
}

#[derive(Clone, Debug)]
struct SeedCut {
    labels: HashMap<usize, u8>,
    cut_hard: usize,
    total_hard: usize,
    objective: f64,
}

pub fn refine_bins_v25(
    contigs: &[Contig],
    coverage: Option<&CoverageTable>,
    markers: Option<&MarkerTable>,
    initial: BinningResult,
    pair_scores: &PairScoreTable,
    cfg: &BridgeBinV21Config,
) -> (BinningResult, PairScoreTable, BridgeBinV25Stats) {
    let mut stats = BridgeBinV25Stats::default();
    let by_id: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| (contig.id.as_str(), index))
        .collect();
    let indexed = index_scores(pair_scores, &by_id);
    let cheap: Vec<CheapFeature> = contigs
        .iter()
        .map(|contig| cheap_feature(contig, coverage))
        .collect();

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for (contig_index, assignment) in initial.assignments.iter().enumerate() {
        if let Some(bin_index) = assignment.bin_index {
            groups.entry(bin_index).or_default().push(contig_index);
        }
    }
    stats.input_bins = groups.len();
    let mut group_ids: Vec<usize> = groups.keys().copied().collect();
    group_ids.sort_unstable();

    let mut new_bin_of = vec![None; contigs.len()];
    let mut next_bin = 0usize;
    let mut consumed_internal_pairs: HashSet<(usize, usize)> = HashSet::new();

    for group_id in group_ids {
        let members = groups.remove(&group_id).unwrap_or_default();
        let member_set: HashSet<usize> = members.iter().copied().collect();
        let mut edges = Vec::new();
        let mut anchors = HashSet::new();
        let mut hard_pairs = 0usize;

        for (&pair, score) in &indexed {
            if !member_set.contains(&pair.0) || !member_set.contains(&pair.1) {
                continue;
            }
            if score.confidence < cfg.min_pair_confidence {
                continue;
            }
            if score.same_probability <= cfg.split_max_same {
                edges.push(SignedEdge {
                    left: pair.0,
                    right: pair.1,
                    weight: hard_weight(score, cfg),
                    hard: true,
                });
                anchors.insert(pair.0);
                anchors.insert(pair.1);
                hard_pairs += 1;
            } else if score.same_probability >= cfg.join_min_same {
                edges.push(SignedEdge {
                    left: pair.0,
                    right: pair.1,
                    weight: positive_weight(score, cfg),
                    hard: false,
                });
                anchors.insert(pair.0);
                anchors.insert(pair.1);
            }
        }

        if hard_pairs < cfg.min_pair_support.max(1) || anchors.len() < 4 {
            assign_one_bin(&members, &mut new_bin_of, &mut next_bin);
            continue;
        }
        stats.biological_conflict_bins += 1;
        stats.biological_anchor_contigs += anchors.len();

        let mut anchor_list: Vec<usize> = anchors.into_iter().collect();
        anchor_list.sort_unstable();
        let Some(seed_cut) = best_seed_cut(&anchor_list, &edges, contigs) else {
            assign_one_bin(&members, &mut new_bin_of, &mut next_bin);
            continue;
        };
        let explained = seed_cut.cut_hard as f64 / seed_cut.total_hard.max(1) as f64;
        let side0_anchors = seed_cut.labels.values().filter(|&&side| side == 0).count();
        let side1_anchors = seed_cut.labels.values().filter(|&&side| side == 1).count();
        if explained < 0.70
            || seed_cut.cut_hard < cfg.min_pair_support.max(2)
            || side0_anchors < 2
            || side1_anchors < 2
        {
            assign_one_bin(&members, &mut new_bin_of, &mut next_bin);
            continue;
        }

        let mut labels = seed_cut.labels;
        let biological_anchors: HashSet<usize> = labels.keys().copied().collect();
        let mut side_bp = [0usize; 2];
        for (&member, &side) in &labels {
            side_bp[usize::from(side)] += contigs[member].seq.len();
        }

        let mut unresolved: Vec<usize> = members
            .iter()
            .copied()
            .filter(|member| !labels.contains_key(member))
            .collect();
        // Resolve the most decisive cheap assignments first. Only biological anchors are
        // used as references, so propagated members can never bootstrap a false identity.
        while !unresolved.is_empty() {
            let mut best_position = 0usize;
            let mut best_margin = f64::NEG_INFINITY;
            let mut best_scores = [0.0f64; 2];
            for (position, &member) in unresolved.iter().enumerate() {
                let scores = propagation_scores(
                    member,
                    &biological_anchors,
                    &labels,
                    &cheap,
                    4,
                );
                let margin = (scores[0] - scores[1]).abs();
                if margin > best_margin + 1e-12
                    || ((margin - best_margin).abs() <= 1e-12
                        && contigs[member].seq.len() > contigs[unresolved[best_position]].seq.len())
                {
                    best_position = position;
                    best_margin = margin;
                    best_scores = scores;
                }
            }
            let member = unresolved.swap_remove(best_position);
            let side = if (best_scores[0] - best_scores[1]).abs() > 1e-6 {
                usize::from(best_scores[1] > best_scores[0])
            } else {
                stats.ambiguous_propagations += 1;
                usize::from(side_bp[0] > side_bp[1])
            } as u8;
            labels.insert(member, side);
            side_bp[usize::from(side)] += contigs[member].seq.len();
            stats.propagated_contigs += 1;
        }

        if side_bp[0] < cfg.min_subbin_bp || side_bp[1] < cfg.min_subbin_bp {
            assign_one_bin(&members, &mut new_bin_of, &mut next_bin);
            continue;
        }

        let left_bin = next_bin;
        let right_bin = next_bin + 1;
        next_bin += 2;
        for member in members {
            new_bin_of[member] = Some(if labels.get(&member).copied().unwrap_or(0) == 0 {
                left_bin
            } else {
                right_bin
            });
        }
        stats.accepted_seed_splits += 1;

        // Internal low-p_same edges were consumed to choose the seed split. Remove only
        // those that now lie inside a retained side so legacy v2.1 cannot recursively
        // fragment a side on calibration noise. Cross-side negatives remain for merge veto.
        for (&pair, score) in &indexed {
            if score.confidence < cfg.min_pair_confidence
                || score.same_probability > cfg.split_max_same
            {
                continue;
            }
            if matches!(
                (new_bin_of[pair.0], new_bin_of[pair.1]),
                (Some(a), Some(b)) if a == b
            ) {
                consumed_internal_pairs.insert(pair);
            }
        }
    }

    let mut prepartitioned = initial;
    for (contig_index, assignment) in prepartitioned.assignments.iter_mut().enumerate() {
        if assignment.bin_index.is_some() {
            assignment.bin_index = new_bin_of[contig_index];
        }
    }
    prepartitioned.bins = summarize_bins(contigs, &prepartitioned);

    let mut filtered_scores = pair_scores.clone();
    filtered_scores.values.retain(|(left, right), _| {
        let (Some(&a), Some(&b)) = (by_id.get(left.as_str()), by_id.get(right.as_str())) else {
            return true;
        };
        !consumed_internal_pairs.contains(&(a.min(b), a.max(b)))
    });

    let (result, v21) = refine_bins_v21(
        contigs,
        markers,
        prepartitioned,
        &filtered_scores,
        cfg,
    );
    stats.v21 = v21;
    (result, filtered_scores, stats)
}

fn assign_one_bin(members: &[usize], new_bin_of: &mut [Option<usize>], next_bin: &mut usize) {
    for &member in members {
        new_bin_of[member] = Some(*next_bin);
    }
    *next_bin += 1;
}

fn index_scores(
    scores: &PairScoreTable,
    by_id: &HashMap<&str, usize>,
) -> HashMap<(usize, usize), PairScore> {
    let mut out = HashMap::new();
    for ((left, right), score) in &scores.values {
        let (Some(&a), Some(&b)) = (by_id.get(left.as_str()), by_id.get(right.as_str())) else {
            continue;
        };
        out.insert((a.min(b), a.max(b)), score.clone());
    }
    out
}

fn hard_weight(score: &PairScore, cfg: &BridgeBinV21Config) -> f64 {
    let scale = cfg.split_max_same.max(1e-6);
    1.0 + ((cfg.split_max_same - score.same_probability) / scale).clamp(0.0, 1.0)
}

fn positive_weight(score: &PairScore, cfg: &BridgeBinV21Config) -> f64 {
    let scale = (1.0 - cfg.join_min_same).max(1e-6);
    1.0 + ((score.same_probability - cfg.join_min_same) / scale).clamp(0.0, 1.0)
}

fn best_seed_cut(
    anchors: &[usize],
    edges: &[SignedEdge],
    contigs: &[Contig],
) -> Option<SeedCut> {
    if anchors.len() < 2 {
        return None;
    }
    if anchors.len() <= 18 {
        return exact_seed_cut(anchors, edges, contigs);
    }
    greedy_seed_cut(anchors, edges, contigs)
}

fn exact_seed_cut(
    anchors: &[usize],
    edges: &[SignedEdge],
    contigs: &[Contig],
) -> Option<SeedCut> {
    let positions: HashMap<usize, usize> = anchors
        .iter()
        .copied()
        .enumerate()
        .map(|(position, member)| (member, position))
        .collect();
    let total_hard = edges.iter().filter(|edge| edge.hard).count();
    let total_masks = 1usize << anchors.len();
    let mut best: Option<SeedCut> = None;

    for mask in 1usize..total_masks {
        if mask & 1 != 0 {
            continue;
        }
        let right_count = mask.count_ones() as usize;
        if right_count == 0 || right_count == anchors.len() {
            continue;
        }
        let mut labels = HashMap::new();
        let mut side_bp = [0usize; 2];
        for (position, &member) in anchors.iter().enumerate() {
            let side = ((mask >> position) & 1) as u8;
            labels.insert(member, side);
            side_bp[usize::from(side)] += contigs[member].seq.len();
        }
        let (objective, cut_hard) = score_cut(&labels, edges);
        let candidate = SeedCut {
            labels,
            cut_hard,
            total_hard,
            objective,
        };
        let replace = match &best {
            None => true,
            Some(existing) => {
                candidate.objective > existing.objective + 1e-12
                    || ((candidate.objective - existing.objective).abs() <= 1e-12
                        && (candidate.cut_hard, side_bp[0].min(side_bp[1]))
                            > (existing.cut_hard, cut_min_bp(&existing.labels, contigs)))
            }
        };
        if replace {
            best = Some(candidate);
        }
    }
    best
}

fn greedy_seed_cut(
    anchors: &[usize],
    edges: &[SignedEdge],
    contigs: &[Contig],
) -> Option<SeedCut> {
    let seed = edges
        .iter()
        .filter(|edge| edge.hard)
        .max_by(|a, b| a.weight.partial_cmp(&b.weight).unwrap_or(Ordering::Equal))?;
    let mut labels = HashMap::from([(seed.left, 0u8), (seed.right, 1u8)]);
    let mut unassigned: HashSet<usize> = anchors.iter().copied().collect();
    unassigned.remove(&seed.left);
    unassigned.remove(&seed.right);
    while !unassigned.is_empty() {
        let mut best_member = None;
        let mut best_evidence = f64::NEG_INFINITY;
        let mut best_scores = [0.0f64; 2];
        for &member in &unassigned {
            let scores = signed_assignment_scores(member, &labels, edges);
            let evidence = scores[0] + scores[1];
            if evidence > best_evidence + 1e-12 {
                best_member = Some(member);
                best_evidence = evidence;
                best_scores = scores;
            }
        }
        let member = best_member?;
        let side = if (best_scores[0] - best_scores[1]).abs() > 1e-12 {
            usize::from(best_scores[1] > best_scores[0])
        } else {
            let bp0 = labels
                .iter()
                .filter(|(_, side)| **side == 0)
                .map(|(&index, _)| contigs[index].seq.len())
                .sum::<usize>();
            let bp1 = labels
                .iter()
                .filter(|(_, side)| **side == 1)
                .map(|(&index, _)| contigs[index].seq.len())
                .sum::<usize>();
            usize::from(bp0 > bp1)
        } as u8;
        labels.insert(member, side);
        unassigned.remove(&member);
    }
    let (objective, cut_hard) = score_cut(&labels, edges);
    Some(SeedCut {
        labels,
        cut_hard,
        total_hard: edges.iter().filter(|edge| edge.hard).count(),
        objective,
    })
}

fn signed_assignment_scores(
    member: usize,
    labels: &HashMap<usize, u8>,
    edges: &[SignedEdge],
) -> [f64; 2] {
    let mut scores = [0.0f64; 2];
    for edge in edges {
        let other = if edge.left == member {
            edge.right
        } else if edge.right == member {
            edge.left
        } else {
            continue;
        };
        let Some(&side) = labels.get(&other) else {
            continue;
        };
        let target = if edge.hard { 1 - side } else { side };
        scores[usize::from(target)] += edge.weight;
    }
    scores
}

fn score_cut(labels: &HashMap<usize, u8>, edges: &[SignedEdge]) -> (f64, usize) {
    let mut objective = 0.0;
    let mut cut_hard = 0usize;
    for edge in edges {
        let Some(&left) = labels.get(&edge.left) else {
            continue;
        };
        let Some(&right) = labels.get(&edge.right) else {
            continue;
        };
        if edge.hard {
            if left != right {
                objective += edge.weight;
                cut_hard += 1;
            }
        } else if left == right {
            objective += edge.weight;
        }
    }
    (objective, cut_hard)
}

fn cut_min_bp(labels: &HashMap<usize, u8>, contigs: &[Contig]) -> usize {
    let mut bp = [0usize; 2];
    for (&member, &side) in labels {
        bp[usize::from(side)] += contigs[member].seq.len();
    }
    bp[0].min(bp[1])
}

fn propagation_scores(
    member: usize,
    biological_anchors: &HashSet<usize>,
    labels: &HashMap<usize, u8>,
    cheap: &[CheapFeature],
    top_k: usize,
) -> [f64; 2] {
    let mut scores = [Vec::new(), Vec::new()];
    for &anchor in biological_anchors {
        let Some(&side) = labels.get(&anchor) else {
            continue;
        };
        scores[usize::from(side)].push(cheap_attraction(&cheap[member], &cheap[anchor]));
    }
    let mut out = [0.0f64; 2];
    for side in 0..2 {
        scores[side].sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
        let keep = top_k.max(1).min(scores[side].len());
        if keep > 0 {
            out[side] = scores[side][..keep].iter().sum::<f64>() / keep as f64;
        }
    }
    out
}

fn cheap_feature(contig: &Contig, coverage: Option<&CoverageTable>) -> CheapFeature {
    CheapFeature {
        gc: gc_fraction(&contig.seq),
        kmer: canonical_5mer_frequency(&contig.seq),
        coverage: coverage
            .and_then(|table| table.values.get(&contig.id))
            .cloned()
            .unwrap_or_default(),
    }
}

fn cheap_attraction(left: &CheapFeature, right: &CheapFeature) -> f64 {
    let composition = (-hellinger(&left.kmer, &right.kmer) / 0.24).exp();
    let gc = (-(left.gc - right.gc).abs() / 0.055).exp();
    let coverage = coverage_similarity(&left.coverage, &right.coverage);
    let mut weighted = 0.42 * composition + 0.05 * gc;
    let mut weight = 0.47;
    if let Some(value) = coverage {
        let coverage_weight = if left.coverage.len() >= 3 { 0.48 } else { 0.34 };
        weighted += coverage_weight * value;
        weight += coverage_weight;
    }
    (weighted / weight.max(1e-12)).clamp(0.0, 1.0)
}

fn coverage_similarity(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.is_empty() || left.len() != right.len() {
        return None;
    }
    let log_left: Vec<f64> = left.iter().map(|value| (value + 0.5).ln()).collect();
    let log_right: Vec<f64> = right.iter().map(|value| (value + 0.5).ln()).collect();
    let mean_abs_ratio = log_left
        .iter()
        .zip(log_right.iter())
        .map(|(a, b)| (a - b).abs())
        .sum::<f64>()
        / left.len() as f64;
    let ratio_score = (-mean_abs_ratio / 0.75).exp();
    if left.len() < 3 {
        return Some(ratio_score);
    }
    let correlation = pearson(&log_left, &log_right);
    Some(
        correlation
            .map(|value| 0.68 * ratio_score + 0.32 * ((value + 1.0) * 0.5))
            .unwrap_or(ratio_score)
            .clamp(0.0, 1.0),
    )
}

fn pearson(left: &[f64], right: &[f64]) -> Option<f64> {
    let mean_left = left.iter().sum::<f64>() / left.len() as f64;
    let mean_right = right.iter().sum::<f64>() / right.len() as f64;
    let mut numerator = 0.0;
    let mut left_sq = 0.0;
    let mut right_sq = 0.0;
    for (&a, &b) in left.iter().zip(right.iter()) {
        let da = a - mean_left;
        let db = b - mean_right;
        numerator += da * db;
        left_sq += da * da;
        right_sq += db * db;
    }
    let denominator = (left_sq * right_sq).sqrt();
    (denominator > 1e-12).then_some((numerator / denominator).clamp(-1.0, 1.0))
}

fn hellinger(left: &[f64; 1024], right: &[f64; 1024]) -> f64 {
    (0.5
        * left
            .iter()
            .zip(right.iter())
            .map(|(a, b)| (a.sqrt() - b.sqrt()).powi(2))
            .sum::<f64>())
    .sqrt()
}

fn gc_fraction(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut valid = 0usize;
    for &base in seq {
        match base.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc += 1;
                valid += 1;
            }
            b'A' | b'T' => valid += 1,
            _ => {}
        }
    }
    if valid == 0 { 0.0 } else { gc as f64 / valid as f64 }
}

fn canonical_5mer_frequency(seq: &[u8]) -> [f64; 1024] {
    let mut counts = [0.0; 1024];
    let mut total = 0.0;
    for window in seq.windows(5) {
        if let (Some(forward), Some(reverse)) = (encode_5mer(window, false), encode_5mer(window, true)) {
            counts[forward.min(reverse)] += 1.0;
            total += 1.0;
        }
    }
    if total > 0.0 {
        for value in &mut counts {
            *value /= total;
        }
    }
    counts
}

fn encode_5mer(window: &[u8], reverse_complement: bool) -> Option<usize> {
    let mut code = 0usize;
    if reverse_complement {
        for &base in window.iter().rev() {
            code = (code << 2) | base_code(base, true)?;
        }
    } else {
        for &base in window {
            code = (code << 2) | base_code(base, false)?;
        }
    }
    Some(code)
}

fn base_code(base: u8, complement: bool) -> Option<usize> {
    let value = match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    Some(if complement { 3 - value } else { value })
}

fn summarize_bins(contigs: &[Contig], result: &BinningResult) -> Vec<BinSummary> {
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for (contig_index, assignment) in result.assignments.iter().enumerate() {
        if let Some(bin_index) = assignment.bin_index {
            groups.entry(bin_index).or_default().push(contig_index);
        }
    }
    let mut ids: Vec<usize> = groups.keys().copied().collect();
    ids.sort_unstable();
    ids.into_iter()
        .map(|bin_index| {
            let members = &groups[&bin_index];
            let total_bp = members.iter().map(|&index| contigs[index].seq.len()).sum::<usize>();
            let gc_weighted = members
                .iter()
                .map(|&index| gc_fraction(&contigs[index].seq) * contigs[index].seq.len() as f64)
                .sum::<f64>();
            BinSummary {
                bin_index,
                contig_count: members.len(),
                total_bp,
                mean_gc: if total_bp == 0 { 0.0 } else { gc_weighted / total_bp as f64 },
            }
        })
        .collect()
}
