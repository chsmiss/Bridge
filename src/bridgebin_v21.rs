use crate::bridgebin::{Assignment, BinSummary, BinningResult, Contig};
use crate::bridgebin_reconcile::MarkerTable;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct BridgeBinV21Config {
    pub min_pair_confidence: f64,
    pub split_max_same: f64,
    pub join_min_same: f64,
    pub rescue_min_same: f64,
    pub rescue_margin: f64,
    pub min_pair_support: usize,
    pub top_pair_support: usize,
    pub min_subbin_bp: usize,
}

impl Default for BridgeBinV21Config {
    fn default() -> Self {
        Self {
            min_pair_confidence: 0.80,
            split_max_same: 0.12,
            join_min_same: 0.88,
            rescue_min_same: 0.84,
            rescue_margin: 0.08,
            min_pair_support: 2,
            top_pair_support: 8,
            min_subbin_bp: 20_000,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct PairScore {
    pub same_probability: f64,
    pub confidence: f64,
    pub source: String,
}

#[derive(Clone, Debug, Default)]
pub struct PairScoreTable {
    pub values: HashMap<(String, String), PairScore>,
}

#[derive(Clone, Debug, Default)]
pub struct BridgeBinV21Stats {
    pub input_bins: usize,
    pub output_bins: usize,
    pub conflicted_input_bins: usize,
    pub split_bins: usize,
    pub hard_negative_pairs: usize,
    pub marker_negative_pairs: usize,
    pub positive_pairs: usize,
    pub rescued_contigs: usize,
    pub ambiguous_residuals: usize,
}

#[derive(Clone, Debug)]
struct RefinedBin {
    members: Vec<usize>,
    bp: usize,
    markers: HashSet<String>,
}

pub fn read_pair_score_table<P: AsRef<Path>>(path: P) -> io::Result<PairScoreTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| invalid("empty pair score table"))?;
    let columns: Vec<&str> = header.trim().split('\t').collect();
    let find = |names: &[&str]| -> Option<usize> {
        columns
            .iter()
            .position(|column| names.iter().any(|name| column.eq_ignore_ascii_case(name)))
    };
    let left_col = find(&["left", "source", "contig_a", "contig1"])
        .ok_or_else(|| invalid("pair score table needs left/source column"))?;
    let right_col = find(&["right", "target", "contig_b", "contig2"])
        .ok_or_else(|| invalid("pair score table needs right/target column"))?;
    let same_col = find(&[
        "p_same",
        "same_probability",
        "same_genome",
        "probability",
        "score",
    ])
    .ok_or_else(|| invalid("pair score table needs p_same/same_probability column"))?;
    let confidence_col = find(&["confidence", "model_confidence", "pair_confidence"]);
    let source_col = find(&["model", "source_name", "evidence_source", "kind"]);

    let mut values: HashMap<(String, String), PairScore> = HashMap::new();
    for (row_index, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let left = fields.get(left_col).copied().unwrap_or("").trim();
        let right = fields.get(right_col).copied().unwrap_or("").trim();
        if left.is_empty() || right.is_empty() || left == right {
            continue;
        }
        let same_probability = parse_probability(
            fields.get(same_col).copied().unwrap_or(""),
            row_index + 2,
            "same probability",
        )?;
        let confidence = match confidence_col.and_then(|column| fields.get(column).copied()) {
            Some(raw) if !raw.trim().is_empty() && raw.trim() != "." => {
                parse_probability(raw, row_index + 2, "pair confidence")?
            }
            _ => 1.0,
        };
        let source = source_col
            .and_then(|column| fields.get(column).copied())
            .unwrap_or("")
            .trim()
            .to_string();
        let key = ordered_pair(left, right);
        match values.get_mut(&key) {
            Some(existing) => {
                let existing_decisiveness = (existing.same_probability - 0.5).abs();
                let new_decisiveness = (same_probability - 0.5).abs();
                if confidence > existing.confidence
                    || ((confidence - existing.confidence).abs() <= 1e-12
                        && new_decisiveness > existing_decisiveness)
                {
                    *existing = PairScore {
                        same_probability,
                        confidence,
                        source,
                    };
                }
            }
            None => {
                values.insert(
                    key,
                    PairScore {
                        same_probability,
                        confidence,
                        source,
                    },
                );
            }
        }
    }
    Ok(PairScoreTable { values })
}

pub fn refine_bins_v21(
    contigs: &[Contig],
    markers: Option<&MarkerTable>,
    initial: BinningResult,
    pair_scores: &PairScoreTable,
    cfg: &BridgeBinV21Config,
) -> (BinningResult, BridgeBinV21Stats) {
    let by_id: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| (contig.id.as_str(), index))
        .collect();
    let indexed_scores = index_scores(pair_scores, &by_id);
    let marker_sets: Vec<HashSet<String>> = contigs
        .iter()
        .map(|contig| {
            markers
                .and_then(|table| table.values.get(&contig.id))
                .cloned()
                .unwrap_or_default()
        })
        .collect();

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut residuals = Vec::new();
    for (contig_index, assignment) in initial.assignments.iter().enumerate() {
        match assignment.bin_index {
            Some(bin_index) => groups.entry(bin_index).or_default().push(contig_index),
            None => residuals.push(contig_index),
        }
    }

    let mut stats = BridgeBinV21Stats {
        input_bins: groups.len(),
        ..Default::default()
    };
    let mut refined = Vec::new();
    let mut group_ids: Vec<usize> = groups.keys().copied().collect();
    group_ids.sort_unstable();

    for group_id in group_ids {
        let members = groups.remove(&group_id).unwrap_or_default();
        let hard = hard_conflicts(&members, &indexed_scores, &marker_sets, cfg, &mut stats);
        if hard.is_empty() {
            refined.push(make_bin(members, contigs, &marker_sets));
            continue;
        }
        stats.conflicted_input_bins += 1;
        let (mut pieces, mut dropped) = split_group(
            &members,
            &hard,
            &indexed_scores,
            contigs,
            &marker_sets,
            cfg,
            &mut stats,
        );
        if pieces.len() > 1 {
            stats.split_bins += 1;
        }
        refined.append(&mut pieces);
        residuals.append(&mut dropped);
    }

    residuals.sort_unstable();
    residuals.dedup();
    for contig_index in residuals {
        let mut candidates = Vec::new();
        for (bin_index, bin) in refined.iter().enumerate() {
            if marker_conflict_with_bin(contig_index, bin, &marker_sets) {
                continue;
            }
            if let Some((score, support)) =
                aggregate_to_bin(contig_index, bin, &indexed_scores, cfg)
            {
                if support >= cfg.min_pair_support {
                    candidates.push((bin_index, score, support));
                }
            }
        }
        candidates.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(Ordering::Equal)
                .then_with(|| b.2.cmp(&a.2))
                .then_with(|| a.0.cmp(&b.0))
        });
        let best = candidates.first().copied();
        let second_score = candidates.get(1).map(|candidate| candidate.1).unwrap_or(0.0);
        match best {
            Some((bin_index, score, _))
                if score >= cfg.rescue_min_same && score - second_score >= cfg.rescue_margin =>
            {
                add_to_bin(&mut refined[bin_index], contig_index, contigs, &marker_sets);
                stats.rescued_contigs += 1;
            }
            _ => {
                stats.ambiguous_residuals += 1;
            }
        }
    }

    refined.sort_by_key(|bin| std::cmp::Reverse(bin.bp));
    stats.output_bins = refined.len();
    let mut assigned_bin = vec![None; contigs.len()];
    let mut assigned_score = vec![0.0; contigs.len()];
    for (bin_index, bin) in refined.iter().enumerate() {
        for &contig_index in &bin.members {
            assigned_bin[contig_index] = Some(bin_index);
            assigned_score[contig_index] = 1.0;
        }
    }

    let assignments = contigs
        .iter()
        .enumerate()
        .map(|(contig_index, contig)| Assignment {
            contig_id: contig.id.clone(),
            bin_index: assigned_bin[contig_index],
            score: assigned_score[contig_index],
            length: contig.seq.len(),
        })
        .collect();
    let bins = refined
        .iter()
        .enumerate()
        .map(|(bin_index, bin)| BinSummary {
            bin_index,
            contig_count: bin.members.len(),
            total_bp: bin.bp,
            mean_gc: weighted_gc(bin, contigs),
        })
        .collect();

    (BinningResult { assignments, bins }, stats)
}

fn hard_conflicts(
    members: &[usize],
    scores: &HashMap<(usize, usize), PairScore>,
    marker_sets: &[HashSet<String>],
    cfg: &BridgeBinV21Config,
    stats: &mut BridgeBinV21Stats,
) -> HashSet<(usize, usize)> {
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let mut hard = HashSet::new();
    for (&(left, right), score) in scores {
        if member_set.contains(&left)
            && member_set.contains(&right)
            && score.confidence >= cfg.min_pair_confidence
            && score.same_probability <= cfg.split_max_same
        {
            hard.insert((left.min(right), left.max(right)));
        }
    }
    stats.hard_negative_pairs += hard.len();

    let mut marker_to_members: HashMap<&str, Vec<usize>> = HashMap::new();
    for &member in members {
        for marker in &marker_sets[member] {
            marker_to_members
                .entry(marker.as_str())
                .or_default()
                .push(member);
        }
    }
    for marker_members in marker_to_members.values() {
        if marker_members.len() < 2 {
            continue;
        }
        for left_pos in 0..marker_members.len() {
            for right_pos in (left_pos + 1)..marker_members.len() {
                let left = marker_members[left_pos];
                let right = marker_members[right_pos];
                if hard.insert((left.min(right), left.max(right))) {
                    stats.marker_negative_pairs += 1;
                }
            }
        }
    }
    hard
}

#[allow(clippy::too_many_arguments)]
fn split_group(
    members: &[usize],
    hard: &HashSet<(usize, usize)>,
    scores: &HashMap<(usize, usize), PairScore>,
    contigs: &[Contig],
    marker_sets: &[HashSet<String>],
    cfg: &BridgeBinV21Config,
    stats: &mut BridgeBinV21Stats,
) -> (Vec<RefinedBin>, Vec<usize>) {
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let mut clusters: Vec<Option<Vec<usize>>> = members
        .iter()
        .map(|member| Some(vec![*member]))
        .collect();
    let mut cluster_of: HashMap<usize, usize> = members
        .iter()
        .enumerate()
        .map(|(cluster, member)| (*member, cluster))
        .collect();

    let mut positive_edges: Vec<(f64, f64, usize, usize)> = scores
        .iter()
        .filter_map(|(&(left, right), score)| {
            (member_set.contains(&left)
                && member_set.contains(&right)
                && score.confidence >= cfg.min_pair_confidence
                && score.same_probability >= cfg.join_min_same)
                .then_some((score.same_probability, score.confidence, left, right))
        })
        .collect();
    positive_edges.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(Ordering::Equal)
            .then_with(|| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal))
    });
    stats.positive_pairs += positive_edges.len();

    for (_, _, left, right) in positive_edges {
        let left_cluster = cluster_of[&left];
        let right_cluster = cluster_of[&right];
        if left_cluster == right_cluster {
            continue;
        }
        let left_members = clusters[left_cluster].as_ref().unwrap();
        let right_members = clusters[right_cluster].as_ref().unwrap();
        if clusters_conflict(left_members, right_members, hard) {
            continue;
        }
        let (keep, remove) = if left_members.len() >= right_members.len() {
            (left_cluster, right_cluster)
        } else {
            (right_cluster, left_cluster)
        };
        let removed = clusters[remove].take().unwrap();
        for member in &removed {
            cluster_of.insert(*member, keep);
        }
        clusters[keep].as_mut().unwrap().extend(removed);
    }

    let mut bins = Vec::new();
    let mut dropped = Vec::new();
    for cluster in clusters.into_iter().flatten() {
        let bp = cluster
            .iter()
            .map(|index| contigs[*index].seq.len())
            .sum::<usize>();
        if bp >= cfg.min_subbin_bp {
            bins.push(make_bin(cluster, contigs, marker_sets));
        } else {
            dropped.extend(cluster);
        }
    }
    bins.sort_by_key(|bin| std::cmp::Reverse(bin.bp));
    (bins, dropped)
}

fn clusters_conflict(left: &[usize], right: &[usize], hard: &HashSet<(usize, usize)>) -> bool {
    left.iter().any(|a| {
        right.iter().any(|b| {
            let key = ((*a).min(*b), (*a).max(*b));
            hard.contains(&key)
        })
    })
}

fn aggregate_to_bin(
    contig_index: usize,
    bin: &RefinedBin,
    scores: &HashMap<(usize, usize), PairScore>,
    cfg: &BridgeBinV21Config,
) -> Option<(f64, usize)> {
    let mut support = Vec::new();
    for &member in &bin.members {
        if member == contig_index {
            continue;
        }
        let key = (contig_index.min(member), contig_index.max(member));
        let Some(score) = scores.get(&key) else {
            continue;
        };
        if score.confidence < cfg.min_pair_confidence {
            continue;
        }
        if score.same_probability <= cfg.split_max_same {
            return None;
        }
        support.push((score.same_probability, score.confidence));
    }
    support.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
    support.truncate(cfg.top_pair_support.max(1));
    if support.is_empty() {
        return None;
    }
    let total_weight = support
        .iter()
        .map(|(_, confidence)| confidence)
        .sum::<f64>();
    if total_weight <= 0.0 {
        return None;
    }
    let weighted = support
        .iter()
        .map(|(same, confidence)| same * confidence)
        .sum::<f64>()
        / total_weight;
    Some((weighted.clamp(0.0, 1.0), support.len()))
}

fn marker_conflict_with_bin(
    contig_index: usize,
    bin: &RefinedBin,
    marker_sets: &[HashSet<String>],
) -> bool {
    !marker_sets[contig_index].is_disjoint(&bin.markers)
}

fn make_bin(
    members: Vec<usize>,
    contigs: &[Contig],
    marker_sets: &[HashSet<String>],
) -> RefinedBin {
    let bp = members
        .iter()
        .map(|index| contigs[*index].seq.len())
        .sum();
    let mut markers = HashSet::new();
    for &member in &members {
        markers.extend(marker_sets[member].iter().cloned());
    }
    RefinedBin {
        members,
        bp,
        markers,
    }
}

fn add_to_bin(
    bin: &mut RefinedBin,
    contig_index: usize,
    contigs: &[Contig],
    marker_sets: &[HashSet<String>],
) {
    bin.bp += contigs[contig_index].seq.len();
    bin.members.push(contig_index);
    bin.markers
        .extend(marker_sets[contig_index].iter().cloned());
}

fn weighted_gc(bin: &RefinedBin, contigs: &[Contig]) -> f64 {
    let mut gc = 0usize;
    let mut valid = 0usize;
    for &index in &bin.members {
        for &base in &contigs[index].seq {
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

fn index_scores(
    table: &PairScoreTable,
    by_id: &HashMap<&str, usize>,
) -> HashMap<(usize, usize), PairScore> {
    table
        .values
        .iter()
        .filter_map(|((left, right), score)| {
            let (&a, &b) = (by_id.get(left.as_str())?, by_id.get(right.as_str())?);
            Some(((a.min(b), a.max(b)), score.clone()))
        })
        .collect()
}

fn parse_probability(raw: &str, row: usize, label: &str) -> io::Result<f64> {
    let value = raw
        .trim()
        .parse::<f64>()
        .map_err(|_| invalid(format!("invalid {label} at row {row}: {raw:?}")))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(invalid(format!("{label} out of range at row {row}")));
    }
    Ok(value)
}

fn ordered_pair(left: &str, right: &str) -> (String, String) {
    if left <= right {
        (left.to_string(), right.to_string())
    } else {
        (right.to_string(), left.to_string())
    }
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(id: &str) -> Contig {
        Contig {
            id: id.to_string(),
            seq: b"ACGT".repeat(100),
        }
    }

    fn initial_one_bin(contigs: &[Contig]) -> BinningResult {
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
                mean_gc: 0.5,
            }],
        }
    }

    #[test]
    fn calibrated_pair_scores_split_coverage_identical_genomes() {
        let contigs = vec![contig("a1"), contig("a2"), contig("b1"), contig("b2")];
        let mut scores = PairScoreTable::default();
        for (left, right, same) in [
            ("a1", "a2", 0.99),
            ("b1", "b2", 0.99),
            ("a1", "b1", 0.01),
            ("a1", "b2", 0.01),
            ("a2", "b1", 0.01),
            ("a2", "b2", 0.01),
        ] {
            scores.values.insert(
                ordered_pair(left, right),
                PairScore {
                    same_probability: same,
                    confidence: 0.99,
                    source: "test".to_string(),
                },
            );
        }
        let cfg = BridgeBinV21Config {
            min_subbin_bp: 1,
            min_pair_support: 1,
            ..Default::default()
        };
        let (result, stats) = refine_bins_v21(
            &contigs,
            None,
            initial_one_bin(&contigs),
            &scores,
            &cfg,
        );
        assert_eq!(stats.split_bins, 1);
        let bins: HashMap<&str, Option<usize>> = result
            .assignments
            .iter()
            .map(|assignment| (assignment.contig_id.as_str(), assignment.bin_index))
            .collect();
        assert_eq!(bins["a1"], bins["a2"]);
        assert_eq!(bins["b1"], bins["b2"]);
        assert_ne!(bins["a1"], bins["b1"]);
    }

    #[test]
    fn duplicated_single_copy_marker_creates_split_pressure() {
        let contigs = vec![contig("a"), contig("b")];
        let markers = MarkerTable {
            values: HashMap::from([
                ("a".to_string(), HashSet::from(["SCG1".to_string()])),
                ("b".to_string(), HashSet::from(["SCG1".to_string()])),
            ]),
        };
        let cfg = BridgeBinV21Config {
            min_subbin_bp: 1,
            min_pair_support: 1,
            ..Default::default()
        };
        let (result, stats) = refine_bins_v21(
            &contigs,
            Some(&markers),
            initial_one_bin(&contigs),
            &PairScoreTable::default(),
            &cfg,
        );
        assert_eq!(stats.marker_negative_pairs, 1);
        assert_ne!(result.assignments[0].bin_index, result.assignments[1].bin_index);
    }
}
