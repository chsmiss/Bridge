use crate::bridgebin::{
    Assignment, BinSummary, BinningResult, BridgeBinConfig, Contig, CoverageTable,
};
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, Default)]
pub struct MarkerTable {
    pub values: HashMap<String, HashSet<String>>,
}

#[derive(Clone, Debug)]
pub struct ReconcileConfig {
    pub enabled: bool,
    pub merge_threshold: f64,
    pub merge_margin: f64,
    pub min_composition_score: f64,
    pub min_coverage_score: f64,
    pub same_coverage_min_composition: f64,
    pub max_gc_delta: f64,
    pub post_rescue_threshold: f64,
    pub post_rescue_margin: f64,
    pub max_merges: usize,
}

impl Default for ReconcileConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            merge_threshold: 0.72,
            merge_margin: 0.015,
            min_composition_score: 0.58,
            min_coverage_score: 0.62,
            same_coverage_min_composition: 0.82,
            max_gc_delta: 0.055,
            post_rescue_threshold: 0.72,
            post_rescue_margin: 0.03,
            max_merges: 256,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct ReconcileStats {
    pub initial_bins: usize,
    pub final_bins: usize,
    pub merges: usize,
    pub rescued_contigs: usize,
    pub marker_blocked_pairs: usize,
}

#[derive(Clone, Debug)]
struct BinNode {
    members: Vec<usize>,
    bp: usize,
    gc_sum: f64,
    kmer_counts: [f64; 1024],
    kmer_total: f64,
    coverage_sum: Vec<f64>,
    markers: HashSet<String>,
}

#[derive(Clone, Copy, Debug)]
struct Similarity {
    combined: f64,
    composition: f64,
    coverage: Option<f64>,
}

pub fn read_marker_table<P: AsRef<Path>>(path: P) -> io::Result<MarkerTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut values: HashMap<String, HashSet<String>> = HashMap::new();
    let mut first_data_row = true;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("marker row {} needs contig and marker columns", line_no + 1),
            ));
        }
        if first_data_row
            && fields[0].eq_ignore_ascii_case("contig")
            && fields[1].to_ascii_lowercase().starts_with("marker")
        {
            first_data_row = false;
            continue;
        }
        first_data_row = false;
        for marker in fields[1].split(',').filter(|m| !m.is_empty()) {
            values
                .entry(fields[0].to_string())
                .or_default()
                .insert(marker.to_string());
        }
    }
    Ok(MarkerTable { values })
}

pub fn reconcile_bins(
    contigs: &[Contig],
    coverage: Option<&CoverageTable>,
    markers: Option<&MarkerTable>,
    initial: BinningResult,
    base_cfg: &BridgeBinConfig,
    cfg: &ReconcileConfig,
) -> (BinningResult, ReconcileStats) {
    if !cfg.enabled || initial.bins.len() <= 1 {
        let bins = initial.bins.len();
        return (
            initial,
            ReconcileStats {
                initial_bins: bins,
                final_bins: bins,
                ..Default::default()
            },
        );
    }

    let contig_index: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| (contig.id.as_str(), index))
        .collect();
    let mut seed_members = vec![Vec::new(); initial.bins.len()];
    for assignment in &initial.assignments {
        if let (Some(bin), Some(&index)) = (
            assignment.bin_index,
            contig_index.get(assignment.contig_id.as_str()),
        ) {
            if let Some(members) = seed_members.get_mut(bin) {
                members.push(index);
            }
        }
    }

    let mut nodes: Vec<BinNode> = seed_members
        .into_iter()
        .filter(|members| !members.is_empty())
        .map(|members| BinNode::from_members(members, contigs, coverage, markers))
        .collect();
    let mut stats = ReconcileStats {
        initial_bins: nodes.len(),
        ..Default::default()
    };

    while nodes.len() > 1 && stats.merges < cfg.max_merges {
        let choices: Vec<Option<(usize, f64, f64)>> = (0..nodes.len())
            .map(|index| best_two(index, &nodes, base_cfg, cfg, &mut stats))
            .collect();
        let mut selected: Option<(usize, usize, f64)> = None;

        for left in 0..nodes.len() {
            let Some((right, score, second)) = choices[left] else {
                continue;
            };
            if right <= left
                || score < cfg.merge_threshold
                || score - second < cfg.merge_margin
            {
                continue;
            }
            let Some((back, back_score, back_second)) = choices[right] else {
                continue;
            };
            if back != left || back_score - back_second < cfg.merge_margin {
                continue;
            }
            let reciprocal_score = score.min(back_score);
            if selected
                .map(|(_, _, previous)| reciprocal_score > previous)
                .unwrap_or(true)
            {
                selected = Some((left, right, reciprocal_score));
            }
        }

        let Some((left, right, _)) = selected else {
            break;
        };
        let other = nodes.remove(right);
        nodes[left].absorb(other);
        stats.merges += 1;
    }

    nodes.sort_by_key(|node| Reverse(node.bp));
    let mut rescued_scores = HashMap::new();
    let mut unbinned: Vec<usize> = initial
        .assignments
        .iter()
        .filter(|assignment| assignment.bin_index.is_none())
        .filter_map(|assignment| contig_index.get(assignment.contig_id.as_str()).copied())
        .filter(|&index| contigs[index].seq.len() >= base_cfg.min_contig_len)
        .collect();
    unbinned.sort_by_key(|&index| Reverse(contigs[index].seq.len()));

    for index in unbinned {
        let probe = BinNode::from_members(vec![index], contigs, coverage, markers);
        let mut best: Option<(usize, f64)> = None;
        let mut second = 0.0;
        for (bin, node) in nodes.iter().enumerate() {
            if marker_conflict(&probe, node) {
                continue;
            }
            let similarity = node_similarity(&probe, node, base_cfg);
            if !compatible(&probe, node, similarity, cfg) {
                continue;
            }
            if best
                .map(|(_, score)| similarity.combined > score)
                .unwrap_or(true)
            {
                second = best.map(|(_, score)| score).unwrap_or(0.0);
                best = Some((bin, similarity.combined));
            } else if similarity.combined > second {
                second = similarity.combined;
            }
        }
        if let Some((bin, score)) = best {
            if score >= cfg.post_rescue_threshold && score - second >= cfg.post_rescue_margin {
                nodes[bin].add_contig(index, &contigs[index], coverage, markers);
                rescued_scores.insert(index, score);
                stats.rescued_contigs += 1;
            }
        }
    }

    nodes.sort_by_key(|node| Reverse(node.bp));
    let mut contig_to_bin = vec![None; contigs.len()];
    for (bin, node) in nodes.iter().enumerate() {
        for &index in &node.members {
            contig_to_bin[index] = Some(bin);
        }
    }

    let initial_scores: HashMap<&str, f64> = initial
        .assignments
        .iter()
        .map(|assignment| (assignment.contig_id.as_str(), assignment.score))
        .collect();
    let assignments = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| Assignment {
            contig_id: contig.id.clone(),
            bin_index: contig_to_bin[index],
            score: rescued_scores
                .get(&index)
                .copied()
                .or_else(|| initial_scores.get(contig.id.as_str()).copied())
                .unwrap_or(0.0),
            length: contig.seq.len(),
        })
        .collect();
    let bins = nodes
        .iter()
        .enumerate()
        .map(|(bin_index, node)| BinSummary {
            bin_index,
            contig_count: node.members.len(),
            total_bp: node.bp,
            mean_gc: node.mean_gc(),
        })
        .collect();

    stats.final_bins = nodes.len();
    (BinningResult { assignments, bins }, stats)
}

fn best_two(
    index: usize,
    nodes: &[BinNode],
    base_cfg: &BridgeBinConfig,
    cfg: &ReconcileConfig,
    stats: &mut ReconcileStats,
) -> Option<(usize, f64, f64)> {
    let mut best: Option<(usize, f64)> = None;
    let mut second = 0.0;
    for other in 0..nodes.len() {
        if other == index {
            continue;
        }
        if marker_conflict(&nodes[index], &nodes[other]) {
            stats.marker_blocked_pairs += 1;
            continue;
        }
        let similarity = node_similarity(&nodes[index], &nodes[other], base_cfg);
        if !compatible(&nodes[index], &nodes[other], similarity, cfg) {
            continue;
        }
        if best
            .map(|(_, score)| similarity.combined > score)
            .unwrap_or(true)
        {
            second = best.map(|(_, score)| score).unwrap_or(0.0);
            best = Some((other, similarity.combined));
        } else if similarity.combined > second {
            second = similarity.combined;
        }
    }
    best.map(|(other, score)| (other, score, second))
}

fn compatible(left: &BinNode, right: &BinNode, score: Similarity, cfg: &ReconcileConfig) -> bool {
    if (left.mean_gc() - right.mean_gc()).abs() > cfg.max_gc_delta
        || score.composition < cfg.min_composition_score
    {
        return false;
    }
    if let Some(coverage) = score.coverage {
        if coverage < cfg.min_coverage_score {
            return false;
        }
        if coverage >= 0.97 && score.composition < cfg.same_coverage_min_composition {
            return false;
        }
    }
    true
}

fn marker_conflict(left: &BinNode, right: &BinNode) -> bool {
    if left.markers.len() <= right.markers.len() {
        left.markers
            .iter()
            .any(|marker| right.markers.contains(marker))
    } else {
        right
            .markers
            .iter()
            .any(|marker| left.markers.contains(marker))
    }
}

impl BinNode {
    fn from_members(
        members: Vec<usize>,
        contigs: &[Contig],
        coverage: Option<&CoverageTable>,
        markers: Option<&MarkerTable>,
    ) -> Self {
        let mut node = Self {
            members: Vec::new(),
            bp: 0,
            gc_sum: 0.0,
            kmer_counts: [0.0; 1024],
            kmer_total: 0.0,
            coverage_sum: Vec::new(),
            markers: HashSet::new(),
        };
        for index in members {
            node.add_contig(index, &contigs[index], coverage, markers);
        }
        node
    }

    fn add_contig(
        &mut self,
        index: usize,
        contig: &Contig,
        coverage: Option<&CoverageTable>,
        markers: Option<&MarkerTable>,
    ) {
        let length = contig.seq.len();
        let weight = length as f64;
        self.members.push(index);
        self.bp += length;
        self.gc_sum += gc_fraction(&contig.seq) * weight;

        let (counts, total) = canonical_5mer_counts(&contig.seq);
        for (target, value) in self.kmer_counts.iter_mut().zip(counts.iter()) {
            *target += *value;
        }
        self.kmer_total += total;

        if let Some(row) = coverage.and_then(|table| table.values.get(&contig.id)) {
            if self.coverage_sum.is_empty() {
                self.coverage_sum.resize(row.len(), 0.0);
            }
            if self.coverage_sum.len() == row.len() {
                for (target, depth) in self.coverage_sum.iter_mut().zip(row.iter()) {
                    *target += *depth * weight;
                }
            }
        }
        if let Some(hits) = markers.and_then(|table| table.values.get(&contig.id)) {
            self.markers.extend(hits.iter().cloned());
        }
    }

    fn absorb(&mut self, other: Self) {
        self.members.extend(other.members);
        self.bp += other.bp;
        self.gc_sum += other.gc_sum;
        self.kmer_total += other.kmer_total;
        for (target, value) in self.kmer_counts.iter_mut().zip(other.kmer_counts.iter()) {
            *target += *value;
        }
        if self.coverage_sum.is_empty() && !other.coverage_sum.is_empty() {
            self.coverage_sum.resize(other.coverage_sum.len(), 0.0);
        }
        if self.coverage_sum.len() == other.coverage_sum.len() {
            for (target, value) in self.coverage_sum.iter_mut().zip(other.coverage_sum.iter()) {
                *target += *value;
            }
        }
        self.markers.extend(other.markers);
    }

    fn mean_gc(&self) -> f64 {
        if self.bp == 0 {
            0.0
        } else {
            self.gc_sum / self.bp as f64
        }
    }
}

fn node_similarity(left: &BinNode, right: &BinNode, cfg: &BridgeBinConfig) -> Similarity {
    let composition = (-hellinger(left, right) / 0.24).exp();
    let gc = (-(left.mean_gc() - right.mean_gc()).abs() / 0.065).exp();
    let mut weighted = cfg.composition_weight * composition + cfg.gc_weight * gc;
    let mut total_weight = cfg.composition_weight + cfg.gc_weight;

    let coverage = if !left.coverage_sum.is_empty()
        && left.coverage_sum.len() == right.coverage_sum.len()
    {
        let left_bp = left.bp.max(1) as f64;
        let right_bp = right.bp.max(1) as f64;
        let distance = left
            .coverage_sum
            .iter()
            .zip(right.coverage_sum.iter())
            .map(|(a, b)| (((a / left_bp) + 0.5) / ((b / right_bp) + 0.5)).ln().abs())
            .sum::<f64>()
            / left.coverage_sum.len() as f64;
        let similarity = (-distance / 0.80).exp();
        weighted += cfg.coverage_weight * similarity;
        total_weight += cfg.coverage_weight;
        Some(similarity)
    } else {
        None
    };

    Similarity {
        combined: if total_weight <= f64::EPSILON {
            0.0
        } else {
            weighted / total_weight
        },
        composition,
        coverage,
    }
}

fn hellinger(left: &BinNode, right: &BinNode) -> f64 {
    if left.kmer_total <= 0.0 || right.kmer_total <= 0.0 {
        return 1.0;
    }
    let distance = left
        .kmer_counts
        .iter()
        .zip(right.kmer_counts.iter())
        .map(|(a, b)| {
            let pa = *a / left.kmer_total;
            let pb = *b / right.kmer_total;
            (pa.sqrt() - pb.sqrt()).powi(2)
        })
        .sum::<f64>();
    (0.5 * distance).sqrt()
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
    if valid == 0 {
        0.0
    } else {
        gc as f64 / valid as f64
    }
}

fn canonical_5mer_counts(seq: &[u8]) -> ([f64; 1024], f64) {
    let mut counts = [0.0; 1024];
    let mut total = 0.0;
    for window in seq.windows(5) {
        if let (Some(forward), Some(reverse)) = (encode_5mer(window, false), encode_5mer(window, true)) {
            counts[forward.min(reverse)] += 1.0;
            total += 1.0;
        }
    }
    (counts, total)
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
    let code = match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    Some(if complement { 3 - code } else { code })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(id: &str, pattern: &str, repeats: usize) -> Contig {
        Contig {
            id: id.to_string(),
            seq: pattern.repeat(repeats).into_bytes(),
        }
    }

    fn singleton_result(contigs: &[Contig]) -> BinningResult {
        BinningResult {
            assignments: contigs
                .iter()
                .enumerate()
                .map(|(index, contig)| Assignment {
                    contig_id: contig.id.clone(),
                    bin_index: Some(index),
                    score: 1.0,
                    length: contig.seq.len(),
                })
                .collect(),
            bins: contigs
                .iter()
                .enumerate()
                .map(|(index, contig)| BinSummary {
                    bin_index: index,
                    contig_count: 1,
                    total_bp: contig.seq.len(),
                    mean_gc: gc_fraction(&contig.seq),
                })
                .collect(),
        }
    }

    #[test]
    fn reciprocal_reconciliation_merges_split_signatures() {
        let contigs = vec![
            contig("a1", "AAAACAAAAGAAAATAAAC", 220),
            contig("a2", "AAAACAAAAGAAAATAAAT", 220),
            contig("b1", "GGGCGGCCGCGGCGCCGC", 220),
            contig("b2", "GGGCGGCCGCGGCGCCGA", 220),
        ];
        let coverage = CoverageTable {
            sample_names: vec!["s1".to_string(), "s2".to_string()],
            values: HashMap::from([
                ("a1".to_string(), vec![20.0, 10.0]),
                ("a2".to_string(), vec![20.5, 10.2]),
                ("b1".to_string(), vec![20.0, 10.0]),
                ("b2".to_string(), vec![19.8, 10.1]),
            ]),
        };
        let config = ReconcileConfig {
            merge_threshold: 0.50,
            min_composition_score: 0.40,
            same_coverage_min_composition: 0.65,
            merge_margin: 0.0,
            ..Default::default()
        };
        let (result, stats) = reconcile_bins(
            &contigs,
            Some(&coverage),
            None,
            singleton_result(&contigs),
            &BridgeBinConfig::default(),
            &config,
        );
        assert_eq!(stats.final_bins, 2);
        let a1 = result
            .assignments
            .iter()
            .find(|assignment| assignment.contig_id == "a1")
            .unwrap()
            .bin_index;
        let a2 = result
            .assignments
            .iter()
            .find(|assignment| assignment.contig_id == "a2")
            .unwrap()
            .bin_index;
        let b1 = result
            .assignments
            .iter()
            .find(|assignment| assignment.contig_id == "b1")
            .unwrap()
            .bin_index;
        assert_eq!(a1, a2);
        assert_ne!(a1, b1);
    }

    #[test]
    fn shared_single_copy_marker_blocks_merge() {
        let contigs = vec![
            contig("x1", "ACGTACGTAAAACCCCGGGG", 220),
            contig("x2", "ACGTACGTAAAACCCCGGGA", 220),
        ];
        let markers = MarkerTable {
            values: HashMap::from([
                (
                    "x1".to_string(),
                    HashSet::from(["SCG001".to_string()]),
                ),
                (
                    "x2".to_string(),
                    HashSet::from(["SCG001".to_string()]),
                ),
            ]),
        };
        let config = ReconcileConfig {
            merge_threshold: 0.1,
            min_composition_score: 0.0,
            min_coverage_score: 0.0,
            merge_margin: 0.0,
            ..Default::default()
        };
        let (_, stats) = reconcile_bins(
            &contigs,
            None,
            Some(&markers),
            singleton_result(&contigs),
            &BridgeBinConfig::default(),
            &config,
        );
        assert_eq!(stats.final_bins, 2);
        assert!(stats.marker_blocked_pairs > 0);
    }
}
