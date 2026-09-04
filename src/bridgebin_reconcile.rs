use crate::bridgebin::{Assignment, BinSummary, BinningResult, BridgeBinConfig, Contig, CoverageTable};
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
    sig_sum: [f64; 1024],
    sig_total: f64,
    cov_sum: Vec<f64>,
    markers: HashSet<String>,
}

#[derive(Clone, Copy, Debug)]
struct SimilarityParts {
    combined: f64,
    composition: f64,
    coverage: Option<f64>,
}

pub fn read_marker_table<P: AsRef<Path>>(path: P) -> io::Result<MarkerTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut values: HashMap<String, HashSet<String>> = HashMap::new();
    let mut saw_data = false;
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("marker row {} needs contig and marker columns", line_no + 1),
            ));
        }
        if !saw_data
            && fields[0].eq_ignore_ascii_case("contig")
            && fields[1].to_ascii_lowercase().starts_with("marker")
        {
            saw_data = true;
            continue;
        }
        saw_data = true;
        for marker in fields[1].split(',').filter(|x| !x.is_empty()) {
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
        let n = initial.bins.len();
        return (
            initial,
            ReconcileStats {
                initial_bins: n,
                final_bins: n,
                ..Default::default()
            },
        );
    }

    let by_id: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.id.as_str(), i))
        .collect();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); initial.bins.len()];
    for a in &initial.assignments {
        if let (Some(bin), Some(&idx)) = (a.bin_index, by_id.get(a.contig_id.as_str())) {
            if let Some(slot) = members.get_mut(bin) {
                slot.push(idx);
            }
        }
    }

    let mut nodes: Vec<BinNode> = members
        .into_iter()
        .filter(|m| !m.is_empty())
        .map(|m| BinNode::from_members(m, contigs, coverage, markers))
        .collect();
    let mut stats = ReconcileStats {
        initial_bins: nodes.len(),
        ..Default::default()
    };

    while stats.merges < cfg.max_merges && nodes.len() > 1 {
        let mut choices = Vec::with_capacity(nodes.len());
        for i in 0..nodes.len() {
            choices.push(best_two_nodes(i, &nodes, base_cfg, cfg, &mut stats));
        }

        let mut best_pair: Option<(usize, usize, f64)> = None;
        for i in 0..nodes.len() {
            let Some((j, score, second)) = choices[i] else {
                continue;
            };
            if j <= i || score < cfg.merge_threshold || score - second < cfg.merge_margin {
                continue;
            }
            let Some((back, back_score, back_second)) = choices[j] else {
                continue;
            };
            if back != i || back_score - back_second < cfg.merge_margin {
                continue;
            }
            let pair_score = score.min(back_score);
            if best_pair.map(|(_, _, s)| pair_score > s).unwrap_or(true) {
                best_pair = Some((i, j, pair_score));
            }
        }

        let Some((i, j, _)) = best_pair else {
            break;
        };
        let other = nodes.remove(j);
        nodes[i].absorb(other);
        stats.merges += 1;
    }

    nodes.sort_by(|a, b| b.bp.cmp(&a.bp));
    let mut contig_to_bin = vec![None; contigs.len()];
    for (bin, node) in nodes.iter().enumerate() {
        for &idx in &node.members {
            contig_to_bin[idx] = Some(bin);
        }
    }

    let mut unbinned: Vec<usize> = initial
        .assignments
        .iter()
        .filter(|a| a.bin_index.is_none())
        .filter_map(|a| by_id.get(a.contig_id.as_str()).copied())
        .filter(|&idx| contigs[idx].seq.len() >= base_cfg.min_contig_len)
        .collect();
    unbinned.sort_by(|&a, &b| contigs[b].seq.len().cmp(&contigs[a].seq.len()));
    let mut rescued_scores: HashMap<usize, f64> = HashMap::new();
    for idx in unbinned {
        if nodes.is_empty() {
            break;
        }
        let probe = BinNode::from_members(vec![idx], contigs, coverage, markers);
        let mut best: Option<(usize, f64)> = None;
        let mut second = 0.0;
        for (bin, node) in nodes.iter().enumerate() {
            if marker_conflict(&probe, node) {
                continue;
            }
            let parts = node_similarity(&probe, node, base_cfg);
            if !hard_compatible(&probe, node, parts, cfg) {
                continue;
            }
            if best.map(|(_, s)| parts.combined > s).unwrap_or(true) {
                second = best.map(|(_, s)| s).unwrap_or(0.0);
                best = Some((bin, parts.combined));
            } else if parts.combined > second {
                second = parts.combined;
            }
        }
        if let Some((bin, score)) = best {
            if score >= cfg.post_rescue_threshold && score - second >= cfg.post_rescue_margin {
                nodes[bin].add_contig(idx, &contigs[idx], coverage, markers);
                contig_to_bin[idx] = Some(bin);
                rescued_scores.insert(idx, score);
                stats.rescued_contigs += 1;
            }
        }
    }

    nodes.sort_by(|a, b| b.bp.cmp(&a.bp));
    contig_to_bin.fill(None);
    for (bin, node) in nodes.iter().enumerate() {
        for &idx in &node.members {
            contig_to_bin[idx] = Some(bin);
        }
    }

    let initial_score: HashMap<&str, f64> = initial
        .assignments
        .iter()
        .map(|a| (a.contig_id.as_str(), a.score))
        .collect();
    let assignments: Vec<Assignment> = contigs
        .iter()
        .enumerate()
        .map(|(idx, c)| Assignment {
            contig_id: c.id.clone(),
            bin_index: contig_to_bin[idx],
            score: rescued_scores
                .get(&idx)
                .copied()
                .or_else(|| initial_score.get(c.id.as_str()).copied())
                .unwrap_or(0.0),
            length: c.seq.len(),
        })
        .collect();

    let bins: Vec<BinSummary> = nodes
        .iter()
        .enumerate()
        .map(|(bin_index, node)| BinSummary {
            bin_index,
            contig_count: node.members.len(),
            total_bp: node.bp,
            mean_gc: if node.bp == 0 {
                0.0
            } else {
                node.gc_sum / node.bp as f64
            },
        })
        .collect();
    stats.final_bins = bins.len();
    (BinningResult { assignments, bins }, stats)
}

fn best_two_nodes(
    i: usize,
    nodes: &[BinNode],
    base_cfg: &BridgeBinConfig,
    cfg: &ReconcileConfig,
    stats: &mut ReconcileStats,
) -> Option<(usize, f64, f64)> {
    let mut best: Option<(usize, f64)> = None;
    let mut second = 0.0;
    for j in 0..nodes.len() {
        if i == j {
            continue;
        }
        if marker_conflict(&nodes[i], &nodes[j]) {
            stats.marker_blocked_pairs += 1;
            continue;
        }
        let parts = node_similarity(&nodes[i], &nodes[j], base_cfg);
        if !hard_compatible(&nodes[i], &nodes[j], parts, cfg) {
            continue;
        }
        if best.map(|(_, s)| parts.combined > s).unwrap_or(true) {
            second = best.map(|(_, s)| s).unwrap_or(0.0);
            best = Some((j, parts.combined));
        } else if parts.combined > second {
            second = parts.combined;
        }
    }
    best.map(|(j, score)| (j, score, second))
}

fn hard_compatible(a: &BinNode, b: &BinNode, s: SimilarityParts, cfg: &ReconcileConfig) -> bool {
    let gc_a = a.gc_sum / a.bp.max(1) as f64;
    let gc_b = b.gc_sum / b.bp.max(1) as f64;
    if (gc_a - gc_b).abs() > cfg.max_gc_delta || s.composition < cfg.min_composition_score {
        return false;
    }
    if let Some(cov) = s.coverage {
        if cov < cfg.min_coverage_score {
            return false;
        }
        if cov >= 0.97 && s.composition < cfg.same_coverage_min_composition {
            return false;
        }
    }
    true
}

fn marker_conflict(a: &BinNode, b: &BinNode) -> bool {
    if a.markers.len() < b.markers.len() {
        a.markers.iter().any(|m| b.markers.contains(m))
    } else {
        b.markers.iter().any(|m| a.markers.contains(m))
    }
}

impl BinNode {
    fn from_members(
        members: Vec<usize>,
        contigs: &[Contig],
        coverage: Option<&CoverageTable>,
        markers: Option<&MarkerTable>,
    ) -> Self {
        let mut out = Self {
            members: Vec::new(),
            bp: 0,
            gc_sum: 0.0,
            sig_sum: [0.0; 1024],
            sig_total: 0.0,
            cov_sum: Vec::new(),
            markers: HashSet::new(),
        };
        for idx in members {
            out.add_contig(idx, &contigs[idx], coverage, markers);
        }
        out
    }

    fn add_contig(
        &mut self,
        idx: usize,
        contig: &Contig,
        coverage: Option<&CoverageTable>,
        markers: Option<&MarkerTable>,
    ) {
        let len = contig.seq.len();
        let w = len as f64;
        self.members.push(idx);
        self.bp += len;
        self.gc_sum += gc_fraction(&contig.seq) * w;
        let (counts, total) = canonical_5mer_counts(&contig.seq);
        for (dst, src) in self.sig_sum.iter_mut().zip(counts.iter()) {
            *dst += *src;
        }
        self.sig_total += total;
        if let Some(row) = coverage.and_then(|t| t.values.get(&contig.id)) {
            if self.cov_sum.is_empty() {
                self.cov_sum.resize(row.len(), 0.0);
            }
            if self.cov_sum.len() == row.len() {
                for (dst, src) in self.cov_sum.iter_mut().zip(row.iter()) {
                    *dst += *src * w;
                }
            }
        }
        if let Some(ms) = markers.and_then(|t| t.values.get(&contig.id)) {
            self.markers.extend(ms.iter().cloned());
        }
    }

    fn absorb(&mut self, other: BinNode) {
        self.members.extend(other.members);
        self.bp += other.bp;
        self.gc_sum += other.gc_sum;
        self.sig_total += other.sig_total;
        for (dst, src) in self.sig_sum.iter_mut().zip(other.sig_sum.iter()) {
            *dst += *src;
        }
        if self.cov_sum.is_empty() && !other.cov_sum.is_empty() {
            self.cov_sum.resize(other.cov_sum.len(), 0.0);
        }
        if self.cov_sum.len() == other.cov_sum.len() {
            for (dst, src) in self.cov_sum.iter_mut().zip(other.cov_sum.iter()) {
                *dst += *src;
            }
        }
        self.markers.extend(other.markers);
    }
}

fn node_similarity(a: &BinNode, b: &BinNode, cfg: &BridgeBinConfig) -> SimilarityParts {
    let composition = (-hellinger_counts(a, b) / 0.24).exp();
    let gc_a = a.gc_sum / a.bp.max(1) as f64;
    let gc_b = b.gc_sum / b.bp.max(1) as f64;
    let gc = (-(gc_a - gc_b).abs() / 0.065).exp();
    let mut score = cfg.composition_weight * composition + cfg.gc_weight * gc;
    let mut weight = cfg.composition_weight + cfg.gc_weight;
    let coverage = if !a.cov_sum.is_empty() && a.cov_sum.len() == b.cov_sum.len() {
        let da = a.bp.max(1) as f64;
        let db = b.bp.max(1) as f64;
        let d = a
            .cov_sum
            .iter()
            .zip(b.cov_sum.iter())
            .map(|(x, y)| (((x / da) + 0.5) / ((y / db) + 0.5)).ln().abs())
            .sum::<f64>()
            / a.cov_sum.len() as f64;
        let cov = (-d / 0.80).exp();
        score += cfg.coverage_weight * cov;
        weight += cfg.coverage_weight;
        Some(cov)
    } else {
        None
    };
    SimilarityParts {
        combined: if weight <= f64::EPSILON { 0.0 } else { score / weight },
        composition,
        coverage,
    }
}

fn hellinger_counts(a: &BinNode, b: &BinNode) -> f64 {
    if a.sig_total <= 0.0 || b.sig_total <= 0.0 {
        return 1.0;
    }
    let sum = a
        .sig_sum
        .iter()
        .zip(b.sig_sum.iter())
        .map(|(x, y)| {
            let px = *x / a.sig_total;
            let py = *y / b.sig_total;
            (px.sqrt() - py.sqrt()).powi(2)
        })
        .sum::<f64>();
    (0.5 * sum).sqrt()
}

fn gc_fraction(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut n = 0usize;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc += 1;
                n += 1;
            }
            b'A' | b'T' => n += 1,
            _ => {}
        }
    }
    if n == 0 { 0.0 } else { gc as f64 / n as f64 }
}

fn canonical_5mer_counts(seq: &[u8]) -> ([f64; 1024], f64) {
    let mut counts = [0.0; 1024];
    let mut total = 0.0;
    for w in seq.windows(5) {
        if let (Some(fwd), Some(rc)) = (encode5(w, false), encode5(w, true)) {
            counts[fwd.min(rc)] += 1.0;
            total += 1.0;
        }
    }
    (counts, total)
}

fn encode5(w: &[u8], reverse_complement: bool) -> Option<usize> {
    let mut code = 0usize;
    if reverse_complement {
        for &b in w.iter().rev() {
            code = (code << 2) | base_code(b, true)?;
        }
    } else {
        for &b in w {
            code = (code << 2) | base_code(b, false)?;
        }
    }
    Some(code)
}

fn base_code(b: u8, complement: bool) -> Option<usize> {
    let x = match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    Some(if complement { 3 - x } else { x })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(id: &str, pattern: &str, n: usize) -> Contig {
        Contig { id: id.into(), seq: pattern.repeat(n).into_bytes() }
    }

    fn singleton_result(contigs: &[Contig]) -> BinningResult {
        BinningResult {
            assignments: contigs
                .iter()
                .enumerate()
                .map(|(i, c)| Assignment { contig_id: c.id.clone(), bin_index: Some(i), score: 1.0, length: c.seq.len() })
                .collect(),
            bins: contigs
                .iter()
                .enumerate()
                .map(|(i, c)| BinSummary { bin_index: i, contig_count: 1, total_bp: c.seq.len(), mean_gc: gc_fraction(&c.seq) })
                .collect(),
        }
    }

    #[test]
    fn reciprocal_reconciliation_merges_split_signatures() {
        let contigs = vec![
            c("a1", "AAAACAAAAGAAAATAAAC", 220),
            c("a2", "AAAGAAAACAAAATAAAA", 220),
            c("b1", "GGGCGGCCGCGGCGCCGC", 220),
            c("b2", "GCCGGCGCGGCCGGCGGC", 220),
        ];
        let coverage = CoverageTable {
            sample_names: vec!["s1".into(), "s2".into()],
            values: HashMap::from([
                ("a1".into(), vec![20.0, 10.0]),
                ("a2".into(), vec![20.5, 10.2]),
                ("b1".into(), vec![20.0, 10.0]),
                ("b2".into(), vec![19.8, 10.1]),
            ]),
        };
        let mut cfg = ReconcileConfig::default();
        cfg.merge_threshold = 0.55;
        cfg.min_composition_score = 0.45;
        cfg.same_coverage_min_composition = 0.70;
        let (result, stats) = reconcile_bins(
            &contigs,
            Some(&coverage),
            None,
            singleton_result(&contigs),
            &BridgeBinConfig::default(),
            &cfg,
        );
        assert_eq!(stats.final_bins, 2);
        let a = result.assignments.iter().find(|x| x.contig_id == "a1").unwrap().bin_index;
        let a2 = result.assignments.iter().find(|x| x.contig_id == "a2").unwrap().bin_index;
        let b = result.assignments.iter().find(|x| x.contig_id == "b1").unwrap().bin_index;
        assert_eq!(a, a2);
        assert_ne!(a, b);
    }

    #[test]
    fn shared_single_copy_marker_blocks_merge() {
        let contigs = vec![
            c("x1", "ACGTACGTAAAACCCCGGGG", 220),
            c("x2", "ACGTACGTAAAACCCCGGGA", 220),
        ];
        let markers = MarkerTable {
            values: HashMap::from([
                ("x1".into(), HashSet::from(["SCG001".into()])),
                ("x2".into(), HashSet::from(["SCG001".into()])),
            ]),
        };
        let mut cfg = ReconcileConfig::default();
        cfg.merge_threshold = 0.1;
        cfg.min_composition_score = 0.0;
        cfg.min_coverage_score = 0.0;
        cfg.merge_margin = 0.0;
        let (_, stats) = reconcile_bins(
            &contigs,
            None,
            Some(&markers),
            singleton_result(&contigs),
            &BridgeBinConfig::default(),
            &cfg,
        );
        assert_eq!(stats.final_bins, 2);
        assert!(stats.marker_blocked_pairs > 0);
    }
}
