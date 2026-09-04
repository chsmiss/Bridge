use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Contig {
    pub id: String,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct CoverageTable {
    pub sample_names: Vec<String>,
    pub values: HashMap<String, Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct BridgeBinConfig {
    pub min_contig_len: usize,
    pub seed_min_len: usize,
    pub join_threshold: f64,
    pub rescue_threshold: f64,
    pub rescue_margin: f64,
    pub composition_weight: f64,
    pub coverage_weight: f64,
    pub gc_weight: f64,
}

impl Default for BridgeBinConfig {
    fn default() -> Self {
        Self {
            min_contig_len: 1_500,
            seed_min_len: 2_500,
            join_threshold: 0.76,
            rescue_threshold: 0.70,
            rescue_margin: 0.025,
            composition_weight: 0.45,
            coverage_weight: 0.50,
            gc_weight: 0.05,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Assignment {
    pub contig_id: String,
    pub bin_index: Option<usize>,
    pub score: f64,
    pub length: usize,
}

#[derive(Clone, Debug)]
pub struct BinSummary {
    pub bin_index: usize,
    pub contig_count: usize,
    pub total_bp: usize,
    pub mean_gc: f64,
}

#[derive(Clone, Debug)]
pub struct BinningResult {
    pub assignments: Vec<Assignment>,
    pub bins: Vec<BinSummary>,
}

#[derive(Clone, Debug)]
struct FeatureVector {
    id: String,
    length: usize,
    gc: f64,
    composition: [f64; 256],
    coverage: Vec<f64>,
}

#[derive(Clone, Debug)]
struct BinState {
    members: Vec<usize>,
    total_bp: usize,
    gc_sum: f64,
    composition_sum: [f64; 256],
    coverage_sum: Vec<f64>,
}

impl BinState {
    fn new(feature_index: usize, f: &FeatureVector) -> Self {
        let weight = f.length as f64;
        let mut composition_sum = [0.0; 256];
        for (dst, src) in composition_sum.iter_mut().zip(f.composition.iter()) {
            *dst = *src * weight;
        }
        let coverage_sum = f.coverage.iter().map(|v| v * weight).collect();
        Self {
            members: vec![feature_index],
            total_bp: f.length,
            gc_sum: f.gc * weight,
            composition_sum,
            coverage_sum,
        }
    }

    fn add(&mut self, feature_index: usize, f: &FeatureVector) {
        let weight = f.length as f64;
        self.members.push(feature_index);
        self.total_bp += f.length;
        self.gc_sum += f.gc * weight;
        for (dst, src) in self.composition_sum.iter_mut().zip(f.composition.iter()) {
            *dst += *src * weight;
        }
        if self.coverage_sum.is_empty() && !f.coverage.is_empty() {
            self.coverage_sum.resize(f.coverage.len(), 0.0);
        }
        if self.coverage_sum.len() == f.coverage.len() {
            for (dst, src) in self.coverage_sum.iter_mut().zip(f.coverage.iter()) {
                *dst += *src * weight;
            }
        }
    }

    fn centroid(&self) -> FeatureVector {
        let denom = self.total_bp.max(1) as f64;
        let mut composition = [0.0; 256];
        for (dst, src) in composition.iter_mut().zip(self.composition_sum.iter()) {
            *dst = *src / denom;
        }
        let coverage = self.coverage_sum.iter().map(|v| *v / denom).collect();
        FeatureVector {
            id: String::new(),
            length: self.total_bp,
            gc: self.gc_sum / denom,
            composition,
            coverage,
        }
    }
}

pub fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<Vec<Contig>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut contigs = Vec::new();
    let mut seen = HashSet::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                contigs.push(Contig {
                    id,
                    seq: std::mem::take(&mut current_seq),
                });
            }
            let id = rest
                .split_whitespace()
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty FASTA header"))?;
            if id.is_empty() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "empty FASTA id"));
            }
            if !seen.insert(id.to_string()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("duplicate FASTA id '{}' at line {}", id, line_no + 1),
                ));
            }
            current_id = Some(id.to_string());
        } else {
            if current_id.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("sequence before first FASTA header at line {}", line_no + 1),
                ));
            }
            current_seq.extend(trimmed.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }

    if let Some(id) = current_id.take() {
        contigs.push(Contig { id, seq: current_seq });
    }
    if contigs.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTA contains no contigs"));
    }
    Ok(contigs)
}

pub fn read_coverage_table<P: AsRef<Path>>(path: P) -> io::Result<CoverageTable> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows: Vec<(String, Vec<f64>)> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();
    let mut expected_cols: Option<usize> = None;
    let mut first_data_seen = false;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("coverage row {} has fewer than 2 columns", line_no + 1),
            ));
        }

        if !first_data_seen && fields[1].parse::<f64>().is_err() {
            sample_names = fields[1..].iter().map(|s| (*s).to_string()).collect();
            expected_cols = Some(sample_names.len());
            first_data_seen = true;
            continue;
        }
        first_data_seen = true;

        let n = fields.len() - 1;
        match expected_cols {
            Some(m) if m != n => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "coverage row {} has {} sample columns; expected {}",
                        line_no + 1,
                        n,
                        m
                    ),
                ));
            }
            None => {
                expected_cols = Some(n);
                sample_names = (1..=n).map(|i| format!("sample{}", i)).collect();
            }
            _ => {}
        }

        let mut values = Vec::with_capacity(n);
        for raw in &fields[1..] {
            let value = raw.parse::<f64>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid coverage value '{}' at line {}", raw, line_no + 1),
                )
            })?;
            if !value.is_finite() || value < 0.0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("coverage must be finite and >= 0 at line {}", line_no + 1),
                ));
            }
            values.push(value);
        }
        rows.push((fields[0].to_string(), values));
    }

    if rows.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "coverage table contains no data rows",
        ));
    }

    let mut map = HashMap::with_capacity(rows.len());
    for (id, values) in rows {
        if map.insert(id.clone(), values).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate contig '{}' in coverage table", id),
            ));
        }
    }

    Ok(CoverageTable {
        sample_names,
        values: map,
    })
}

pub fn bin_contigs(
    contigs: &[Contig],
    coverage: Option<&CoverageTable>,
    config: &BridgeBinConfig,
) -> BinningResult {
    let features: Vec<FeatureVector> = contigs
        .iter()
        .filter(|c| c.seq.len() >= config.min_contig_len)
        .map(|c| feature_for_contig(c, coverage))
        .collect();

    if features.is_empty() {
        return BinningResult {
            assignments: contigs
                .iter()
                .map(|c| Assignment {
                    contig_id: c.id.clone(),
                    bin_index: None,
                    score: 0.0,
                    length: c.seq.len(),
                })
                .collect(),
            bins: Vec::new(),
        };
    }

    let mut seed_indices: Vec<usize> = features
        .iter()
        .enumerate()
        .filter_map(|(i, f)| (f.length >= config.seed_min_len).then_some(i))
        .collect();
    if seed_indices.is_empty() {
        seed_indices = (0..features.len()).collect();
    }
    seed_indices.sort_by(|&a, &b| features[b].length.cmp(&features[a].length));

    let seed_set: HashSet<usize> = seed_indices.iter().copied().collect();
    let mut bins: Vec<BinState> = Vec::new();
    let mut assignment_by_feature: Vec<Option<(usize, f64)>> = vec![None; features.len()];

    for &idx in &seed_indices {
        let (best_bin, best_score, _) = best_two_bins(&features[idx], &bins, config);
        if let Some(bin_idx) = best_bin.filter(|_| best_score >= config.join_threshold) {
            bins[bin_idx].add(idx, &features[idx]);
            assignment_by_feature[idx] = Some((bin_idx, best_score));
        } else {
            let bin_idx = bins.len();
            bins.push(BinState::new(idx, &features[idx]));
            assignment_by_feature[idx] = Some((bin_idx, 1.0));
        }
    }

    let mut rescue_indices: Vec<usize> = (0..features.len())
        .filter(|i| !seed_set.contains(i))
        .collect();
    rescue_indices.sort_by(|&a, &b| features[b].length.cmp(&features[a].length));
    for idx in rescue_indices {
        let (best_bin, best_score, second_score) = best_two_bins(&features[idx], &bins, config);
        let margin_ok = bins.len() <= 1 || best_score - second_score >= config.rescue_margin;
        if let Some(bin_idx) = best_bin.filter(|_| best_score >= config.rescue_threshold && margin_ok) {
            bins[bin_idx].add(idx, &features[idx]);
            assignment_by_feature[idx] = Some((bin_idx, best_score));
        }
    }

    let mut order: Vec<usize> = (0..bins.len()).collect();
    order.sort_by(|&a, &b| bins[b].total_bp.cmp(&bins[a].total_bp).then_with(|| a.cmp(&b)));
    let mut remap = vec![0usize; bins.len()];
    for (new_idx, old_idx) in order.iter().copied().enumerate() {
        remap[old_idx] = new_idx;
    }

    let feature_index_by_id: HashMap<&str, usize> = features
        .iter()
        .enumerate()
        .map(|(i, f)| (f.id.as_str(), i))
        .collect();

    let assignments = contigs
        .iter()
        .map(|c| {
            if let Some(&fi) = feature_index_by_id.get(c.id.as_str()) {
                if let Some((old_bin, score)) = assignment_by_feature[fi] {
                    Assignment {
                        contig_id: c.id.clone(),
                        bin_index: Some(remap[old_bin]),
                        score,
                        length: c.seq.len(),
                    }
                } else {
                    Assignment {
                        contig_id: c.id.clone(),
                        bin_index: None,
                        score: 0.0,
                        length: c.seq.len(),
                    }
                }
            } else {
                Assignment {
                    contig_id: c.id.clone(),
                    bin_index: None,
                    score: 0.0,
                    length: c.seq.len(),
                }
            }
        })
        .collect();

    let mut summaries = Vec::with_capacity(bins.len());
    for old_idx in order {
        let state = &bins[old_idx];
        let centroid = state.centroid();
        summaries.push(BinSummary {
            bin_index: remap[old_idx],
            contig_count: state.members.len(),
            total_bp: state.total_bp,
            mean_gc: centroid.gc,
        });
    }
    summaries.sort_by_key(|b| b.bin_index);

    BinningResult {
        assignments,
        bins: summaries,
    }
}

pub fn write_outputs<P: AsRef<Path>>(
    contigs: &[Contig],
    result: &BinningResult,
    out_dir: P,
    emit_unbinned: bool,
) -> io::Result<()> {
    let out_dir = out_dir.as_ref();
    fs::create_dir_all(out_dir)?;
    let bins_dir = out_dir.join("bins");
    fs::create_dir_all(&bins_dir)?;

    let mut assignment_map: HashMap<&str, Option<usize>> = HashMap::new();
    for a in &result.assignments {
        assignment_map.insert(a.contig_id.as_str(), a.bin_index);
    }

    let mut writers: HashMap<usize, BufWriter<File>> = HashMap::new();
    for summary in &result.bins {
        let path = bins_dir.join(format!("bin_{:04}.fa", summary.bin_index + 1));
        writers.insert(summary.bin_index, BufWriter::new(File::create(path)?));
    }
    let mut unbinned = if emit_unbinned {
        Some(BufWriter::new(File::create(out_dir.join("unbinned.fa"))?))
    } else {
        None
    };

    for contig in contigs {
        let target = assignment_map.get(contig.id.as_str()).copied().flatten();
        match target {
            Some(bin_idx) => {
                if let Some(writer) = writers.get_mut(&bin_idx) {
                    write_fasta_record(writer, contig)?;
                }
            }
            None => {
                if let Some(writer) = unbinned.as_mut() {
                    write_fasta_record(writer, contig)?;
                }
            }
        }
    }

    let mut assignments = BufWriter::new(File::create(out_dir.join("assignments.tsv"))?);
    writeln!(assignments, "contig\tbin\tlength\tscore")?;
    for a in &result.assignments {
        let bin_name = a
            .bin_index
            .map(|i| format!("bin_{:04}", i + 1))
            .unwrap_or_else(|| "unbinned".to_string());
        writeln!(assignments, "{}\t{}\t{}\t{:.6}", a.contig_id, bin_name, a.length, a.score)?;
    }

    let mut summary = BufWriter::new(File::create(out_dir.join("bins.tsv"))?);
    writeln!(summary, "bin\tcontigs\ttotal_bp\tmean_gc")?;
    for b in &result.bins {
        writeln!(
            summary,
            "bin_{:04}\t{}\t{}\t{:.6}",
            b.bin_index + 1,
            b.contig_count,
            b.total_bp,
            b.mean_gc
        )?;
    }
    Ok(())
}

fn write_fasta_record<W: Write>(writer: &mut W, contig: &Contig) -> io::Result<()> {
    writeln!(writer, ">{}", contig.id)?;
    for chunk in contig.seq.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn feature_for_contig(contig: &Contig, coverage: Option<&CoverageTable>) -> FeatureVector {
    FeatureVector {
        id: contig.id.clone(),
        length: contig.seq.len(),
        gc: gc_fraction(&contig.seq),
        composition: canonical_tnf(&contig.seq),
        coverage: coverage
            .and_then(|t| t.values.get(&contig.id))
            .cloned()
            .unwrap_or_default(),
    }
}

fn best_two_bins(
    feature: &FeatureVector,
    bins: &[BinState],
    config: &BridgeBinConfig,
) -> (Option<usize>, f64, f64) {
    let mut best_bin = None;
    let mut best_score = f64::NEG_INFINITY;
    let mut second_score = f64::NEG_INFINITY;
    for (bin_idx, bin) in bins.iter().enumerate() {
        let centroid = bin.centroid();
        let score = similarity(feature, &centroid, config);
        if score > best_score {
            second_score = best_score;
            best_score = score;
            best_bin = Some(bin_idx);
        } else if score > second_score {
            second_score = score;
        }
    }
    if best_bin.is_none() {
        (None, 0.0, 0.0)
    } else {
        let second = if second_score.is_finite() { second_score } else { 0.0 };
        (best_bin, best_score, second)
    }
}

fn similarity(a: &FeatureVector, b: &FeatureVector, config: &BridgeBinConfig) -> f64 {
    let comp_distance = hellinger(&a.composition, &b.composition);
    let comp_similarity = (-comp_distance / 0.30).exp();
    let gc_similarity = (-(a.gc - b.gc).abs() / 0.08).exp();

    let mut weighted = config.composition_weight * comp_similarity + config.gc_weight * gc_similarity;
    let mut total_weight = config.composition_weight + config.gc_weight;

    if !a.coverage.is_empty() && a.coverage.len() == b.coverage.len() {
        let cov_distance = a
            .coverage
            .iter()
            .zip(b.coverage.iter())
            .map(|(x, y)| ((x + 0.5) / (y + 0.5)).ln().abs())
            .sum::<f64>()
            / a.coverage.len() as f64;
        let cov_similarity = (-cov_distance / 0.85).exp();
        weighted += config.coverage_weight * cov_similarity;
        total_weight += config.coverage_weight;
    }

    if total_weight <= f64::EPSILON {
        0.0
    } else {
        weighted / total_weight
    }
}

fn hellinger(a: &[f64; 256], b: &[f64; 256]) -> f64 {
    let sum = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x.sqrt() - y.sqrt()).powi(2))
        .sum::<f64>();
    (0.5 * sum).sqrt()
}

fn gc_fraction(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut valid = 0usize;
    for &b in seq {
        match b.to_ascii_uppercase() {
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

fn canonical_tnf(seq: &[u8]) -> [f64; 256] {
    let mut counts = [0.0f64; 256];
    let mut total = 0.0f64;
    if seq.len() < 4 {
        return counts;
    }
    for window in seq.windows(4) {
        if let (Some(code), Some(rc)) = (encode_4mer(window), encode_revcomp_4mer(window)) {
            counts[code.min(rc)] += 1.0;
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

fn encode_4mer(window: &[u8]) -> Option<usize> {
    let mut code = 0usize;
    for &b in window {
        code = (code << 2) | base_code(b)?;
    }
    Some(code)
}

fn encode_revcomp_4mer(window: &[u8]) -> Option<usize> {
    let mut code = 0usize;
    for &b in window.iter().rev() {
        code = (code << 2) | complement_code(b)?;
    }
    Some(code)
}

fn base_code(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn complement_code(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(3),
        b'C' => Some(2),
        b'G' => Some(1),
        b'T' => Some(0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn contig(id: &str, pattern: &str, repeats: usize) -> Contig {
        Contig {
            id: id.to_string(),
            seq: pattern.repeat(repeats).into_bytes(),
        }
    }

    #[test]
    fn canonical_tnf_is_reverse_complement_invariant() {
        let x = canonical_tnf(b"AAAACGTTCCGA");
        let rc = canonical_tnf(b"TCGGAACGTTTT");
        for i in 0..256 {
            assert!((x[i] - rc[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn coverage_table_accepts_header() {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "bridgebin_cov_{}_{}.tsv",
            std::process::id(),
            stamp
        ));
        fs::write(&path, "contig\ts1\ts2\na\t10\t2\nb\t9.5\t2.5\n").unwrap();
        let table = read_coverage_table(&path).unwrap();
        fs::remove_file(path).ok();
        assert_eq!(table.sample_names, vec!["s1", "s2"]);
        assert_eq!(table.values["a"], vec![10.0, 2.0]);
    }

    #[test]
    fn separates_two_synthetic_genome_signatures() {
        let contigs = vec![
            contig("a1", "AAAACAAAAGAAAATAAAC", 180),
            contig("a2", "AAAGAAAACAAAATAAAA", 180),
            contig("b1", "GGGCGGCCGCGGCGCCGC", 180),
            contig("b2", "GCCGGCGCGGCCGGCGGC", 180),
        ];
        let coverage = CoverageTable {
            sample_names: vec!["s1".into(), "s2".into()],
            values: HashMap::from([
                ("a1".into(), vec![30.0, 5.0]),
                ("a2".into(), vec![29.0, 5.5]),
                ("b1".into(), vec![7.0, 25.0]),
                ("b2".into(), vec![7.5, 24.0]),
            ]),
        };
        let mut cfg = BridgeBinConfig::default();
        cfg.min_contig_len = 1_000;
        cfg.seed_min_len = 2_000;
        let result = bin_contigs(&contigs, Some(&coverage), &cfg);
        let bins: HashMap<&str, usize> = result
            .assignments
            .iter()
            .filter_map(|a| a.bin_index.map(|b| (a.contig_id.as_str(), b)))
            .collect();
        assert_eq!(bins["a1"], bins["a2"]);
        assert_eq!(bins["b1"], bins["b2"]);
        assert_ne!(bins["a1"], bins["b1"]);
        assert_eq!(result.bins.len(), 2);
    }

    #[test]
    fn short_contigs_remain_unbinned() {
        let contigs = vec![contig("tiny", "ACGT", 10), contig("long", "ACGT", 800)];
        let cfg = BridgeBinConfig::default();
        let result = bin_contigs(&contigs, None, &cfg);
        let tiny = result
            .assignments
            .iter()
            .find(|a| a.contig_id == "tiny")
            .unwrap();
        assert!(tiny.bin_index.is_none());
    }
}
