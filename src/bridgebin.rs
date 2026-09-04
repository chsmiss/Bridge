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
struct Feature {
    id: String,
    len: usize,
    gc: f64,
    tnf: [f64; 256],
    cov: Vec<f64>,
}

#[derive(Clone, Debug)]
struct BinState {
    members: Vec<usize>,
    bp: usize,
    gc_sum: f64,
    tnf_sum: [f64; 256],
    cov_sum: Vec<f64>,
}

impl BinState {
    fn new(index: usize, f: &Feature) -> Self {
        let w = f.len as f64;
        let mut tnf_sum = [0.0; 256];
        for (dst, src) in tnf_sum.iter_mut().zip(f.tnf.iter()) {
            *dst = *src * w;
        }
        Self {
            members: vec![index],
            bp: f.len,
            gc_sum: f.gc * w,
            tnf_sum,
            cov_sum: f.cov.iter().map(|v| v * w).collect(),
        }
    }

    fn add(&mut self, index: usize, f: &Feature) {
        let w = f.len as f64;
        self.members.push(index);
        self.bp += f.len;
        self.gc_sum += f.gc * w;
        for (dst, src) in self.tnf_sum.iter_mut().zip(f.tnf.iter()) {
            *dst += *src * w;
        }
        if self.cov_sum.is_empty() && !f.cov.is_empty() {
            self.cov_sum.resize(f.cov.len(), 0.0);
        }
        if self.cov_sum.len() == f.cov.len() {
            for (dst, src) in self.cov_sum.iter_mut().zip(f.cov.iter()) {
                *dst += *src * w;
            }
        }
    }

    fn centroid(&self) -> Feature {
        let d = self.bp.max(1) as f64;
        let mut tnf = [0.0; 256];
        for (dst, src) in tnf.iter_mut().zip(self.tnf_sum.iter()) {
            *dst = *src / d;
        }
        Feature {
            id: String::new(),
            len: self.bp,
            gc: self.gc_sum / d,
            tnf,
            cov: self.cov_sum.iter().map(|v| *v / d).collect(),
        }
    }
}

pub fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<Vec<Contig>> {
    let reader = BufReader::new(File::open(path)?);
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let mut id: Option<String> = None;
    let mut seq = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() {
            continue;
        }
        if let Some(header) = s.strip_prefix('>') {
            if let Some(old_id) = id.take() {
                out.push(Contig {
                    id: old_id,
                    seq: std::mem::take(&mut seq),
                });
            }
            let new_id = header.split_whitespace().next().unwrap_or("");
            if new_id.is_empty() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "empty FASTA id"));
            }
            if !seen.insert(new_id.to_string()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("duplicate FASTA id '{}' at line {}", new_id, line_no + 1),
                ));
            }
            id = Some(new_id.to_string());
        } else {
            if id.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sequence before first FASTA header",
                ));
            }
            seq.extend(s.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(last_id) = id {
        out.push(Contig { id: last_id, seq });
    }
    if out.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTA contains no contigs"));
    }
    Ok(out)
}

pub fn read_coverage_table<P: AsRef<Path>>(path: P) -> io::Result<CoverageTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut sample_names = Vec::new();
    let mut expected: Option<usize> = None;
    let mut values = HashMap::new();
    let mut first = true;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "coverage row needs >=2 columns"));
        }
        if first && fields[1].parse::<f64>().is_err() {
            sample_names = fields[1..].iter().map(|x| (*x).to_string()).collect();
            expected = Some(sample_names.len());
            first = false;
            continue;
        }
        first = false;
        let n = fields.len() - 1;
        if let Some(m) = expected {
            if m != n {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("coverage row {} has {} columns; expected {}", line_no + 1, n, m),
                ));
            }
        } else {
            expected = Some(n);
            sample_names = (1..=n).map(|i| format!("sample{}", i)).collect();
        }
        let mut row = Vec::with_capacity(n);
        for raw in &fields[1..] {
            let v = raw.parse::<f64>().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, format!("invalid coverage '{}'", raw))
            })?;
            if !v.is_finite() || v < 0.0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "coverage must be finite and >=0"));
            }
            row.push(v);
        }
        if values.insert(fields[0].to_string(), row).is_some() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "duplicate coverage contig"));
        }
    }
    if values.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "coverage table has no data"));
    }
    Ok(CoverageTable {
        sample_names,
        values,
    })
}

pub fn bin_contigs(
    contigs: &[Contig],
    coverage: Option<&CoverageTable>,
    cfg: &BridgeBinConfig,
) -> BinningResult {
    let features: Vec<Feature> = contigs
        .iter()
        .filter(|c| c.seq.len() >= cfg.min_contig_len)
        .map(|c| feature(c, coverage))
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

    let mut seeds: Vec<usize> = features
        .iter()
        .enumerate()
        .filter_map(|(i, f)| (f.len >= cfg.seed_min_len).then_some(i))
        .collect();
    if seeds.is_empty() {
        seeds = (0..features.len()).collect();
    }
    seeds.sort_by(|&a, &b| features[b].len.cmp(&features[a].len));
    let seed_set: HashSet<usize> = seeds.iter().copied().collect();

    let mut bins: Vec<BinState> = Vec::new();
    let mut assigned: Vec<Option<(usize, f64)>> = vec![None; features.len()];
    for &i in &seeds {
        let (best, score, _) = best_two(&features[i], &bins, cfg);
        if let Some(bin) = best.filter(|_| score >= cfg.join_threshold) {
            bins[bin].add(i, &features[i]);
            assigned[i] = Some((bin, score));
        } else {
            let bin = bins.len();
            bins.push(BinState::new(i, &features[i]));
            assigned[i] = Some((bin, 1.0));
        }
    }

    let mut rescue: Vec<usize> = (0..features.len())
        .filter(|i| !seed_set.contains(i))
        .collect();
    rescue.sort_by(|&a, &b| features[b].len.cmp(&features[a].len));
    for i in rescue {
        let (best, score, second) = best_two(&features[i], &bins, cfg);
        let margin_ok = bins.len() <= 1 || score - second >= cfg.rescue_margin;
        if let Some(bin) = best.filter(|_| score >= cfg.rescue_threshold && margin_ok) {
            bins[bin].add(i, &features[i]);
            assigned[i] = Some((bin, score));
        }
    }

    let mut order: Vec<usize> = (0..bins.len()).collect();
    order.sort_by(|&a, &b| bins[b].bp.cmp(&bins[a].bp).then_with(|| a.cmp(&b)));
    let mut remap = vec![0usize; bins.len()];
    for (new, old) in order.iter().copied().enumerate() {
        remap[old] = new;
    }
    let by_id: HashMap<&str, usize> = features
        .iter()
        .enumerate()
        .map(|(i, f)| (f.id.as_str(), i))
        .collect();

    let assignments = contigs
        .iter()
        .map(|c| match by_id.get(c.id.as_str()).and_then(|&i| assigned[i]) {
            Some((old, score)) => Assignment {
                contig_id: c.id.clone(),
                bin_index: Some(remap[old]),
                score,
                length: c.seq.len(),
            },
            None => Assignment {
                contig_id: c.id.clone(),
                bin_index: None,
                score: 0.0,
                length: c.seq.len(),
            },
        })
        .collect();

    let mut summaries = Vec::with_capacity(bins.len());
    for old in order {
        let c = bins[old].centroid();
        summaries.push(BinSummary {
            bin_index: remap[old],
            contig_count: bins[old].members.len(),
            total_bp: bins[old].bp,
            mean_gc: c.gc,
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
    let out = out_dir.as_ref();
    let bins_dir = out.join("bins");
    fs::create_dir_all(&bins_dir)?;
    let assignment: HashMap<&str, Option<usize>> = result
        .assignments
        .iter()
        .map(|a| (a.contig_id.as_str(), a.bin_index))
        .collect();
    let mut writers = HashMap::new();
    for b in &result.bins {
        writers.insert(
            b.bin_index,
            BufWriter::new(File::create(bins_dir.join(format!("bin_{:04}.fa", b.bin_index + 1)))?),
        );
    }
    let mut unbinned = if emit_unbinned {
        Some(BufWriter::new(File::create(out.join("unbinned.fa"))?))
    } else {
        None
    };
    for c in contigs {
        match assignment.get(c.id.as_str()).copied().flatten() {
            Some(bin) => write_fasta(writers.get_mut(&bin).expect("known bin"), c)?,
            None => {
                if let Some(w) = unbinned.as_mut() {
                    write_fasta(w, c)?;
                }
            }
        }
    }
    let mut aout = BufWriter::new(File::create(out.join("assignments.tsv"))?);
    writeln!(aout, "contig\tbin\tlength\tscore")?;
    for a in &result.assignments {
        let name = a
            .bin_index
            .map(|i| format!("bin_{:04}", i + 1))
            .unwrap_or_else(|| "unbinned".to_string());
        writeln!(aout, "{}\t{}\t{}\t{:.6}", a.contig_id, name, a.length, a.score)?;
    }
    let mut bout = BufWriter::new(File::create(out.join("bins.tsv"))?);
    writeln!(bout, "bin\tcontigs\ttotal_bp\tmean_gc")?;
    for b in &result.bins {
        writeln!(
            bout,
            "bin_{:04}\t{}\t{}\t{:.6}",
            b.bin_index + 1,
            b.contig_count,
            b.total_bp,
            b.mean_gc
        )?;
    }
    Ok(())
}

fn write_fasta<W: Write>(w: &mut W, c: &Contig) -> io::Result<()> {
    writeln!(w, ">{}", c.id)?;
    for chunk in c.seq.chunks(80) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

fn feature(c: &Contig, coverage: Option<&CoverageTable>) -> Feature {
    Feature {
        id: c.id.clone(),
        len: c.seq.len(),
        gc: gc_fraction(&c.seq),
        tnf: canonical_tnf(&c.seq),
        cov: coverage
            .and_then(|x| x.values.get(&c.id))
            .cloned()
            .unwrap_or_default(),
    }
}

fn best_two(f: &Feature, bins: &[BinState], cfg: &BridgeBinConfig) -> (Option<usize>, f64, f64) {
    let mut best = None;
    let mut s1 = f64::NEG_INFINITY;
    let mut s2 = f64::NEG_INFINITY;
    for (i, bin) in bins.iter().enumerate() {
        let s = similarity(f, &bin.centroid(), cfg);
        if s > s1 {
            s2 = s1;
            s1 = s;
            best = Some(i);
        } else if s > s2 {
            s2 = s;
        }
    }
    if best.is_none() {
        (None, 0.0, 0.0)
    } else {
        (best, s1, if s2.is_finite() { s2 } else { 0.0 })
    }
}

fn similarity(a: &Feature, b: &Feature, cfg: &BridgeBinConfig) -> f64 {
    let comp = (-hellinger(&a.tnf, &b.tnf) / 0.30).exp();
    let gc = (-(a.gc - b.gc).abs() / 0.08).exp();
    let mut score = cfg.composition_weight * comp + cfg.gc_weight * gc;
    let mut weight = cfg.composition_weight + cfg.gc_weight;
    if !a.cov.is_empty() && a.cov.len() == b.cov.len() {
        let d = a
            .cov
            .iter()
            .zip(b.cov.iter())
            .map(|(x, y)| ((x + 0.5) / (y + 0.5)).ln().abs())
            .sum::<f64>()
            / a.cov.len() as f64;
        score += cfg.coverage_weight * (-d / 0.85).exp();
        weight += cfg.coverage_weight;
    }
    if weight <= f64::EPSILON {
        0.0
    } else {
        score / weight
    }
}

fn hellinger(a: &[f64; 256], b: &[f64; 256]) -> f64 {
    (0.5
        * a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x.sqrt() - y.sqrt()).powi(2))
            .sum::<f64>())
    .sqrt()
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
    if n == 0 {
        0.0
    } else {
        gc as f64 / n as f64
    }
}

fn canonical_tnf(seq: &[u8]) -> [f64; 256] {
    let mut counts = [0.0; 256];
    let mut total = 0.0;
    for w in seq.windows(4) {
        if let (Some(fwd), Some(rc)) = (encode4(w, false), encode4(w, true)) {
            counts[fwd.min(rc)] += 1.0;
            total += 1.0;
        }
    }
    if total > 0.0 {
        for x in &mut counts {
            *x /= total;
        }
    }
    counts
}

fn encode4(w: &[u8], reverse_complement: bool) -> Option<usize> {
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
    use std::time::{SystemTime, UNIX_EPOCH};

    fn contig(id: &str, pattern: &str, repeats: usize) -> Contig {
        Contig {
            id: id.to_string(),
            seq: pattern.repeat(repeats).into_bytes(),
        }
    }

    #[test]
    fn tnf_is_reverse_complement_invariant() {
        let a = canonical_tnf(b"AAAACGTTCCGA");
        let b = canonical_tnf(b"TCGGAACGTTTT");
        for i in 0..256 {
            assert!((a[i] - b[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn coverage_header_is_supported() {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("bridgebin_cov_{}_{}.tsv", process_id(), stamp));
        fs::write(&path, "contig\ts1\ts2\na\t10\t2\nb\t9.5\t2.5\n").unwrap();
        let table = read_coverage_table(&path).unwrap();
        fs::remove_file(path).ok();
        assert_eq!(table.sample_names, vec!["s1", "s2"]);
        assert_eq!(table.values["a"], vec![10.0, 2.0]);
    }

    fn process_id() -> u32 {
        std::process::id()
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
        cfg.join_threshold = 0.60;
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
        let result = bin_contigs(&contigs, None, &BridgeBinConfig::default());
        assert!(result
            .assignments
            .iter()
            .find(|a| a.contig_id == "tiny")
            .unwrap()
            .bin_index
            .is_none());
    }
}
