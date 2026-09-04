use crate::bridgebin::{Assignment, BinSummary, BinningResult, Contig, CoverageTable};
use crate::bridgebin_reconcile::MarkerTable;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct BridgeBinV2Config {
    pub min_contig_len: usize,
    pub max_neighbors: usize,
    pub min_component_bp: usize,
    pub core_min_attraction: f64,
    pub core_max_repulsion: f64,
    pub component_min_attraction: f64,
    pub rescue_min_attraction: f64,
    pub rescue_margin: f64,
    pub max_gc_delta: f64,
    pub min_component_composition: f64,
    pub min_component_coverage: f64,
    pub taxonomy_confidence: f64,
    pub hard_marker_veto: bool,
}

impl Default for BridgeBinV2Config {
    fn default() -> Self {
        Self {
            min_contig_len: 1_500,
            max_neighbors: 64,
            min_component_bp: 20_000,
            core_min_attraction: 0.80,
            core_max_repulsion: 0.20,
            component_min_attraction: 0.72,
            rescue_min_attraction: 0.76,
            rescue_margin: 0.05,
            max_gc_delta: 0.075,
            min_component_composition: 0.50,
            min_component_coverage: 0.48,
            taxonomy_confidence: 0.90,
            hard_marker_veto: true,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct BioFeature {
    pub taxonomy: Option<String>,
    pub taxonomy_confidence: f64,
    pub gene_profile: Vec<f64>,
    pub gene_confidence: f64,
    pub protein_embedding: Vec<f64>,
    pub protein_confidence: f64,
}

#[derive(Clone, Debug, Default)]
pub struct BioFeatureTable {
    pub values: HashMap<String, BioFeature>,
}

#[derive(Clone, Debug, Default)]
pub struct LinkEvidence {
    pub must_link: f64,
    pub cannot_link: f64,
    pub source: String,
}

#[derive(Clone, Debug, Default)]
pub struct LinkTable {
    pub values: HashMap<(String, String), LinkEvidence>,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct PairEvidence {
    pub attraction: f64,
    pub repulsion: f64,
    pub composition: f64,
    pub coverage: Option<f64>,
    pub gc: f64,
    pub gene: Option<f64>,
    pub protein: Option<f64>,
    pub marker_conflict: bool,
    pub taxonomy_conflict: bool,
    pub external_cannot_link: bool,
}

impl PairEvidence {
    fn hard_conflict(&self) -> bool {
        self.marker_conflict || self.taxonomy_conflict || self.external_cannot_link
    }
}

#[derive(Clone, Debug, Default)]
pub struct BridgeBinV2Stats {
    pub eligible_contigs: usize,
    pub candidate_edges: usize,
    pub accepted_core_edges: usize,
    pub marker_blocked_edges: usize,
    pub taxonomy_blocked_edges: usize,
    pub external_blocked_edges: usize,
    pub core_bins: usize,
    pub rescued_components: usize,
    pub unbinned_contigs: usize,
}

#[derive(Clone, Debug)]
struct Feature {
    contig_index: usize,
    id: String,
    len: usize,
    gc: f64,
    kmer: [f64; 1024],
    coverage: Vec<f64>,
    markers: HashSet<String>,
    bio: BioFeature,
}

#[derive(Clone, Debug)]
struct CandidateEdge {
    left: usize,
    right: usize,
    evidence: PairEvidence,
}

#[derive(Clone, Debug)]
struct Component {
    members: Vec<usize>,
    bp: usize,
    gc_sum: f64,
    kmer_sum: [f64; 1024],
    coverage_sum: Vec<f64>,
    markers: HashSet<String>,
    taxonomies: Vec<(String, f64)>,
    gene_sum: Vec<f64>,
    gene_weight: f64,
    protein_sum: Vec<f64>,
    protein_weight: f64,
}

pub fn read_bio_feature_table<P: AsRef<Path>>(path: P) -> io::Result<BioFeatureTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty bio feature table"))?;
    let columns: Vec<&str> = header.trim().split('\t').collect();
    let find = |names: &[&str]| -> Option<usize> {
        columns
            .iter()
            .position(|column| names.iter().any(|name| column.eq_ignore_ascii_case(name)))
    };
    let contig_col = find(&["contig", "contig_id", "sequence"])
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "bio feature table needs contig column"))?;
    let taxonomy_col = find(&["taxonomy", "taxon", "lineage"]);
    let taxonomy_conf_col = find(&["taxonomy_confidence", "tax_confidence", "tax_conf"]);
    let gene_col = find(&["gene_profile", "gene_embedding", "gene"]);
    let gene_conf_col = find(&["gene_confidence", "gene_conf"]);
    let protein_col = find(&["protein_embedding", "esm_embedding", "esmc_embedding", "protein", "esm"]);
    let protein_conf_col = find(&["protein_confidence", "esm_confidence", "protein_conf", "esm_conf"]);

    let mut values = HashMap::new();
    for (line_number, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |index: Option<usize>| -> &str {
            index.and_then(|i| fields.get(i).copied()).unwrap_or("")
        };
        let id = fields.get(contig_col).copied().unwrap_or("").trim();
        if id.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("empty contig id in bio feature row {}", line_number + 2),
            ));
        }
        let taxonomy = nonempty(get(taxonomy_col)).map(str::to_string);
        let taxonomy_confidence = parse_probability(get(taxonomy_conf_col), 0.0, line_number + 2)?;
        let gene_profile = parse_vector(get(gene_col), line_number + 2)?;
        let gene_confidence = parse_probability(
            get(gene_conf_col),
            if gene_profile.is_empty() { 0.0 } else { 1.0 },
            line_number + 2,
        )?;
        let protein_embedding = parse_vector(get(protein_col), line_number + 2)?;
        let protein_confidence = parse_probability(
            get(protein_conf_col),
            if protein_embedding.is_empty() { 0.0 } else { 1.0 },
            line_number + 2,
        )?;
        let feature = BioFeature {
            taxonomy,
            taxonomy_confidence,
            gene_profile,
            gene_confidence,
            protein_embedding,
            protein_confidence,
        };
        if values.insert(id.to_string(), feature).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate contig '{}' in bio feature table", id),
            ));
        }
    }
    Ok(BioFeatureTable { values })
}

pub fn read_link_table<P: AsRef<Path>>(path: P) -> io::Result<LinkTable> {
    let reader = BufReader::new(File::open(path)?);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty link table"))?;
    let columns: Vec<&str> = header.trim().split('\t').collect();
    let find = |names: &[&str]| -> Option<usize> {
        columns
            .iter()
            .position(|column| names.iter().any(|name| column.eq_ignore_ascii_case(name)))
    };
    let left_col = find(&["left", "source", "contig_a", "contig1"])
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "link table needs left/source column"))?;
    let right_col = find(&["right", "target", "contig_b", "contig2"])
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "link table needs right/target column"))?;
    let must_col = find(&["must_link", "positive", "support"]);
    let cannot_col = find(&["cannot_link", "negative", "conflict"]);
    let source_col = find(&["evidence_source", "link_source", "kind"]);

    let mut values = HashMap::new();
    for (line_number, line) in lines.enumerate() {
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
        let must_link = parse_probability(
            must_col.and_then(|i| fields.get(i).copied()).unwrap_or(""),
            0.0,
            line_number + 2,
        )?;
        let cannot_link = parse_probability(
            cannot_col.and_then(|i| fields.get(i).copied()).unwrap_or(""),
            0.0,
            line_number + 2,
        )?;
        let source = source_col
            .and_then(|i| fields.get(i).copied())
            .unwrap_or("")
            .trim()
            .to_string();
        let key = ordered_pair(left, right);
        let entry = values.entry(key).or_insert_with(LinkEvidence::default);
        entry.must_link = entry.must_link.max(must_link);
        entry.cannot_link = entry.cannot_link.max(cannot_link);
        if !source.is_empty() {
            if !entry.source.is_empty() {
                entry.source.push(',');
            }
            entry.source.push_str(&source);
        }
    }
    Ok(LinkTable { values })
}

pub fn bin_contigs_v2(
    contigs: &[Contig],
    coverage: Option<&CoverageTable>,
    markers: Option<&MarkerTable>,
    bio: Option<&BioFeatureTable>,
    links: Option<&LinkTable>,
    cfg: &BridgeBinV2Config,
) -> (BinningResult, BridgeBinV2Stats) {
    let mut stats = BridgeBinV2Stats::default();
    let features: Vec<Feature> = contigs
        .iter()
        .enumerate()
        .filter(|(_, contig)| contig.seq.len() >= cfg.min_contig_len)
        .map(|(contig_index, contig)| feature(contig_index, contig, coverage, markers, bio))
        .collect();
    stats.eligible_contigs = features.len();

    if features.is_empty() {
        return (
            empty_result(contigs),
            BridgeBinV2Stats {
                unbinned_contigs: contigs.len(),
                ..stats
            },
        );
    }

    let feature_by_id: HashMap<&str, usize> = features
        .iter()
        .enumerate()
        .map(|(index, f)| (f.id.as_str(), index))
        .collect();
    let explicit_links = indexed_links(links, &feature_by_id);
    let hard_cannot = hard_cannot_neighbors(features.len(), &explicit_links);
    let candidates = candidate_pairs(&features, &explicit_links, cfg.max_neighbors);
    stats.candidate_edges = candidates.len();

    let mut edges: Vec<CandidateEdge> = candidates
        .into_iter()
        .map(|(left, right)| CandidateEdge {
            left,
            right,
            evidence: pair_evidence(
                &features[left],
                &features[right],
                explicit_links.get(&(left.min(right), left.max(right))),
                cfg,
            ),
        })
        .collect();
    for edge in &edges {
        if edge.evidence.marker_conflict {
            stats.marker_blocked_edges += 1;
        }
        if edge.evidence.taxonomy_conflict {
            stats.taxonomy_blocked_edges += 1;
        }
        if edge.evidence.external_cannot_link {
            stats.external_blocked_edges += 1;
        }
    }
    edges.sort_by(|a, b| {
        edge_priority(&b.evidence)
            .partial_cmp(&edge_priority(&a.evidence))
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.left.cmp(&b.left))
            .then_with(|| a.right.cmp(&b.right))
    });

    let mut components: Vec<Option<Component>> = features
        .iter()
        .map(|feature| Some(Component::from_feature(feature)))
        .collect();
    let mut component_of: Vec<usize> = (0..features.len()).collect();

    for edge in edges {
        if edge.evidence.hard_conflict()
            || edge.evidence.attraction < cfg.core_min_attraction
            || edge.evidence.repulsion > cfg.core_max_repulsion
        {
            continue;
        }
        let left = component_of[edge.left];
        let right = component_of[edge.right];
        if left == right {
            continue;
        }
        if !components_compatible(
            left,
            right,
            &components,
            &component_of,
            &hard_cannot,
            cfg,
        ) {
            continue;
        }
        merge_component_pair(left, right, &mut components, &mut component_of);
        stats.accepted_core_edges += 1;
    }

    let mut core_roots: Vec<usize> = components
        .iter()
        .enumerate()
        .filter_map(|(root, component)| {
            component
                .as_ref()
                .filter(|component| component.bp >= cfg.min_component_bp)
                .map(|_| root)
        })
        .collect();
    core_roots.sort_by_key(|root| std::cmp::Reverse(components[*root].as_ref().unwrap().bp));
    stats.core_bins = core_roots.len();

    let mut small_roots: Vec<usize> = components
        .iter()
        .enumerate()
        .filter_map(|(root, component)| {
            component
                .as_ref()
                .filter(|component| component.bp < cfg.min_component_bp)
                .map(|_| root)
        })
        .collect();
    small_roots.sort_by_key(|root| std::cmp::Reverse(components[*root].as_ref().unwrap().bp));

    for root in small_roots {
        if components[root].is_none() || core_roots.is_empty() {
            continue;
        }
        let mut scored = Vec::new();
        for &core in &core_roots {
            if components[core].is_none()
                || !components_compatible(
                    root,
                    core,
                    &components,
                    &component_of,
                    &hard_cannot,
                    cfg,
                )
            {
                continue;
            }
            let evidence = component_evidence(
                components[root].as_ref().unwrap(),
                components[core].as_ref().unwrap(),
                cfg,
            );
            if !evidence.hard_conflict() && evidence.repulsion <= cfg.core_max_repulsion {
                scored.push((core, evidence.attraction));
            }
        }
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        let Some(&(best, best_score)) = scored.first() else {
            continue;
        };
        let second = scored.get(1).map(|x| x.1).unwrap_or(0.0);
        if best_score >= cfg.rescue_min_attraction && best_score - second >= cfg.rescue_margin {
            merge_component_pair(best, root, &mut components, &mut component_of);
            stats.rescued_components += 1;
        }
    }

    let mut retained_roots: Vec<usize> = core_roots
        .into_iter()
        .filter(|root| components[*root].is_some())
        .collect();
    retained_roots.sort_by_key(|root| std::cmp::Reverse(components[*root].as_ref().unwrap().bp));
    let root_to_bin: HashMap<usize, usize> = retained_roots
        .iter()
        .enumerate()
        .map(|(bin, root)| (*root, bin))
        .collect();

    let eligible_by_contig: HashMap<usize, usize> = features
        .iter()
        .enumerate()
        .map(|(feature_index, f)| (f.contig_index, feature_index))
        .collect();
    let assignments: Vec<Assignment> = contigs
        .iter()
        .enumerate()
        .map(|(contig_index, contig)| {
            let bin_index = eligible_by_contig.get(&contig_index).and_then(|feature_index| {
                let root = component_of[*feature_index];
                root_to_bin.get(&root).copied()
            });
            Assignment {
                contig_id: contig.id.clone(),
                bin_index,
                score: if bin_index.is_some() { 1.0 } else { 0.0 },
                length: contig.seq.len(),
            }
        })
        .collect();
    stats.unbinned_contigs = assignments.iter().filter(|a| a.bin_index.is_none()).count();

    let bins = retained_roots
        .iter()
        .enumerate()
        .map(|(bin_index, root)| {
            let component = components[*root].as_ref().unwrap();
            BinSummary {
                bin_index,
                contig_count: component.members.len(),
                total_bp: component.bp,
                mean_gc: component.mean_gc(),
            }
        })
        .collect();

    (BinningResult { assignments, bins }, stats)
}

fn empty_result(contigs: &[Contig]) -> BinningResult {
    BinningResult {
        assignments: contigs
            .iter()
            .map(|contig| Assignment {
                contig_id: contig.id.clone(),
                bin_index: None,
                score: 0.0,
                length: contig.seq.len(),
            })
            .collect(),
        bins: Vec::new(),
    }
}

fn feature(
    contig_index: usize,
    contig: &Contig,
    coverage: Option<&CoverageTable>,
    markers: Option<&MarkerTable>,
    bio: Option<&BioFeatureTable>,
) -> Feature {
    Feature {
        contig_index,
        id: contig.id.clone(),
        len: contig.seq.len(),
        gc: gc_fraction(&contig.seq),
        kmer: canonical_5mer_frequency(&contig.seq),
        coverage: coverage
            .and_then(|table| table.values.get(&contig.id))
            .cloned()
            .unwrap_or_default(),
        markers: markers
            .and_then(|table| table.values.get(&contig.id))
            .cloned()
            .unwrap_or_default(),
        bio: bio
            .and_then(|table| table.values.get(&contig.id))
            .cloned()
            .unwrap_or_default(),
    }
}

fn candidate_pairs(
    features: &[Feature],
    explicit_links: &HashMap<(usize, usize), LinkEvidence>,
    max_neighbors: usize,
) -> HashSet<(usize, usize)> {
    let mut candidates = HashSet::new();
    let neighbors = max_neighbors.max(1).min(features.len().saturating_sub(1));
    for left in 0..features.len() {
        let mut scored = Vec::new();
        for right in 0..features.len() {
            if left == right {
                continue;
            }
            let gc_delta = (features[left].gc - features[right].gc).abs();
            if gc_delta > 0.14 {
                continue;
            }
            let gc_score = (-gc_delta / 0.055).exp();
            let coverage = coverage_similarity(&features[left].coverage, &features[right].coverage);
            if coverage.is_some_and(|score| score < 0.12) {
                continue;
            }
            let score = coverage.map(|value| 0.75 * value + 0.25 * gc_score).unwrap_or(gc_score);
            scored.push((score, right));
        }
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
        for (_, right) in scored.into_iter().take(neighbors) {
            candidates.insert((left.min(right), left.max(right)));
        }
    }
    candidates.extend(explicit_links.keys().copied());
    candidates
}

fn pair_evidence(
    left: &Feature,
    right: &Feature,
    link: Option<&LinkEvidence>,
    cfg: &BridgeBinV2Config,
) -> PairEvidence {
    let composition = (-hellinger(&left.kmer, &right.kmer) / 0.24).exp();
    let gc_delta = (left.gc - right.gc).abs();
    let gc = (-gc_delta / 0.055).exp();
    let coverage = coverage_similarity(&left.coverage, &right.coverage);
    let gene = vector_similarity(&left.bio.gene_profile, &right.bio.gene_profile);
    let protein = vector_similarity(&left.bio.protein_embedding, &right.bio.protein_embedding);
    let marker_conflict = !left.markers.is_disjoint(&right.markers);
    let taxonomy_conflict = taxonomy_conflict(&left.bio, &right.bio, cfg.taxonomy_confidence);
    let external_cannot_link = link.is_some_and(|value| value.cannot_link >= 0.95);

    let mut weighted = 0.42 * composition + 0.05 * gc;
    let mut weight = 0.47;
    if let Some(value) = coverage {
        let coverage_weight = if left.coverage.len() >= 3 { 0.48 } else { 0.34 };
        weighted += coverage_weight * value;
        weight += coverage_weight;
    }
    if let Some(value) = gene {
        let confidence = left.bio.gene_confidence.min(right.bio.gene_confidence);
        let w = 0.18 * confidence;
        weighted += w * value;
        weight += w;
    }
    if let Some(value) = protein {
        let confidence = left
            .bio
            .protein_confidence
            .min(right.bio.protein_confidence);
        let w = 0.10 * confidence;
        weighted += w * value;
        weight += w;
    }
    if let Some(value) = taxonomy_similarity(&left.bio, &right.bio) {
        let confidence = left
            .bio
            .taxonomy_confidence
            .min(right.bio.taxonomy_confidence);
        let w = 0.14 * confidence;
        weighted += w * value;
        weight += w;
    }
    let mut attraction = if weight > 0.0 { weighted / weight } else { 0.0 };
    if let Some(link) = link {
        attraction = 1.0 - (1.0 - attraction) * (1.0 - 0.75 * link.must_link);
    }

    let mut repulsion: f64 = 0.0;
    if gc_delta > cfg.max_gc_delta {
        repulsion = repulsion.max(((gc_delta - cfg.max_gc_delta) / 0.08).clamp(0.0, 1.0));
    }
    if composition < 0.45 {
        repulsion = repulsion.max(((0.45 - composition) / 0.45).clamp(0.0, 0.85));
    }
    if let Some(value) = coverage {
        if value < 0.38 {
            repulsion = repulsion.max(((0.38 - value) / 0.38).clamp(0.0, 0.90));
        }
    }
    if let Some(value) = gene {
        let confidence = left.bio.gene_confidence.min(right.bio.gene_confidence);
        if value < 0.30 && confidence >= 0.7 {
            repulsion = repulsion.max((0.30 - value) / 0.30 * confidence * 0.8);
        }
    }
    if let Some(value) = protein {
        let confidence = left
            .bio
            .protein_confidence
            .min(right.bio.protein_confidence);
        if value < 0.20 && confidence >= 0.7 {
            repulsion = repulsion.max((0.20 - value) / 0.20 * confidence * 0.5);
        }
    }
    if marker_conflict && cfg.hard_marker_veto {
        repulsion = 1.0;
    }
    if taxonomy_conflict {
        repulsion = 1.0;
    }
    if let Some(link) = link {
        repulsion = repulsion.max(link.cannot_link);
    }

    PairEvidence {
        attraction: attraction.clamp(0.0, 1.0),
        repulsion: repulsion.clamp(0.0, 1.0),
        composition,
        coverage,
        gc,
        gene,
        protein,
        marker_conflict: marker_conflict && cfg.hard_marker_veto,
        taxonomy_conflict,
        external_cannot_link,
    }
}

fn edge_priority(evidence: &PairEvidence) -> f64 {
    evidence.attraction - 1.4 * evidence.repulsion
}

impl Component {
    fn from_feature(feature: &Feature) -> Self {
        let weight = feature.len as f64;
        let mut kmer_sum = [0.0; 1024];
        for (target, value) in kmer_sum.iter_mut().zip(feature.kmer.iter()) {
            *target = *value * weight;
        }
        let coverage_sum = feature.coverage.iter().map(|value| value * weight).collect();
        let mut component = Self {
            members: vec![feature.contig_index],
            bp: feature.len,
            gc_sum: feature.gc * weight,
            kmer_sum,
            coverage_sum,
            markers: feature.markers.clone(),
            taxonomies: feature
                .bio
                .taxonomy
                .as_ref()
                .map(|taxonomy| vec![(taxonomy.clone(), feature.bio.taxonomy_confidence)])
                .unwrap_or_default(),
            gene_sum: Vec::new(),
            gene_weight: 0.0,
            protein_sum: Vec::new(),
            protein_weight: 0.0,
        };
        component.add_bio_vector(
            &feature.bio.gene_profile,
            feature.bio.gene_confidence * weight,
            true,
        );
        component.add_bio_vector(
            &feature.bio.protein_embedding,
            feature.bio.protein_confidence * weight,
            false,
        );
        component
    }

    fn add_bio_vector(&mut self, values: &[f64], weight: f64, gene: bool) {
        if values.is_empty() || weight <= 0.0 {
            return;
        }
        let (sum, total_weight) = if gene {
            (&mut self.gene_sum, &mut self.gene_weight)
        } else {
            (&mut self.protein_sum, &mut self.protein_weight)
        };
        if sum.is_empty() {
            sum.resize(values.len(), 0.0);
        }
        if sum.len() != values.len() {
            return;
        }
        for (target, value) in sum.iter_mut().zip(values.iter()) {
            *target += *value * weight;
        }
        *total_weight += weight;
    }

    fn absorb(&mut self, mut other: Component) {
        self.members.append(&mut other.members);
        self.bp += other.bp;
        self.gc_sum += other.gc_sum;
        for (target, value) in self.kmer_sum.iter_mut().zip(other.kmer_sum.iter()) {
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
        self.taxonomies.append(&mut other.taxonomies);
        merge_vector_sum(
            &mut self.gene_sum,
            &mut self.gene_weight,
            &other.gene_sum,
            other.gene_weight,
        );
        merge_vector_sum(
            &mut self.protein_sum,
            &mut self.protein_weight,
            &other.protein_sum,
            other.protein_weight,
        );
    }

    fn mean_gc(&self) -> f64 {
        self.gc_sum / self.bp.max(1) as f64
    }

    fn kmer_frequency(&self) -> [f64; 1024] {
        let mut out = [0.0; 1024];
        let bp = self.bp.max(1) as f64;
        for (target, value) in out.iter_mut().zip(self.kmer_sum.iter()) {
            *target = *value / bp;
        }
        renormalize(&mut out);
        out
    }

    fn coverage_mean(&self) -> Vec<f64> {
        let bp = self.bp.max(1) as f64;
        self.coverage_sum.iter().map(|value| value / bp).collect()
    }

    fn gene_mean(&self) -> Vec<f64> {
        if self.gene_weight <= 0.0 {
            Vec::new()
        } else {
            self.gene_sum
                .iter()
                .map(|value| value / self.gene_weight)
                .collect()
        }
    }

    fn protein_mean(&self) -> Vec<f64> {
        if self.protein_weight <= 0.0 {
            Vec::new()
        } else {
            self.protein_sum
                .iter()
                .map(|value| value / self.protein_weight)
                .collect()
        }
    }
}

fn merge_vector_sum(
    target: &mut Vec<f64>,
    target_weight: &mut f64,
    source: &[f64],
    source_weight: f64,
) {
    if source.is_empty() || source_weight <= 0.0 {
        return;
    }
    if target.is_empty() {
        target.resize(source.len(), 0.0);
    }
    if target.len() != source.len() {
        return;
    }
    for (left, right) in target.iter_mut().zip(source.iter()) {
        *left += *right;
    }
    *target_weight += source_weight;
}

fn component_evidence(left: &Component, right: &Component, cfg: &BridgeBinV2Config) -> PairEvidence {
    let left_kmer = left.kmer_frequency();
    let right_kmer = right.kmer_frequency();
    let composition = (-hellinger(&left_kmer, &right_kmer) / 0.24).exp();
    let gc_delta = (left.mean_gc() - right.mean_gc()).abs();
    let gc = (-gc_delta / 0.055).exp();
    let coverage = coverage_similarity(&left.coverage_mean(), &right.coverage_mean());
    let gene = vector_similarity(&left.gene_mean(), &right.gene_mean());
    let protein = vector_similarity(&left.protein_mean(), &right.protein_mean());
    let marker_conflict = !left.markers.is_disjoint(&right.markers) && cfg.hard_marker_veto;
    let taxonomy_conflict = component_taxonomy_conflict(left, right, cfg.taxonomy_confidence);

    let mut weighted = 0.45 * composition + 0.05 * gc;
    let mut weight = 0.50;
    if let Some(value) = coverage {
        weighted += 0.45 * value;
        weight += 0.45;
    }
    if let Some(value) = gene {
        weighted += 0.16 * value;
        weight += 0.16;
    }
    if let Some(value) = protein {
        weighted += 0.08 * value;
        weight += 0.08;
    }
    let attraction = weighted / weight;
    let mut repulsion: f64 = 0.0;
    if gc_delta > cfg.max_gc_delta {
        repulsion = repulsion.max(((gc_delta - cfg.max_gc_delta) / 0.08).clamp(0.0, 1.0));
    }
    if composition < cfg.min_component_composition {
        repulsion = repulsion.max(
            ((cfg.min_component_composition - composition) / cfg.min_component_composition)
                .clamp(0.0, 0.9),
        );
    }
    if let Some(value) = coverage {
        if value < cfg.min_component_coverage {
            repulsion = repulsion.max(
                ((cfg.min_component_coverage - value) / cfg.min_component_coverage).clamp(0.0, 0.9),
            );
        }
    }
    if marker_conflict || taxonomy_conflict {
        repulsion = 1.0;
    }
    PairEvidence {
        attraction,
        repulsion,
        composition,
        coverage,
        gc,
        gene,
        protein,
        marker_conflict,
        taxonomy_conflict,
        external_cannot_link: false,
    }
}

fn components_compatible(
    left_root: usize,
    right_root: usize,
    components: &[Option<Component>],
    component_of: &[usize],
    hard_cannot: &[HashSet<usize>],
    cfg: &BridgeBinV2Config,
) -> bool {
    let Some(left) = components[left_root].as_ref() else {
        return false;
    };
    let Some(right) = components[right_root].as_ref() else {
        return false;
    };
    if cfg.hard_marker_veto && !left.markers.is_disjoint(&right.markers) {
        return false;
    }
    if component_taxonomy_conflict(left, right, cfg.taxonomy_confidence) {
        return false;
    }
    let (smaller, other_root) = if left.members.len() <= right.members.len() {
        (&left.members, right_root)
    } else {
        (&right.members, left_root)
    };
    for &contig_index in smaller {
        for &blocked in &hard_cannot[contig_index] {
            if component_of[blocked] == other_root {
                return false;
            }
        }
    }
    let evidence = component_evidence(left, right, cfg);
    evidence.attraction >= cfg.component_min_attraction
        && evidence.repulsion <= cfg.core_max_repulsion
        && evidence.composition >= cfg.min_component_composition
        && evidence
            .coverage
            .map(|value| value >= cfg.min_component_coverage)
            .unwrap_or(true)
        && (left.mean_gc() - right.mean_gc()).abs() <= cfg.max_gc_delta
}

fn merge_component_pair(
    left_root: usize,
    right_root: usize,
    components: &mut [Option<Component>],
    component_of: &mut [usize],
) -> usize {
    if left_root == right_root {
        return left_root;
    }
    let left_bp = components[left_root].as_ref().map(|c| c.bp).unwrap_or(0);
    let right_bp = components[right_root].as_ref().map(|c| c.bp).unwrap_or(0);
    let (keep, remove) = if left_bp >= right_bp {
        (left_root, right_root)
    } else {
        (right_root, left_root)
    };
    let other = components[remove].take().unwrap();
    let moved = other.members.clone();
    components[keep].as_mut().unwrap().absorb(other);
    for member in moved {
        component_of[member] = keep;
    }
    keep
}

fn indexed_links(
    links: Option<&LinkTable>,
    feature_by_id: &HashMap<&str, usize>,
) -> HashMap<(usize, usize), LinkEvidence> {
    let mut out = HashMap::new();
    let Some(links) = links else {
        return out;
    };
    for ((left, right), evidence) in &links.values {
        let (Some(&a), Some(&b)) = (
            feature_by_id.get(left.as_str()),
            feature_by_id.get(right.as_str()),
        ) else {
            continue;
        };
        out.insert((a.min(b), a.max(b)), evidence.clone());
    }
    out
}

fn hard_cannot_neighbors(
    n: usize,
    links: &HashMap<(usize, usize), LinkEvidence>,
) -> Vec<HashSet<usize>> {
    let mut out = vec![HashSet::new(); n];
    for (&(left, right), evidence) in links {
        if evidence.cannot_link >= 0.95 {
            out[left].insert(right);
            out[right].insert(left);
        }
    }
    out
}

fn taxonomy_conflict(left: &BioFeature, right: &BioFeature, min_confidence: f64) -> bool {
    if left.taxonomy_confidence < min_confidence || right.taxonomy_confidence < min_confidence {
        return false;
    }
    match (&left.taxonomy, &right.taxonomy) {
        (Some(a), Some(b)) => lineages_conflict(a, b),
        _ => false,
    }
}

fn component_taxonomy_conflict(left: &Component, right: &Component, min_confidence: f64) -> bool {
    left.taxonomies.iter().any(|(a, a_conf)| {
        *a_conf >= min_confidence
            && right.taxonomies.iter().any(|(b, b_conf)| {
                *b_conf >= min_confidence && lineages_conflict(a, b)
            })
    })
}

fn taxonomy_similarity(left: &BioFeature, right: &BioFeature) -> Option<f64> {
    let (Some(a), Some(b)) = (&left.taxonomy, &right.taxonomy) else {
        return None;
    };
    let a: Vec<&str> = lineage_parts(a);
    let b: Vec<&str> = lineage_parts(b);
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let shared = a.iter().zip(b.iter()).take_while(|(x, y)| x == y).count();
    Some(shared as f64 / a.len().min(b.len()) as f64)
}

fn lineages_conflict(left: &str, right: &str) -> bool {
    let left = lineage_parts(left);
    let right = lineage_parts(right);
    if left.is_empty() || right.is_empty() {
        return false;
    }
    left.iter()
        .zip(right.iter())
        .any(|(a, b)| !a.eq_ignore_ascii_case(b))
}

fn lineage_parts(lineage: &str) -> Vec<&str> {
    lineage
        .split([';', '|'])
        .map(str::trim)
        .filter(|value| !value.is_empty() && *value != "." && *value != "NA")
        .collect()
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
    match correlation {
        Some(correlation) => Some((0.68 * ratio_score + 0.32 * ((correlation + 1.0) * 0.5)).clamp(0.0, 1.0)),
        None => Some(ratio_score),
    }
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

fn vector_similarity(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.is_empty() || left.len() != right.len() {
        return None;
    }
    let mut dot = 0.0;
    let mut left_norm = 0.0;
    let mut right_norm = 0.0;
    for (&a, &b) in left.iter().zip(right.iter()) {
        dot += a * b;
        left_norm += a * a;
        right_norm += b * b;
    }
    let denominator = (left_norm * right_norm).sqrt();
    if denominator <= 1e-12 {
        None
    } else {
        Some(((dot / denominator).clamp(-1.0, 1.0) + 1.0) * 0.5)
    }
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

fn renormalize(values: &mut [f64; 1024]) {
    let total = values.iter().sum::<f64>();
    if total > 0.0 {
        for value in values {
            *value /= total;
        }
    }
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

fn hellinger(left: &[f64; 1024], right: &[f64; 1024]) -> f64 {
    (0.5
        * left
            .iter()
            .zip(right.iter())
            .map(|(a, b)| (a.sqrt() - b.sqrt()).powi(2))
            .sum::<f64>())
    .sqrt()
}

fn parse_vector(raw: &str, line_number: usize) -> io::Result<Vec<f64>> {
    let raw = raw.trim();
    if raw.is_empty() || raw == "." || raw.eq_ignore_ascii_case("NA") {
        return Ok(Vec::new());
    }
    raw.split(',')
        .map(|value| {
            let parsed = value.trim().parse::<f64>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid vector value '{}' at row {}", value, line_number),
                )
            })?;
            if !parsed.is_finite() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("non-finite vector value at row {}", line_number),
                ));
            }
            Ok(parsed)
        })
        .collect()
}

fn parse_probability(raw: &str, default: f64, line_number: usize) -> io::Result<f64> {
    let raw = raw.trim();
    if raw.is_empty() || raw == "." || raw.eq_ignore_ascii_case("NA") {
        return Ok(default);
    }
    let value = raw.parse::<f64>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid probability '{}' at row {}", raw, line_number),
        )
    })?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("probability out of range at row {}", line_number),
        ));
    }
    Ok(value)
}

fn nonempty(value: &str) -> Option<&str> {
    let value = value.trim();
    (!value.is_empty() && value != "." && !value.eq_ignore_ascii_case("NA")).then_some(value)
}

fn ordered_pair(left: &str, right: &str) -> (String, String) {
    if left <= right {
        (left.to_string(), right.to_string())
    } else {
        (right.to_string(), left.to_string())
    }
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

    fn tiny_config() -> BridgeBinV2Config {
        BridgeBinV2Config {
            min_contig_len: 100,
            min_component_bp: 100,
            max_neighbors: 8,
            core_min_attraction: 0.70,
            component_min_attraction: 0.68,
            ..Default::default()
        }
    }

    #[test]
    fn duplicate_marker_is_a_hard_cannot_link() {
        let contigs = vec![
            contig("a", "AAAACCCCGGGGTTTT", 40),
            contig("b", "AAAACCCCGGGGTTTT", 40),
        ];
        let markers = MarkerTable {
            values: HashMap::from([
                ("a".to_string(), HashSet::from(["m1".to_string()])),
                ("b".to_string(), HashSet::from(["m1".to_string()])),
            ]),
        };
        let (result, _) = bin_contigs_v2(
            &contigs,
            None,
            Some(&markers),
            None,
            None,
            &tiny_config(),
        );
        assert_ne!(result.assignments[0].bin_index, result.assignments[1].bin_index);
    }

    #[test]
    fn external_cannot_link_survives_transitive_merges() {
        let contigs = vec![
            contig("a", "AAAACCCCGGGGTTTT", 40),
            contig("b", "AAAACCCCGGGGTTTT", 40),
            contig("c", "AAAACCCCGGGGTTTT", 40),
        ];
        let links = LinkTable {
            values: HashMap::from([(
                ordered_pair("a", "c"),
                LinkEvidence {
                    must_link: 0.0,
                    cannot_link: 1.0,
                    source: "test".to_string(),
                },
            )]),
        };
        let (result, _) = bin_contigs_v2(
            &contigs,
            None,
            None,
            None,
            Some(&links),
            &tiny_config(),
        );
        assert_ne!(result.assignments[0].bin_index, result.assignments[2].bin_index);
    }

    #[test]
    fn high_confidence_taxonomy_is_a_veto_but_prefix_is_compatible() {
        let a = BioFeature {
            taxonomy: Some("Bacteria;Proteobacteria;Enterobacterales;Escherichia".to_string()),
            taxonomy_confidence: 0.99,
            ..Default::default()
        };
        let b = BioFeature {
            taxonomy: Some("Bacteria;Proteobacteria;Enterobacterales;Escherichia;coli".to_string()),
            taxonomy_confidence: 0.99,
            ..Default::default()
        };
        let c = BioFeature {
            taxonomy: Some("Bacteria;Proteobacteria;Enterobacterales;Salmonella".to_string()),
            taxonomy_confidence: 0.99,
            ..Default::default()
        };
        assert!(!taxonomy_conflict(&a, &b, 0.90));
        assert!(taxonomy_conflict(&a, &c, 0.90));
    }

    #[test]
    fn gene_profiles_can_resolve_coverage_matched_groups() {
        let contigs = vec![
            contig("a1", "AAAACCCCGGGGTTTT", 40),
            contig("a2", "AAAACCCCGGGGTTTT", 40),
            contig("b1", "AAAACCCCGGGGTTTT", 40),
            contig("b2", "AAAACCCCGGGGTTTT", 40),
        ];
        let coverage = CoverageTable {
            sample_names: vec!["s1".to_string(), "s2".to_string(), "s3".to_string()],
            values: contigs
                .iter()
                .map(|contig| (contig.id.clone(), vec![20.0, 8.0, 30.0]))
                .collect(),
        };
        let bio = BioFeatureTable {
            values: HashMap::from([
                ("a1".to_string(), BioFeature { gene_profile: vec![1.0, 0.0], gene_confidence: 1.0, ..Default::default() }),
                ("a2".to_string(), BioFeature { gene_profile: vec![1.0, 0.0], gene_confidence: 1.0, ..Default::default() }),
                ("b1".to_string(), BioFeature { gene_profile: vec![0.0, 1.0], gene_confidence: 1.0, ..Default::default() }),
                ("b2".to_string(), BioFeature { gene_profile: vec![0.0, 1.0], gene_confidence: 1.0, ..Default::default() }),
            ]),
        };
        let mut cfg = tiny_config();
        cfg.core_min_attraction = 0.76;
        let (result, _) = bin_contigs_v2(
            &contigs,
            Some(&coverage),
            None,
            Some(&bio),
            None,
            &cfg,
        );
        let bins: HashMap<&str, Option<usize>> = result
            .assignments
            .iter()
            .map(|assignment| (assignment.contig_id.as_str(), assignment.bin_index))
            .collect();
        assert_eq!(bins["a1"], bins["a2"]);
        assert_eq!(bins["b1"], bins["b2"]);
        assert_ne!(bins["a1"], bins["b1"]);
    }
}
