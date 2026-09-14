use crate::dna::{base_bits, canonical_kmers, KmerKey, OrientedKmer, MAX_K};
use crate::fastq::for_each_pair;
use crate::graph::{compact_unitigs, summarize, GraphSummary, RawGraph, UnitigGraph};
use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct MultiKConfig {
    pub read1: PathBuf,
    pub read2: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub ks: Vec<usize>,
    pub min_count: u32,
    pub min_fragment_support: u32,
    pub min_mean_quality: f32,
    pub max_pairs: Option<usize>,
    pub max_rescue_bases: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct ObservationEvidence {
    count: u32,
    fragment_count: u32,
    quality_sum: u64,
}

impl ObservationEvidence {
    fn mean_quality(self, k: usize) -> f32 {
        if self.count == 0 || k == 0 {
            0.0
        } else {
            self.quality_sum as f32 / (self.count as f32 * k as f32)
        }
    }
}

#[derive(Debug)]
struct PendingLayer {
    k: usize,
    node_evidence: FxHashMap<KmerKey, ObservationEvidence>,
    edge_evidence: FxHashMap<KmerKey, ObservationEvidence>,
    observations: u64,
    edge_observations: u64,
}

impl PendingLayer {
    fn new(k: usize) -> Self {
        Self {
            k,
            node_evidence: FxHashMap::default(),
            edge_evidence: FxHashMap::default(),
            observations: 0,
            edge_observations: 0,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct MultiKLayerSummary {
    pub k: usize,
    pub observations: u64,
    pub edge_observations: u64,
    pub distinct_kmers: usize,
    pub retained_kmers: usize,
    pub distinct_kplus1: usize,
    pub retained_directed_edges: usize,
    pub graph: GraphSummary,
}

#[derive(Debug)]
pub struct MultiKLayer {
    pub k: usize,
    pub raw_graph: RawGraph,
    pub unitig_graph: UnitigGraph,
    pub summary: MultiKLayerSummary,
}

#[derive(Clone, Debug)]
pub struct ProjectionPath {
    pub high_unitig: u32,
    pub low_unitigs: Vec<u32>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ProjectionSummary {
    pub high_k: usize,
    pub low_k: usize,
    pub high_unitigs: usize,
    pub exact_unitigs: usize,
    pub missing_node_unitigs: usize,
    pub missing_edge_unitigs: usize,
    pub projected_low_unitig_visits: usize,
}

#[derive(Debug)]
pub struct ProjectionPair {
    pub high_k: usize,
    pub low_k: usize,
    pub exact_paths: Vec<ProjectionPath>,
    pub summary: ProjectionSummary,
}

#[derive(Clone, Debug, Serialize)]
pub struct AdaptiveRescueSummary {
    pub high_k: usize,
    pub low_k: usize,
    pub dead_end_unitigs: usize,
    pub low_anchor_missing: usize,
    pub branch_stops: usize,
    pub dead_stops: usize,
    pub cycle_stops: usize,
    pub maxed_out: usize,
    pub reanchored: usize,
    pub total_extension_bases: usize,
    pub max_extension_bases: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct AdaptiveRescueCandidate {
    pub high_k: usize,
    pub low_k: usize,
    pub high_unitig: u32,
    pub reanchor_state: u32,
    pub extension_bases: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct MultiKTimingSummary {
    pub count_seconds: f64,
    pub finalize_seconds: f64,
    pub projection_seconds: f64,
    pub rescue_seconds: f64,
    pub total_seconds: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct MultiKSummary {
    pub version: String,
    pub read_pairs: usize,
    pub ks: Vec<usize>,
    pub min_count: u32,
    pub min_fragment_support: u32,
    pub min_mean_quality: f32,
    pub max_rescue_bases: usize,
    pub layers: Vec<MultiKLayerSummary>,
    pub projections: Vec<ProjectionSummary>,
    pub adaptive_rescues: Vec<AdaptiveRescueSummary>,
    pub timings_seconds: MultiKTimingSummary,
}

#[derive(Debug)]
pub struct MultiKGraph {
    pub layers: Vec<MultiKLayer>,
    pub projection_pairs: Vec<ProjectionPair>,
    pub rescue_candidates: Vec<AdaptiveRescueCandidate>,
    pub summary: MultiKSummary,
}

#[derive(Clone, Copy, Debug)]
struct RollerState {
    k: usize,
    forward: KmerKey,
    reverse: KmerKey,
    valid: usize,
}

impl RollerState {
    fn new(k: usize) -> Self {
        Self {
            k,
            forward: KmerKey::ZERO,
            reverse: KmerKey::ZERO,
            valid: 0,
        }
    }

    fn reset(&mut self) {
        self.forward = KmerKey::ZERO;
        self.reverse = KmerKey::ZERO;
        self.valid = 0;
    }
}

fn for_each_multi_canonical_kmer<F>(
    sequence: &[u8],
    orders: &[usize],
    mut callback: F,
) -> Result<usize>
where
    F: FnMut(usize, OrientedKmer),
{
    if orders.is_empty() {
        bail!("at least one k-mer order is required");
    }
    let mut previous = None;
    for &k in orders {
        if k == 0 || k > MAX_K {
            bail!("k-mer order must be in 1..={MAX_K}");
        }
        if previous.is_some_and(|value| k <= value) {
            bail!("multi-k orders must be strictly increasing");
        }
        previous = Some(k);
    }

    let mut states: Vec<RollerState> = orders.iter().copied().map(RollerState::new).collect();
    let mut emitted = 0_usize;

    for (index, &base) in sequence.iter().enumerate() {
        let Some(bits) = base_bits(base) else {
            for state in &mut states {
                state.reset();
            }
            continue;
        };

        for state in &mut states {
            state.forward.shift_left_append(bits, state.k);
            state.valid += 1;
            if state.valid < state.k {
                continue;
            }
            if state.valid == state.k {
                state.reverse = state.forward.reverse_complement(state.k);
            } else {
                state
                    .reverse
                    .shift_right_prepend_complement(bits, state.k);
            }

            let position = index + 1 - state.k;
            if state.reverse < state.forward {
                callback(
                    state.k,
                    OrientedKmer {
                        key: state.reverse,
                        reverse: true,
                        position,
                    },
                );
            } else {
                callback(
                    state.k,
                    OrientedKmer {
                        key: state.forward,
                        reverse: false,
                        position,
                    },
                );
            }
            emitted += 1;
        }
    }
    Ok(emitted)
}

fn scan_record(
    sequence: &[u8],
    quality: &[u8],
    orders: &[usize],
    pending: &mut [PendingLayer],
    node_seen: &mut [FxHashSet<KmerKey>],
    edge_seen: &mut [FxHashSet<KmerKey>],
) -> Result<()> {
    let mut quality_prefix = vec![0_u64; quality.len() + 1];
    for (index, value) in quality.iter().enumerate() {
        quality_prefix[index + 1] = quality_prefix[index] + u64::from(value.saturating_sub(33));
    }

    for_each_multi_canonical_kmer(sequence, orders, |order, item| {
        for (layer_index, layer) in pending.iter_mut().enumerate() {
            if order == layer.k {
                let window_quality =
                    quality_prefix[item.position + order] - quality_prefix[item.position];
                let entry = layer.node_evidence.entry(item.key).or_default();
                entry.count = entry.count.saturating_add(1);
                entry.quality_sum = entry.quality_sum.saturating_add(window_quality);
                layer.observations = layer.observations.saturating_add(1);
                node_seen[layer_index].insert(item.key);
            } else if order == layer.k + 1 {
                let entry = layer.edge_evidence.entry(item.key).or_default();
                entry.count = entry.count.saturating_add(1);
                layer.edge_observations = layer.edge_observations.saturating_add(1);
                edge_seen[layer_index].insert(item.key);
            }
        }
    })?;
    Ok(())
}

fn count_multi_k(
    read1: &Path,
    read2: Option<&Path>,
    ks: &[usize],
    max_pairs: Option<usize>,
) -> Result<(Vec<PendingLayer>, usize)> {
    let mut orders = Vec::with_capacity(ks.len() * 2);
    for &k in ks {
        orders.push(k);
        orders.push(k + 1);
    }
    orders.sort_unstable();
    orders.dedup();

    let mut pending: Vec<PendingLayer> = ks.iter().copied().map(PendingLayer::new).collect();
    let mut node_seen: Vec<FxHashSet<KmerKey>> =
        (0..ks.len()).map(|_| FxHashSet::default()).collect();
    let mut edge_seen: Vec<FxHashSet<KmerKey>> =
        (0..ks.len()).map(|_| FxHashSet::default()).collect();

    let read_pairs = for_each_pair(read1, read2, max_pairs, |_pair_index, left, right| {
        for values in &mut node_seen {
            values.clear();
        }
        for values in &mut edge_seen {
            values.clear();
        }

        scan_record(
            &left.sequence,
            &left.quality,
            &orders,
            &mut pending,
            &mut node_seen,
            &mut edge_seen,
        )?;
        if let Some(right) = right {
            scan_record(
                &right.sequence,
                &right.quality,
                &orders,
                &mut pending,
                &mut node_seen,
                &mut edge_seen,
            )?;
        }

        for layer_index in 0..pending.len() {
            for key in node_seen[layer_index].drain() {
                if let Some(entry) = pending[layer_index].node_evidence.get_mut(&key) {
                    entry.fragment_count = entry.fragment_count.saturating_add(1);
                }
            }
            for key in edge_seen[layer_index].drain() {
                if let Some(entry) = pending[layer_index].edge_evidence.get_mut(&key) {
                    entry.fragment_count = entry.fragment_count.saturating_add(1);
                }
            }
        }
        Ok(())
    })?;

    Ok((pending, read_pairs))
}

fn oriented_state(
    sequence: &[u8],
    k: usize,
    index: &FxHashMap<KmerKey, u32>,
) -> Result<Option<u32>> {
    let key = KmerKey::from_sequence(sequence)?;
    let (canonical, reverse) = key.canonical(k);
    Ok(index
        .get(&canonical)
        .copied()
        .map(|node| node * 2 + u32::from(reverse)))
}

fn finalize_layer(
    pending: PendingLayer,
    min_count: u32,
    min_fragment_support: u32,
    min_mean_quality: f32,
) -> Result<MultiKLayer> {
    let distinct_kmers = pending.node_evidence.len();
    let distinct_kplus1 = pending.edge_evidence.len();
    let mut keys: Vec<KmerKey> = pending
        .node_evidence
        .iter()
        .filter_map(|(key, value)| {
            (value.count >= min_count
                && value.fragment_count >= min_fragment_support
                && value.mean_quality(pending.k) >= min_mean_quality)
                .then_some(*key)
        })
        .collect();
    keys.sort_unstable();
    if keys.len() > (u32::MAX as usize) / 2 {
        bail!("k={} graph has too many canonical nodes", pending.k);
    }

    let index: FxHashMap<KmerKey, u32> = keys
        .iter()
        .enumerate()
        .map(|(node_id, key)| (*key, node_id as u32))
        .collect();
    let counts: Vec<u32> = keys
        .iter()
        .map(|key| pending.node_evidence.get(key).map_or(0, |value| value.count))
        .collect();

    let mut edges = Vec::new();
    let mut singleton_solid_edges = 0_usize;
    for (edge_key, evidence) in &pending.edge_evidence {
        let sequence = edge_key.to_sequence(pending.k + 1);
        let Some(source) = oriented_state(&sequence[..pending.k], pending.k, &index)? else {
            continue;
        };
        let Some(target) = oriented_state(&sequence[1..], pending.k, &index)? else {
            continue;
        };
        edges.push((source, target));
        edges.push((RawGraph::reverse_state(target), RawGraph::reverse_state(source)));
        if evidence.count < min_count {
            singleton_solid_edges += 2;
        }
    }
    edges.sort_unstable();
    edges.dedup();

    let state_count = keys.len() * 2;
    let mut out_offsets = vec![0_u64; state_count + 1];
    let mut indegree = vec![0_u32; state_count];
    for &(source, target) in &edges {
        out_offsets[source as usize + 1] += 1;
        indegree[target as usize] = indegree[target as usize].saturating_add(1);
    }
    for offset in 1..out_offsets.len() {
        out_offsets[offset] += out_offsets[offset - 1];
    }
    let out_targets: Vec<u32> = edges.into_iter().map(|(_, target)| target).collect();
    let raw_graph = RawGraph {
        k: pending.k,
        keys,
        counts,
        out_offsets,
        out_targets,
        indegree,
        singleton_solid_edges,
    };
    let unitig_graph = compact_unitigs(&raw_graph);
    let graph_summary = summarize(&raw_graph, &unitig_graph);
    let summary = MultiKLayerSummary {
        k: pending.k,
        observations: pending.observations,
        edge_observations: pending.edge_observations,
        distinct_kmers,
        retained_kmers: raw_graph.keys.len(),
        distinct_kplus1,
        retained_directed_edges: raw_graph.out_targets.len(),
        graph: graph_summary,
    };
    Ok(MultiKLayer {
        k: pending.k,
        raw_graph,
        unitig_graph,
        summary,
    })
}

fn project_high_unitig_to_low(
    high_sequence: &[u8],
    low: &MultiKLayer,
) -> Result<std::result::Result<Vec<u32>, bool>> {
    let kmers = canonical_kmers(high_sequence, low.k)?;
    let mut states = Vec::with_capacity(kmers.len());
    for item in kmers {
        let Ok(node) = low.raw_graph.keys.binary_search(&item.key) else {
            return Ok(Err(false));
        };
        states.push(node as u32 * 2 + u32::from(item.reverse));
    }

    let mut low_unitigs = Vec::new();
    for pair in states.windows(2) {
        let Some(&unitig) = low.unitig_graph.edge_to_unitig.get(&(pair[0], pair[1])) else {
            return Ok(Err(true));
        };
        if low_unitigs.last().copied() != Some(unitig) {
            low_unitigs.push(unitig);
        }
    }
    Ok(Ok(low_unitigs))
}

fn build_projection_pair(high: &MultiKLayer, low: &MultiKLayer) -> Result<ProjectionPair> {
    let mut exact_paths = Vec::new();
    let mut missing_node_unitigs = 0_usize;
    let mut missing_edge_unitigs = 0_usize;
    let mut projected_low_unitig_visits = 0_usize;

    for unitig in &high.unitig_graph.unitigs {
        match project_high_unitig_to_low(&unitig.sequence, low)? {
            Ok(low_unitigs) => {
                projected_low_unitig_visits += low_unitigs.len();
                exact_paths.push(ProjectionPath {
                    high_unitig: unitig.id,
                    low_unitigs,
                });
            }
            Err(true) => missing_edge_unitigs += 1,
            Err(false) => missing_node_unitigs += 1,
        }
    }

    let summary = ProjectionSummary {
        high_k: high.k,
        low_k: low.k,
        high_unitigs: high.unitig_graph.unitigs.len(),
        exact_unitigs: exact_paths.len(),
        missing_node_unitigs,
        missing_edge_unitigs,
        projected_low_unitig_visits,
    };
    Ok(ProjectionPair {
        high_k: high.k,
        low_k: low.k,
        exact_paths,
        summary,
    })
}

fn probe_adaptive_rescue(
    high: &MultiKLayer,
    low: &MultiKLayer,
    max_bases: usize,
) -> Result<(AdaptiveRescueSummary, Vec<AdaptiveRescueCandidate>)> {
    let mut summary = AdaptiveRescueSummary {
        high_k: high.k,
        low_k: low.k,
        dead_end_unitigs: 0,
        low_anchor_missing: 0,
        branch_stops: 0,
        dead_stops: 0,
        cycle_stops: 0,
        maxed_out: 0,
        reanchored: 0,
        total_extension_bases: 0,
        max_extension_bases: 0,
    };
    let mut candidates = Vec::new();

    for unitig in &high.unitig_graph.unitigs {
        if high.unitig_graph.outdegree(unitig.id) != 0 || unitig.sequence.len() < high.k {
            continue;
        }
        summary.dead_end_unitigs += 1;

        let low_kmers = canonical_kmers(&unitig.sequence, low.k)?;
        let Some(anchor) = low_kmers.last().copied() else {
            summary.low_anchor_missing += 1;
            continue;
        };
        let Ok(low_node) = low.raw_graph.keys.binary_search(&anchor.key) else {
            summary.low_anchor_missing += 1;
            continue;
        };
        let mut current = low_node as u32 * 2 + u32::from(anchor.reverse);
        let mut visited = FxHashSet::default();
        visited.insert(current);
        let own_states: FxHashSet<u32> = unitig.states.iter().copied().collect();
        let context_start = unitig.sequence.len().saturating_sub(high.k - 1);
        let mut context = unitig.sequence[context_start..].to_vec();
        let mut extension_bases = 0_usize;
        let mut terminal_reason = None;

        while extension_bases < max_bases {
            let range = low.raw_graph.out_range(current);
            if range.is_empty() {
                terminal_reason = Some(0_u8);
                break;
            }
            if range.len() != 1 {
                terminal_reason = Some(1_u8);
                break;
            }
            let next = low.raw_graph.out_targets[range.start];
            if !visited.insert(next) {
                terminal_reason = Some(2_u8);
                break;
            }
            let sequence = low.raw_graph.state_sequence(next);
            let Some(&base) = sequence.last() else {
                terminal_reason = Some(0_u8);
                break;
            };
            context.push(base);
            extension_bases += 1;
            current = next;

            if context.len() >= high.k {
                let start = context.len() - high.k;
                let key = KmerKey::from_sequence(&context[start..])?;
                let (canonical, reverse) = key.canonical(high.k);
                if let Ok(node) = high.raw_graph.keys.binary_search(&canonical) {
                    let high_state = node as u32 * 2 + u32::from(reverse);
                    if !own_states.contains(&high_state) {
                        summary.reanchored += 1;
                        summary.total_extension_bases += extension_bases;
                        summary.max_extension_bases =
                            summary.max_extension_bases.max(extension_bases);
                        candidates.push(AdaptiveRescueCandidate {
                            high_k: high.k,
                            low_k: low.k,
                            high_unitig: unitig.id,
                            reanchor_state: high_state,
                            extension_bases,
                        });
                        terminal_reason = Some(3_u8);
                        break;
                    }
                }
            }
        }

        match terminal_reason {
            Some(0) => summary.dead_stops += 1,
            Some(1) => summary.branch_stops += 1,
            Some(2) => summary.cycle_stops += 1,
            Some(3) => {}
            None => summary.maxed_out += 1,
            Some(_) => unreachable!(),
        }
    }

    Ok((summary, candidates))
}

pub fn build_multik_graph(config: &MultiKConfig) -> Result<MultiKGraph> {
    if config.ks.len() < 2 {
        bail!("multi-k graph requires at least two k values");
    }
    if !(0.0..=60.0).contains(&config.min_mean_quality) {
        bail!("minimum mean quality must be in 0..=60");
    }
    let mut ks = config.ks.clone();
    ks.sort_unstable();
    ks.dedup();
    if ks.len() < 2 {
        bail!("multi-k graph requires at least two distinct k values");
    }
    if ks.iter().any(|&k| k == 0 || k >= MAX_K) {
        bail!("multi-k values must be in 1..{} so k+1 remains representable", MAX_K);
    }

    let total_started = Instant::now();
    let started = Instant::now();
    let (pending, read_pairs) = count_multi_k(
        &config.read1,
        config.read2.as_deref(),
        &ks,
        config.max_pairs,
    )?;
    let count_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let mut layers = Vec::with_capacity(pending.len());
    for layer in pending {
        layers.push(finalize_layer(
            layer,
            config.min_count,
            config.min_fragment_support,
            config.min_mean_quality,
        )?);
    }
    let finalize_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let mut projection_pairs = Vec::new();
    for pair in layers.windows(2) {
        projection_pairs.push(build_projection_pair(&pair[1], &pair[0])?);
    }
    let projection_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let mut adaptive_rescues = Vec::new();
    let mut rescue_candidates = Vec::new();
    for pair in layers.windows(2) {
        let (summary, mut candidates) =
            probe_adaptive_rescue(&pair[1], &pair[0], config.max_rescue_bases)?;
        adaptive_rescues.push(summary);
        rescue_candidates.append(&mut candidates);
    }
    let rescue_seconds = started.elapsed().as_secs_f64();

    let summary = MultiKSummary {
        version: "stage34-layered-multik-v0".to_string(),
        read_pairs,
        ks,
        min_count: config.min_count,
        min_fragment_support: config.min_fragment_support,
        min_mean_quality: config.min_mean_quality,
        max_rescue_bases: config.max_rescue_bases,
        layers: layers.iter().map(|layer| layer.summary.clone()).collect(),
        projections: projection_pairs
            .iter()
            .map(|pair| pair.summary.clone())
            .collect(),
        adaptive_rescues,
        timings_seconds: MultiKTimingSummary {
            count_seconds,
            finalize_seconds,
            projection_seconds,
            rescue_seconds,
            total_seconds: total_started.elapsed().as_secs_f64(),
        },
    };

    Ok(MultiKGraph {
        layers,
        projection_pairs,
        rescue_candidates,
        summary,
    })
}

pub fn write_multik_outputs(graph: &MultiKGraph, output_dir: &Path) -> Result<()> {
    fs::create_dir_all(output_dir)
        .with_context(|| format!("failed to create {}", output_dir.display()))?;

    let summary_path = output_dir.join("multik_summary.json");
    let summary_file = File::create(&summary_path)
        .with_context(|| format!("failed to create {}", summary_path.display()))?;
    serde_json::to_writer_pretty(BufWriter::new(summary_file), &graph.summary)?;

    let mut layers = BufWriter::new(File::create(output_dir.join("layer_summary.tsv"))?);
    writeln!(
        layers,
        "k\tobservations\tedge_observations\tdistinct_kmers\tretained_kmers\tdistinct_kplus1\tdirected_edges\tunitigs\tunitig_n50\tbranching_unitigs"
    )?;
    for layer in &graph.summary.layers {
        writeln!(
            layers,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            layer.k,
            layer.observations,
            layer.edge_observations,
            layer.distinct_kmers,
            layer.retained_kmers,
            layer.distinct_kplus1,
            layer.retained_directed_edges,
            layer.graph.unitigs,
            layer.graph.unitig_n50,
            layer.graph.branching_unitigs
        )?;
    }

    let mut projections = BufWriter::new(File::create(output_dir.join("projection_summary.tsv"))?);
    writeln!(
        projections,
        "high_k\tlow_k\thigh_unitigs\texact_unitigs\tmissing_node_unitigs\tmissing_edge_unitigs\tprojected_low_unitig_visits"
    )?;
    for projection in &graph.summary.projections {
        writeln!(
            projections,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            projection.high_k,
            projection.low_k,
            projection.high_unitigs,
            projection.exact_unitigs,
            projection.missing_node_unitigs,
            projection.missing_edge_unitigs,
            projection.projected_low_unitig_visits
        )?;
    }

    let mut rescues = BufWriter::new(File::create(output_dir.join("adaptive_rescue.tsv"))?);
    writeln!(
        rescues,
        "high_k\tlow_k\tdead_end_unitigs\tlow_anchor_missing\tbranch_stops\tdead_stops\tcycle_stops\tmaxed_out\treanchored\ttotal_extension_bases\tmax_extension_bases"
    )?;
    for rescue in &graph.summary.adaptive_rescues {
        writeln!(
            rescues,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rescue.high_k,
            rescue.low_k,
            rescue.dead_end_unitigs,
            rescue.low_anchor_missing,
            rescue.branch_stops,
            rescue.dead_stops,
            rescue.cycle_stops,
            rescue.maxed_out,
            rescue.reanchored,
            rescue.total_extension_bases,
            rescue.max_extension_bases
        )?;
    }

    let mut candidates = BufWriter::new(File::create(output_dir.join("adaptive_candidates.tsv"))?);
    writeln!(
        candidates,
        "high_k\tlow_k\thigh_unitig\treanchor_state\textension_bases"
    )?;
    for candidate in &graph.rescue_candidates {
        writeln!(
            candidates,
            "{}\t{}\t{}\t{}\t{}",
            candidate.high_k,
            candidate.low_k,
            candidate.high_unitig,
            candidate.reanchor_state,
            candidate.extension_bases
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn simultaneous_roller_matches_single_order_extraction() {
        let sequence = b"ACGTTGCANNTGCAACGTACGATCGTACGTTAGC";
        let orders = [3_usize, 5, 7, 11];
        let mut observed: FxHashMap<usize, Vec<OrientedKmer>> = FxHashMap::default();
        for_each_multi_canonical_kmer(sequence, &orders, |k, item| {
            observed.entry(k).or_default().push(item);
        })
        .unwrap();
        for &k in &orders {
            assert_eq!(observed.get(&k).unwrap(), &canonical_kmers(sequence, k).unwrap());
        }
    }

    #[test]
    fn builds_layered_graph_and_cross_k_projection() {
        let dir = tempfile::tempdir().unwrap();
        let reads = dir.path().join("reads.fastq");
        let mut handle = File::create(&reads).unwrap();
        for index in 0..3 {
            writeln!(
                handle,
                "@r{index}\nACGTTGCAACGTCAGTACGATCGTAGCTAACGTTGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
            )
            .unwrap();
        }
        drop(handle);

        let config = MultiKConfig {
            read1: reads,
            read2: None,
            output_dir: dir.path().join("out"),
            ks: vec![5, 9],
            min_count: 1,
            min_fragment_support: 1,
            min_mean_quality: 0.0,
            max_pairs: None,
            max_rescue_bases: 20,
        };
        let graph = build_multik_graph(&config).unwrap();
        assert_eq!(graph.layers.len(), 2);
        assert!(graph.layers.iter().all(|layer| !layer.raw_graph.keys.is_empty()));
        assert_eq!(graph.projection_pairs.len(), 1);
        assert!(graph.projection_pairs[0].summary.exact_unitigs > 0);
    }
}
