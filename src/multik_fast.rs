use crate::dna::{base_bits, canonical_kmers, KmerKey};
use crate::fastq::for_each_pair;
use crate::graph::{compact_unitigs, summarize, RawGraph};
use crate::multik::{
    AdaptiveRescueCandidate, AdaptiveRescueSummary, MultiKConfig, MultiKGraph, MultiKLayer,
    MultiKLayerSummary, MultiKSummary, MultiKTimingSummary, ProjectionPair, ProjectionPath,
    ProjectionSummary,
};
use anyhow::{bail, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;
use std::sync::{mpsc, Arc};
use std::time::Instant;

const MAX_COMPACT_K: usize = 63;
const BATCH_PAIRS: usize = 4096;
const PROGRESS_PAIRS: usize = 1_000_000;

#[derive(Clone, Copy, Debug, Default)]
struct NodeEvidence {
    count: u32,
    fragment_count: u32,
    quality_sum: u64,
    last_fragment: usize,
}

impl NodeEvidence {
    fn mean_quality(self, k: usize) -> f32 {
        if self.count == 0 {
            0.0
        } else {
            self.quality_sum as f32 / (self.count as f32 * k as f32)
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct EdgeEvidence {
    count: u32,
    fragment_count: u32,
    last_fragment: usize,
}

#[derive(Debug)]
struct PendingLayer {
    k: usize,
    nodes: FxHashMap<u128, NodeEvidence>,
    edges: FxHashMap<u128, EdgeEvidence>,
    observations: u64,
    edge_observations: u64,
}

impl PendingLayer {
    fn new(k: usize) -> Self {
        Self {
            k,
            nodes: FxHashMap::default(),
            edges: FxHashMap::default(),
            observations: 0,
            edge_observations: 0,
        }
    }
}

#[derive(Clone, Debug)]
struct PairRecord {
    fragment_id: usize,
    left_sequence: Vec<u8>,
    left_quality: Vec<u8>,
    right: Option<(Vec<u8>, Vec<u8>)>,
}

type PairBatch = Arc<Vec<PairRecord>>;

#[derive(Clone, Copy, Debug)]
struct Roller {
    order: usize,
    mask: u128,
    reverse_shift: usize,
    forward: u128,
    reverse: u128,
    valid: usize,
}

impl Roller {
    fn new(order: usize) -> Self {
        Self {
            order,
            mask: packed_mask(order),
            reverse_shift: 2 * (order - 1),
            forward: 0,
            reverse: 0,
            valid: 0,
        }
    }

    #[inline]
    fn reset(&mut self) {
        self.forward = 0;
        self.reverse = 0;
        self.valid = 0;
    }

    #[inline]
    fn push(&mut self, bits: u8) -> Option<u128> {
        self.forward = ((self.forward << 2) | u128::from(bits)) & self.mask;
        self.reverse = (self.reverse >> 2)
            | (u128::from(3 - (bits & 0b11)) << self.reverse_shift);
        self.valid += 1;
        (self.valid >= self.order).then_some(self.forward.min(self.reverse))
    }
}

#[inline]
fn packed_mask(k: usize) -> u128 {
    (1_u128 << (2 * k)) - 1
}

fn reverse_complement_packed(mut value: u128, k: usize) -> u128 {
    let mut output = 0_u128;
    for _ in 0..k {
        let bits = (value & 0b11) as u8;
        value >>= 2;
        output = (output << 2) | u128::from(3 - bits);
    }
    output
}

fn canonical_packed(value: u128, k: usize) -> (u128, bool) {
    let reverse = reverse_complement_packed(value, k);
    if reverse < value {
        (reverse, true)
    } else {
        (value, false)
    }
}

fn packed_to_kmer(value: u128) -> KmerKey {
    KmerKey {
        words: [value as u64, (value >> 64) as u64, 0, 0, 0],
    }
}

#[inline]
fn observe_node(
    pending: &mut PendingLayer,
    key: u128,
    fragment_id: usize,
    quality_sum: u64,
) {
    let entry = pending.nodes.entry(key).or_default();
    entry.count = entry.count.saturating_add(1);
    entry.quality_sum = entry.quality_sum.saturating_add(quality_sum);
    if entry.fragment_count == 0 || entry.last_fragment != fragment_id {
        entry.fragment_count = entry.fragment_count.saturating_add(1);
        entry.last_fragment = fragment_id;
    }
    pending.observations = pending.observations.saturating_add(1);
}

#[inline]
fn observe_edge(pending: &mut PendingLayer, key: u128, fragment_id: usize) {
    let entry = pending.edges.entry(key).or_default();
    entry.count = entry.count.saturating_add(1);
    if entry.fragment_count == 0 || entry.last_fragment != fragment_id {
        entry.fragment_count = entry.fragment_count.saturating_add(1);
        entry.last_fragment = fragment_id;
    }
    pending.edge_observations = pending.edge_observations.saturating_add(1);
}

fn scan_record_layer(
    sequence: &[u8],
    quality: &[u8],
    fragment_id: usize,
    pending: &mut PendingLayer,
) {
    let k = pending.k;
    let mut node = Roller::new(k);
    let mut edge = Roller::new(k + 1);
    let mut quality_sum = 0_u64;
    let mut valid = 0_usize;

    for (index, (&base, &quality_byte)) in sequence.iter().zip(quality.iter()).enumerate() {
        let Some(bits) = base_bits(base) else {
            node.reset();
            edge.reset();
            quality_sum = 0;
            valid = 0;
            continue;
        };

        let q = u64::from(quality_byte.saturating_sub(33));
        quality_sum = quality_sum.saturating_add(q);
        valid += 1;
        if valid > k {
            quality_sum = quality_sum.saturating_sub(u64::from(quality[index - k].saturating_sub(33)));
        }

        if let Some(key) = node.push(bits) {
            observe_node(pending, key, fragment_id, quality_sum);
        }
        if let Some(key) = edge.push(bits) {
            observe_edge(pending, key, fragment_id);
        }
    }
}

fn process_batch_for_layer(batch: &PairBatch, pending: &mut PendingLayer) {
    for pair in batch.iter() {
        scan_record_layer(
            &pair.left_sequence,
            &pair.left_quality,
            pair.fragment_id,
            pending,
        );
        if let Some((sequence, quality)) = &pair.right {
            scan_record_layer(sequence, quality, pair.fragment_id, pending);
        }
    }
}

fn count_multi_k_parallel(
    read1: &Path,
    read2: Option<&Path>,
    ks: &[usize],
    max_pairs: Option<usize>,
) -> Result<(Vec<PendingLayer>, usize)> {
    let mut layers: Vec<Option<PendingLayer>> = (0..ks.len()).map(|_| None).collect();
    let mut read_pairs = 0_usize;

    std::thread::scope(|scope| -> Result<()> {
        let mut senders = Vec::with_capacity(ks.len());
        let mut handles = Vec::with_capacity(ks.len());
        for &k in ks {
            let (sender, receiver) = mpsc::sync_channel::<PairBatch>(2);
            senders.push(sender);
            handles.push(scope.spawn(move || {
                let mut pending = PendingLayer::new(k);
                while let Ok(batch) = receiver.recv() {
                    process_batch_for_layer(&batch, &mut pending);
                }
                pending
            }));
        }

        let mut batch = Vec::with_capacity(BATCH_PAIRS);
        read_pairs = for_each_pair(read1, read2, max_pairs, |pair_index, left, right| {
            let fragment_id = pair_index + 1;
            batch.push(PairRecord {
                fragment_id,
                left_sequence: left.sequence,
                left_quality: left.quality,
                right: right.map(|record| (record.sequence, record.quality)),
            });

            if batch.len() >= BATCH_PAIRS {
                let shared = Arc::new(std::mem::replace(
                    &mut batch,
                    Vec::with_capacity(BATCH_PAIRS),
                ));
                for sender in &senders {
                    sender
                        .send(Arc::clone(&shared))
                        .map_err(|_| anyhow::anyhow!("multi-k layer worker stopped unexpectedly"))?;
                }
            }
            if fragment_id % PROGRESS_PAIRS == 0 {
                eprintln!("stage34 count progress: {fragment_id} read pairs");
            }
            Ok(())
        })?;

        if !batch.is_empty() {
            let shared = Arc::new(batch);
            for sender in &senders {
                sender
                    .send(Arc::clone(&shared))
                    .map_err(|_| anyhow::anyhow!("multi-k layer worker stopped unexpectedly"))?;
            }
        }
        drop(senders);

        for (index, handle) in handles.into_iter().enumerate() {
            layers[index] = Some(
                handle
                    .join()
                    .map_err(|_| anyhow::anyhow!("multi-k layer worker panicked"))?,
            );
        }
        Ok(())
    })?;

    Ok((
        layers
            .into_iter()
            .map(|layer| layer.expect("all layer workers joined"))
            .collect(),
        read_pairs,
    ))
}

fn finalize_layer(
    pending: PendingLayer,
    min_count: u32,
    min_fragment_support: u32,
    min_mean_quality: f32,
) -> Result<MultiKLayer> {
    let distinct_kmers = pending.nodes.len();
    let distinct_kplus1 = pending.edges.len();
    let mut packed_sorted: Vec<u128> = pending
        .nodes
        .iter()
        .filter_map(|(key, value)| {
            (value.count >= min_count
                && value.fragment_count >= min_fragment_support
                && value.mean_quality(pending.k) >= min_mean_quality)
                .then_some(*key)
        })
        .collect();
    if packed_sorted.len() > (u32::MAX as usize) / 2 {
        bail!("k={} graph has too many canonical nodes", pending.k);
    }
    packed_sorted.sort_unstable();

    let packed_index: FxHashMap<u128, u32> = packed_sorted
        .iter()
        .enumerate()
        .map(|(node, key)| (*key, node as u32))
        .collect();
    let keys: Vec<KmerKey> = packed_sorted.iter().copied().map(packed_to_kmer).collect();
    let counts: Vec<u32> = packed_sorted
        .iter()
        .map(|key| pending.nodes.get(key).map_or(0, |value| value.count))
        .collect();

    let target_mask = packed_mask(pending.k);
    let mut edges = Vec::new();
    let mut singleton_solid_edges = 0_usize;
    for (edge_key, evidence) in &pending.edges {
        let source_forward = edge_key >> 2;
        let target_forward = edge_key & target_mask;
        let (source_key, source_reverse) = canonical_packed(source_forward, pending.k);
        let (target_key, target_reverse) = canonical_packed(target_forward, pending.k);
        let (Some(&source_node), Some(&target_node)) =
            (packed_index.get(&source_key), packed_index.get(&target_key))
        else {
            continue;
        };
        let source = source_node * 2 + u32::from(source_reverse);
        let target = target_node * 2 + u32::from(target_reverse);
        edges.push((source, target));
        edges.push((
            RawGraph::reverse_state(target),
            RawGraph::reverse_state(source),
        ));
        if evidence.count < min_count || evidence.fragment_count < min_fragment_support {
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
    for index in 1..out_offsets.len() {
        out_offsets[index] += out_offsets[index - 1];
    }
    let out_targets = edges.into_iter().map(|(_, target)| target).collect();
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

pub fn build_multik_graph_fast(config: &MultiKConfig) -> Result<MultiKGraph> {
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
    if ks.iter().any(|&k| k == 0 || k >= MAX_COMPACT_K) {
        bail!(
            "Stage34 compact multi-k values must be in 1..{} so k+1 fits in u128",
            MAX_COMPACT_K
        );
    }

    let total_started = Instant::now();
    let started = Instant::now();
    let (pending, read_pairs) = count_multi_k_parallel(
        &config.read1,
        config.read2.as_deref(),
        &ks,
        config.max_pairs,
    )?;
    let count_seconds = started.elapsed().as_secs_f64();
    eprintln!("stage34 count complete: {read_pairs} pairs in {count_seconds:.3}s");

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
        version: "stage34-layered-multik-v1-fast".to_string(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn physical_fragment_support_counts_pair_once() {
        let mut pending = PendingLayer::new(3);
        scan_record_layer(b"AAAAAA", b"IIIIII", 7, &mut pending);
        scan_record_layer(b"AAAAAA", b"IIIIII", 7, &mut pending);
        let value = pending.nodes.values().next().unwrap();
        assert!(value.count > 1);
        assert_eq!(value.fragment_count, 1);
        scan_record_layer(b"AAAAAA", b"IIIIII", 8, &mut pending);
        let value = pending.nodes.values().next().unwrap();
        assert_eq!(value.fragment_count, 2);
    }

    #[test]
    fn fast_builder_constructs_layers() {
        let dir = tempfile::tempdir().unwrap();
        let reads = dir.path().join("reads.fastq");
        let sequence = "ACGTTGCAACGTCAGTACGATCGTAGCTAACGTTGCA";
        let mut handle = File::create(&reads).unwrap();
        for index in 0..4 {
            writeln!(handle, "@r{index}\n{sequence}\n+\n{}", "I".repeat(sequence.len()))
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
        let graph = build_multik_graph_fast(&config).unwrap();
        assert_eq!(graph.layers.len(), 2);
        assert_eq!(graph.summary.version, "stage34-layered-multik-v1-fast");
        assert!(graph.layers.iter().all(|layer| !layer.raw_graph.keys.is_empty()));
    }
}
