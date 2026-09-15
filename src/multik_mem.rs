use crate::dna::{base_bits, canonical_kmers, KmerKey};
use crate::fastq::for_each_pair;
use crate::graph::{compact_unitigs, summarize, RawGraph};
use crate::multik::{
    AdaptiveRescueCandidate, AdaptiveRescueSummary, MultiKConfig, MultiKGraph, MultiKLayer,
    MultiKLayerSummary, MultiKSummary, MultiKTimingSummary, ProjectionPair, ProjectionPath,
    ProjectionSummary,
};
use crate::multik_fast::build_multik_graph_fast;
use anyhow::{bail, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;
use std::time::Instant;

const MAX_COMPACT_K: usize = 63;
const BLOOM_BITS: usize = 1 << 30;
const BLOOM_WORDS: usize = BLOOM_BITS / 64;
const PROGRESS_PAIRS: usize = 1_000_000;
const HLL_P: u32 = 16;
const HLL_REGISTERS: usize = 1 << HLL_P;

#[derive(Clone, Copy, Debug, Default)]
struct NodeEvidence {
    count: u32,
    fragment_count: u32,
    quality_sum: u64,
    last_fragment: u32,
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

#[derive(Debug)]
struct Bloom {
    words: Vec<u64>,
}

impl Bloom {
    fn new() -> Self {
        Self {
            words: vec![0; BLOOM_WORDS],
        }
    }

    #[inline]
    fn positions(key: u128) -> [usize; 3] {
        let folded = (key as u64) ^ ((key >> 64) as u64).rotate_left(17);
        let h1 = mix64(folded ^ 0x9e3779b97f4a7c15);
        let h2 = mix64(folded ^ 0xbf58476d1ce4e5b9);
        let h3 = mix64(folded ^ 0x94d049bb133111eb);
        let mask = BLOOM_BITS - 1;
        [
            (h1 as usize) & mask,
            (h2 as usize) & mask,
            (h3 as usize) & mask,
        ]
    }

    #[inline]
    fn contains(&self, key: u128) -> bool {
        Self::positions(key).into_iter().all(|bit| {
            let word = bit >> 6;
            let mask = 1_u64 << (bit & 63);
            self.words[word] & mask != 0
        })
    }

    #[inline]
    fn insert(&mut self, key: u128) {
        for bit in Self::positions(key) {
            let word = bit >> 6;
            let mask = 1_u64 << (bit & 63);
            self.words[word] |= mask;
        }
    }
}

#[derive(Debug)]
struct Hll {
    registers: Vec<u8>,
}

impl Hll {
    fn new() -> Self {
        Self {
            registers: vec![0; HLL_REGISTERS],
        }
    }

    #[inline]
    fn insert(&mut self, key: u128) {
        let folded = (key as u64) ^ ((key >> 64) as u64).rotate_left(23);
        let hash = mix64(folded ^ 0xd6e8feb86659fd93);
        let index = (hash >> (64 - HLL_P)) as usize;
        let shifted = hash << HLL_P;
        let rank = shifted.leading_zeros().saturating_add(1).min(64) as u8;
        self.registers[index] = self.registers[index].max(rank);
    }

    fn estimate(&self) -> usize {
        let m = HLL_REGISTERS as f64;
        let alpha = 0.7213 / (1.0 + 1.079 / m);
        let mut reciprocal_sum = 0.0_f64;
        let mut zeros = 0_usize;
        for &register in &self.registers {
            reciprocal_sum += 2_f64.powi(-(register as i32));
            if register == 0 {
                zeros += 1;
            }
        }
        let raw = alpha * m * m / reciprocal_sum;
        let estimate = if raw <= 2.5 * m && zeros > 0 {
            m * (m / zeros as f64).ln()
        } else {
            raw
        };
        estimate.max(0.0).round() as usize
    }
}

#[derive(Debug)]
struct DiscoveryLayer {
    k: usize,
    seen: Bloom,
    repeated: Bloom,
    node_hll: Hll,
    edge_hll: Hll,
    observations: u64,
    edge_observations: u64,
}

impl DiscoveryLayer {
    fn new(k: usize) -> Self {
        Self {
            k,
            seen: Bloom::new(),
            repeated: Bloom::new(),
            node_hll: Hll::new(),
            edge_hll: Hll::new(),
            observations: 0,
            edge_observations: 0,
        }
    }
}

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
    fn push(&mut self, bits: u8) -> Option<(u128, bool)> {
        self.forward = ((self.forward << 2) | u128::from(bits)) & self.mask;
        self.reverse =
            (self.reverse >> 2) | (u128::from(3 - (bits & 0b11)) << self.reverse_shift);
        self.valid += 1;
        if self.valid < self.order {
            None
        } else if self.reverse < self.forward {
            Some((self.reverse, true))
        } else {
            Some((self.forward, false))
        }
    }
}

#[inline]
fn mix64(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

#[inline]
fn packed_mask(k: usize) -> u128 {
    (1_u128 << (2 * k)) - 1
}

fn packed_to_kmer(value: u128) -> KmerKey {
    KmerKey {
        words: [value as u64, (value >> 64) as u64, 0, 0, 0],
    }
}

fn scan_discovery_record(sequence: &[u8], layer: &mut DiscoveryLayer, fragment_keys: &mut Vec<u128>) {
    let mut node = Roller::new(layer.k);
    let mut edge = Roller::new(layer.k + 1);
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            node.reset();
            edge.reset();
            continue;
        };
        if let Some((key, _)) = node.push(bits) {
            fragment_keys.push(key);
            layer.node_hll.insert(key);
            layer.observations = layer.observations.saturating_add(1);
        }
        if let Some((key, _)) = edge.push(bits) {
            layer.edge_hll.insert(key);
            layer.edge_observations = layer.edge_observations.saturating_add(1);
        }
    }
}

fn discover_repeated(
    read1: &Path,
    read2: Option<&Path>,
    ks: &[usize],
    max_pairs: Option<usize>,
) -> Result<(Vec<DiscoveryLayer>, usize)> {
    let mut layers: Vec<DiscoveryLayer> = ks.iter().copied().map(DiscoveryLayer::new).collect();
    let mut per_layer_keys: Vec<Vec<u128>> = ks.iter().map(|_| Vec::with_capacity(512)).collect();
    let read_pairs = for_each_pair(read1, read2, max_pairs, |pair_index, left, right| {
        for (layer, keys) in layers.iter_mut().zip(per_layer_keys.iter_mut()) {
            keys.clear();
            scan_discovery_record(&left.sequence, layer, keys);
            if let Some(right) = right.as_ref() {
                scan_discovery_record(&right.sequence, layer, keys);
            }
            keys.sort_unstable();
            keys.dedup();
            for &key in keys.iter() {
                if layer.seen.contains(key) {
                    layer.repeated.insert(key);
                } else {
                    layer.seen.insert(key);
                }
            }
        }
        let pairs = pair_index + 1;
        if pairs % PROGRESS_PAIRS == 0 {
            eprintln!("stage34 mem discovery progress: {pairs} read pairs");
        }
        Ok(())
    })?;
    Ok((layers, read_pairs))
}

fn exact_node_count(
    read1: &Path,
    read2: Option<&Path>,
    k: usize,
    repeated: &Bloom,
    max_pairs: Option<usize>,
) -> Result<FxHashMap<u128, NodeEvidence>> {
    let mut evidence = FxHashMap::default();
    let mut scan_record = |sequence: &[u8], quality: &[u8], fragment_id: u32| {
        let mut roller = Roller::new(k);
        let mut quality_sum = 0_u64;
        let mut valid = 0_usize;
        for (index, (&base, &quality_byte)) in sequence.iter().zip(quality.iter()).enumerate() {
            let Some(bits) = base_bits(base) else {
                roller.reset();
                quality_sum = 0;
                valid = 0;
                continue;
            };
            quality_sum = quality_sum.saturating_add(u64::from(quality_byte.saturating_sub(33)));
            valid += 1;
            if valid > k {
                quality_sum =
                    quality_sum.saturating_sub(u64::from(quality[index - k].saturating_sub(33)));
            }
            let Some((key, _)) = roller.push(bits) else {
                continue;
            };
            if !repeated.contains(key) {
                continue;
            }
            let entry: &mut NodeEvidence = evidence.entry(key).or_default();
            entry.count = entry.count.saturating_add(1);
            entry.quality_sum = entry.quality_sum.saturating_add(quality_sum);
            if entry.fragment_count == 0 || entry.last_fragment != fragment_id {
                entry.fragment_count = entry.fragment_count.saturating_add(1);
                entry.last_fragment = fragment_id;
            }
        }
    };

    for_each_pair(read1, read2, max_pairs, |pair_index, left, right| {
        let fragment_id = u32::try_from(pair_index + 1)
            .map_err(|_| anyhow::anyhow!("too many read pairs for u32 fragment ids"))?;
        scan_record(&left.sequence, &left.quality, fragment_id);
        if let Some(right) = right.as_ref() {
            scan_record(&right.sequence, &right.quality, fragment_id);
        }
        let pairs = pair_index + 1;
        if pairs % PROGRESS_PAIRS == 0 {
            eprintln!("stage34 mem exact k={k} progress: {pairs} read pairs");
        }
        Ok(())
    })?;
    Ok(evidence)
}

#[inline]
fn pack_edge(source: u32, target: u32) -> u64 {
    (u64::from(source) << 32) | u64::from(target)
}

#[inline]
fn unpack_edge(edge: u64) -> (u32, u32) {
    ((edge >> 32) as u32, edge as u32)
}

fn add_record_edges(
    sequence: &[u8],
    k: usize,
    index: &FxHashMap<u128, u32>,
    edge_counts: &mut FxHashMap<u64, u32>,
) {
    let mut roller = Roller::new(k);
    let mut previous: Option<(u128, bool)> = None;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous = None;
            continue;
        };
        let Some(current) = roller.push(bits) else {
            continue;
        };
        if let Some((left_key, left_reverse)) = previous {
            let (right_key, right_reverse) = current;
            if let (Some(&left_node), Some(&right_node)) =
                (index.get(&left_key), index.get(&right_key))
            {
                let source = left_node * 2 + u32::from(left_reverse);
                let target = right_node * 2 + u32::from(right_reverse);
                let forward = edge_counts.entry(pack_edge(source, target)).or_insert(0);
                *forward = forward.saturating_add(1);
                let reverse_source = RawGraph::reverse_state(target);
                let reverse_target = RawGraph::reverse_state(source);
                let reverse = edge_counts
                    .entry(pack_edge(reverse_source, reverse_target))
                    .or_insert(0);
                *reverse = reverse.saturating_add(1);
            }
        }
        previous = Some(current);
    }
}

fn build_layer_memory_bounded(
    read1: &Path,
    read2: Option<&Path>,
    discovery: DiscoveryLayer,
    min_count: u32,
    min_fragment_support: u32,
    min_mean_quality: f32,
    max_pairs: Option<usize>,
) -> Result<MultiKLayer> {
    let DiscoveryLayer {
        k,
        seen,
        repeated,
        node_hll,
        edge_hll,
        observations,
        edge_observations,
    } = discovery;
    drop(seen);

    let exact_started = Instant::now();
    let evidence = exact_node_count(read1, read2, k, &repeated, max_pairs)?;
    drop(repeated);
    let candidate_kmers = evidence.len();
    eprintln!(
        "stage34 mem k={k}: exact repeated-candidate table {candidate_kmers} entries in {:.3}s",
        exact_started.elapsed().as_secs_f64()
    );

    let mut packed_sorted: Vec<u128> = evidence
        .iter()
        .filter_map(|(key, value)| {
            (value.count >= min_count
                && value.fragment_count >= min_fragment_support
                && value.mean_quality(k) >= min_mean_quality)
                .then_some(*key)
        })
        .collect();
    if packed_sorted.len() > (u32::MAX as usize) / 2 {
        bail!("k={k} graph has too many canonical nodes");
    }
    packed_sorted.sort_unstable();
    let counts: Vec<u32> = packed_sorted
        .iter()
        .map(|key| evidence.get(key).map_or(0, |value| value.count))
        .collect();
    drop(evidence);

    let packed_index: FxHashMap<u128, u32> = packed_sorted
        .iter()
        .enumerate()
        .map(|(node, key)| (*key, node as u32))
        .collect();
    eprintln!(
        "stage34 mem k={k}: retained {} nodes; rejected candidate evidence released before edge scan",
        packed_sorted.len()
    );

    let mut edge_counts: FxHashMap<u64, u32> = FxHashMap::default();
    for_each_pair(read1, read2, max_pairs, |pair_index, left, right| {
        add_record_edges(&left.sequence, k, &packed_index, &mut edge_counts);
        if let Some(right) = right.as_ref() {
            add_record_edges(&right.sequence, k, &packed_index, &mut edge_counts);
        }
        let pairs = pair_index + 1;
        if pairs % PROGRESS_PAIRS == 0 {
            eprintln!("stage34 mem edge k={k} progress: {pairs} read pairs");
        }
        Ok(())
    })?;

    let mut singleton_solid_edges = 0_usize;
    let mut edges = Vec::with_capacity(edge_counts.len());
    for (edge, support) in edge_counts {
        if support < min_count {
            singleton_solid_edges += 1;
        }
        edges.push(unpack_edge(edge));
    }
    edges.sort_unstable();
    drop(packed_index);

    let state_count = packed_sorted.len() * 2;
    let mut out_offsets = vec![0_u64; state_count + 1];
    let mut indegree = vec![0_u32; state_count];
    for &(source, target) in &edges {
        out_offsets[source as usize + 1] += 1;
        indegree[target as usize] = indegree[target as usize].saturating_add(1);
    }
    for index in 1..out_offsets.len() {
        out_offsets[index] += out_offsets[index - 1];
    }
    let out_targets: Vec<u32> = edges.into_iter().map(|(_, target)| target).collect();
    let keys: Vec<KmerKey> = packed_sorted.into_iter().map(packed_to_kmer).collect();
    let raw_graph = RawGraph {
        k,
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
        k,
        observations,
        edge_observations,
        distinct_kmers: node_hll.estimate(),
        retained_kmers: raw_graph.keys.len(),
        distinct_kplus1: edge_hll.estimate(),
        retained_directed_edges: raw_graph.out_targets.len(),
        graph: graph_summary,
    };
    Ok(MultiKLayer {
        k,
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

pub fn build_multik_graph_memory_bounded(config: &MultiKConfig) -> Result<MultiKGraph> {
    if config.min_count <= 1 || config.min_fragment_support <= 1 {
        eprintln!(
            "stage34 memory-bounded prefilter requires min-count>=2 and min-fragment-support>=2; using v1 exact path"
        );
        return build_multik_graph_fast(config);
    }
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
    let count_started = Instant::now();
    let (mut discovery, read_pairs) = discover_repeated(
        &config.read1,
        config.read2.as_deref(),
        &ks,
        config.max_pairs,
    )?;
    let discovery_seconds = count_started.elapsed().as_secs_f64();
    eprintln!(
        "stage34 mem repeat discovery complete: {read_pairs} pairs in {discovery_seconds:.3}s"
    );

    let finalize_started = Instant::now();
    let mut layers = Vec::with_capacity(discovery.len());
    for layer in discovery.drain(..) {
        let k = layer.k;
        let started = Instant::now();
        layers.push(build_layer_memory_bounded(
            &config.read1,
            config.read2.as_deref(),
            layer,
            config.min_count,
            config.min_fragment_support,
            config.min_mean_quality,
            config.max_pairs,
        )?);
        eprintln!(
            "stage34 mem k={k} layer complete in {:.3}s",
            started.elapsed().as_secs_f64()
        );
    }
    let finalize_seconds = finalize_started.elapsed().as_secs_f64();

    let projection_started = Instant::now();
    let mut projection_pairs = Vec::new();
    for pair in layers.windows(2) {
        projection_pairs.push(build_projection_pair(&pair[1], &pair[0])?);
    }
    let projection_seconds = projection_started.elapsed().as_secs_f64();

    let rescue_started = Instant::now();
    let mut adaptive_rescues = Vec::new();
    let mut rescue_candidates = Vec::new();
    for pair in layers.windows(2) {
        let (summary, mut candidates) =
            probe_adaptive_rescue(&pair[1], &pair[0], config.max_rescue_bases)?;
        adaptive_rescues.push(summary);
        rescue_candidates.append(&mut candidates);
    }
    let rescue_seconds = rescue_started.elapsed().as_secs_f64();

    let summary = MultiKSummary {
        version: "stage34-layered-multik-v2-memory-bounded".to_string(),
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
            count_seconds: discovery_seconds,
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
    fn bloom_keeps_inserted_keys() {
        let mut bloom = Bloom::new();
        for key in [1_u128, 7, 123456789, u64::MAX as u128 + 99] {
            bloom.insert(key);
            assert!(bloom.contains(key));
        }
    }

    #[test]
    fn memory_bounded_matches_fast_retained_graph_on_repeated_fixture() {
        let dir = tempfile::tempdir().unwrap();
        let reads = dir.path().join("reads.fastq");
        let sequence = "ACGTTGCAACGTCAGTACGATCGTAGCTAACGTTGCA";
        let mut handle = File::create(&reads).unwrap();
        for index in 0..4 {
            writeln!(
                handle,
                "@r{index}\n{sequence}\n+\n{}",
                "I".repeat(sequence.len())
            )
            .unwrap();
        }
        drop(handle);
        let config = MultiKConfig {
            read1: reads,
            read2: None,
            output_dir: dir.path().join("out"),
            ks: vec![5, 9],
            min_count: 2,
            min_fragment_support: 2,
            min_mean_quality: 20.0,
            max_pairs: None,
            max_rescue_bases: 20,
        };
        let fast = build_multik_graph_fast(&config).unwrap();
        let bounded = build_multik_graph_memory_bounded(&config).unwrap();
        assert_eq!(bounded.layers.len(), fast.layers.len());
        for (left, right) in bounded.layers.iter().zip(fast.layers.iter()) {
            assert_eq!(left.raw_graph.keys, right.raw_graph.keys);
            assert_eq!(left.raw_graph.out_offsets, right.raw_graph.out_offsets);
            assert_eq!(left.raw_graph.out_targets, right.raw_graph.out_targets);
        }
    }
}
