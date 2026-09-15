use crate::dna::base_bits;
use crate::fastq::{for_each_pair, FastqRecord};
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

const BACKBONE_MAX_K: usize = 62;
const WIDE_MAX_K: usize = 255;
const WIDE_WORDS: usize = 8;
const BLOOM_BITS: usize = 1 << 30;
const BLOOM_WORDS: usize = BLOOM_BITS / 64;
const BATCH_PAIRS: usize = 8192;
const PREFIX_BITS: usize = 23;
const NO_NEIGHBORHOOD: u64 = 0;

#[derive(Clone, Debug)]
pub struct HybridConfig {
    pub read1: PathBuf,
    pub read2: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub backbone_k: usize,
    pub local_start_k: usize,
    pub max_local_k: usize,
    pub min_count: u32,
    pub min_fragment_support: u32,
    pub min_mean_quality: f32,
    pub neighborhood_radius: usize,
    pub max_neighborhood_states: usize,
    pub max_neighborhoods: usize,
    pub max_pairs: Option<usize>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalKSummary {
    pub k: usize,
    pub retained_kmers: usize,
    pub directed_edges: usize,
    pub ambiguous_states: usize,
    pub dead_end_states: usize,
    pub projection_missing_nodes: usize,
    pub projection_missing_edges: usize,
    pub resolved: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct HybridNeighborhoodSummary {
    pub id: u32,
    pub seed_state: u32,
    pub states: usize,
    pub ambiguous_states: usize,
    pub routed_fragments: usize,
    pub supported_max_k: usize,
    pub selected_k: Option<usize>,
    pub resolved: bool,
    pub candidates: Vec<LocalKSummary>,
}

#[derive(Clone, Debug, Serialize)]
pub struct HybridSummary {
    pub version: String,
    pub read_pairs: usize,
    pub threads: usize,
    pub backbone_k: usize,
    pub backbone_retained_kmers: usize,
    pub backbone_directed_edges: usize,
    pub backbone_ambiguous_states: usize,
    pub neighborhoods: usize,
    pub routed_fragment_copies: usize,
    pub local_start_k: usize,
    pub max_local_k: usize,
    pub backbone_seconds: f64,
    pub routing_seconds: f64,
    pub local_refinement_seconds: f64,
    pub total_seconds: f64,
    pub local: Vec<HybridNeighborhoodSummary>,
}

#[derive(Clone, Debug, Serialize)]
pub struct HybridRunStats {
    pub summary: HybridSummary,
}

#[derive(Clone, Copy, Debug, Default)]
struct ExactEvidence {
    count: u32,
    fragment_count: u32,
    quality_sum: u64,
}

impl ExactEvidence {
    fn mean_quality(self, k: usize) -> f32 {
        if self.count == 0 {
            0.0
        } else {
            self.quality_sum as f32 / (self.count as f32 * k as f32)
        }
    }

    fn add(&mut self, other: Self) {
        self.count = self.count.saturating_add(other.count);
        self.fragment_count = self.fragment_count.saturating_add(other.fragment_count);
        self.quality_sum = self.quality_sum.saturating_add(other.quality_sum);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct LocalEvidence {
    exact: ExactEvidence,
    last_fragment: u64,
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

    fn release(&mut self) {
        self.words = Vec::new();
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

#[derive(Clone, Copy, Debug)]
struct Roller128 {
    k: usize,
    mask: u128,
    reverse_shift: usize,
    forward: u128,
    reverse: u128,
    valid: usize,
}

impl Roller128 {
    fn new(k: usize) -> Self {
        Self {
            k,
            mask: (1_u128 << (2 * k)) - 1,
            reverse_shift: 2 * (k - 1),
            forward: 0,
            reverse: 0,
            valid: 0,
        }
    }

    fn reset(&mut self) {
        self.forward = 0;
        self.reverse = 0;
        self.valid = 0;
    }

    #[inline]
    fn push(&mut self, bits: u8) -> Option<(u128, bool)> {
        self.forward = ((self.forward << 2) | u128::from(bits)) & self.mask;
        self.reverse = (self.reverse >> 2) | (u128::from(3 - (bits & 0b11)) << self.reverse_shift);
        self.valid += 1;
        if self.valid < self.k {
            None
        } else if self.reverse < self.forward {
            Some((self.reverse, true))
        } else {
            Some((self.forward, false))
        }
    }
}

#[derive(Debug)]
struct PairLite {
    fragment_id: u64,
    left_sequence: Vec<u8>,
    left_quality: Vec<u8>,
    right_sequence: Option<Vec<u8>>,
    right_quality: Option<Vec<u8>>,
}

impl PairLite {
    fn new(index: usize, left: FastqRecord, right: Option<FastqRecord>) -> Self {
        let (right_sequence, right_quality) = match right {
            Some(record) => (Some(record.sequence), Some(record.quality)),
            None => (None, None),
        };
        Self {
            fragment_id: index as u64 + 1,
            left_sequence: left.sequence,
            left_quality: left.quality,
            right_sequence,
            right_quality,
        }
    }
}

fn process_pair_batches<F>(
    read1: &Path,
    read2: Option<&Path>,
    max_pairs: Option<usize>,
    mut process: F,
) -> Result<usize>
where
    F: FnMut(&[PairLite]) -> Result<()>,
{
    let mut batch = Vec::with_capacity(BATCH_PAIRS);
    let read_pairs = for_each_pair(read1, read2, max_pairs, |index, left, right| {
        batch.push(PairLite::new(index, left, right));
        if batch.len() == BATCH_PAIRS {
            process(&batch)?;
            batch.clear();
        }
        Ok(())
    })?;
    if !batch.is_empty() {
        process(&batch)?;
    }
    Ok(read_pairs)
}

#[derive(Debug)]
struct PrefixIndex {
    shift: usize,
    mask: usize,
    offsets: Vec<u32>,
}

impl PrefixIndex {
    fn build(keys: &[u128], k: usize) -> Self {
        let bits = PREFIX_BITS.min(2 * k);
        let bucket_count = 1_usize << bits;
        let shift = 2 * k - bits;
        let mask = bucket_count - 1;
        let mut offsets = vec![0_u32; bucket_count + 1];
        for &key in keys {
            let prefix = ((key >> shift) as usize) & mask;
            offsets[prefix + 1] += 1;
        }
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        Self {
            shift,
            mask,
            offsets,
        }
    }

    #[inline]
    fn find(&self, keys: &[u128], key: u128) -> Option<u32> {
        let prefix = ((key >> self.shift) as usize) & self.mask;
        let start = self.offsets[prefix] as usize;
        let end = self.offsets[prefix + 1] as usize;
        keys[start..end]
            .binary_search(&key)
            .ok()
            .map(|local| (start + local) as u32)
    }
}

#[derive(Debug)]
struct BackboneGraph {
    k: usize,
    keys: Vec<u128>,
    prefix: PrefixIndex,
    out_offsets: Vec<u32>,
    out_targets: Vec<u32>,
    indegree: Vec<u8>,
}

impl BackboneGraph {
    fn state_count(&self) -> usize {
        self.keys.len() * 2
    }

    fn out_range(&self, state: u32) -> std::ops::Range<usize> {
        let state = state as usize;
        self.out_offsets[state] as usize..self.out_offsets[state + 1] as usize
    }

    fn outdegree(&self, state: u32) -> usize {
        self.out_range(state).len()
    }

    fn find_node(&self, key: u128) -> Option<u32> {
        self.prefix.find(&self.keys, key)
    }

    fn has_edge(&self, source: u32, target: u32) -> bool {
        self.out_targets[self.out_range(source)]
            .binary_search(&target)
            .is_ok()
    }

    fn state_first_base_bits(&self, state: u32) -> u8 {
        let key = self.keys[(state / 2) as usize];
        if state & 1 == 0 {
            ((key >> (2 * (self.k - 1))) & 0b11) as u8
        } else {
            3 - (key & 0b11) as u8
        }
    }

    fn target_for_base(&self, source: u32, base: u8) -> Option<u32> {
        let key = self.keys[(source / 2) as usize];
        let key_rc = reverse_complement_u128(key, self.k);
        let oriented = if source & 1 == 0 { key } else { key_rc };
        let next = ((oriented << 2) | u128::from(base)) & ((1_u128 << (2 * self.k)) - 1);
        let next_rc = reverse_complement_u128(next, self.k);
        let (canonical, reverse) = if next_rc < next {
            (next_rc, true)
        } else {
            (next, false)
        };
        self.find_node(canonical)
            .map(|node| node * 2 + u32::from(reverse))
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

fn reverse_complement_u128(mut value: u128, k: usize) -> u128 {
    let mut output = 0_u128;
    for _ in 0..k {
        let bits = (value & 0b11) as u8;
        value >>= 2;
        output = (output << 2) | u128::from(3 - bits);
    }
    output
}

fn discover_backbone(config: &HybridConfig) -> Result<(Bloom, usize)> {
    let mut seen = Bloom::new();
    let mut repeated = Bloom::new();
    let k = config.backbone_k;
    let read_pairs = for_each_pair(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |_index, left, right| {
            let mut keys = Vec::with_capacity(512);
            collect_u128_keys(&left.sequence, k, &mut keys);
            if let Some(right) = right {
                collect_u128_keys(&right.sequence, k, &mut keys);
            }
            keys.sort_unstable();
            keys.dedup();
            for key in keys {
                if seen.contains(key) {
                    repeated.insert(key);
                } else {
                    seen.insert(key);
                }
            }
            Ok(())
        },
    )?;
    seen.release();
    Ok((repeated, read_pairs))
}

fn collect_u128_keys(sequence: &[u8], k: usize, output: &mut Vec<u128>) {
    let mut roller = Roller128::new(k);
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            continue;
        };
        if let Some((key, _)) = roller.push(bits) {
            output.push(key);
        }
    }
}

fn scan_exact_u128(
    sequence: &[u8],
    quality: &[u8],
    k: usize,
    repeated: &Bloom,
    fragment_id: u64,
    evidence: &mut FxHashMap<u128, LocalEvidence>,
) {
    let mut roller = Roller128::new(k);
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
        let entry = evidence.entry(key).or_default();
        entry.exact.count = entry.exact.count.saturating_add(1);
        entry.exact.quality_sum = entry.exact.quality_sum.saturating_add(quality_sum);
        if entry.last_fragment != fragment_id {
            entry.last_fragment = fragment_id;
            entry.exact.fragment_count = entry.exact.fragment_count.saturating_add(1);
        }
    }
}

fn merge_local_maps(
    mut left: FxHashMap<u128, LocalEvidence>,
    right: FxHashMap<u128, LocalEvidence>,
) -> FxHashMap<u128, LocalEvidence> {
    if left.len() < right.len() {
        return merge_local_maps(right, left);
    }
    for (key, value) in right {
        left.entry(key).or_default().exact.add(value.exact);
    }
    left
}

fn exact_backbone(config: &HybridConfig, repeated: &Bloom, pool: &ThreadPool) -> Result<Vec<u128>> {
    let mut global: FxHashMap<u128, ExactEvidence> = FxHashMap::default();
    process_pair_batches(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |batch| {
            let workers = pool.current_num_threads().max(1);
            let chunk = batch.len().div_ceil(workers).max(1);
            let local = pool.install(|| {
                batch
                    .par_chunks(chunk)
                    .map(|pairs| {
                        let mut local = FxHashMap::default();
                        for pair in pairs {
                            scan_exact_u128(
                                &pair.left_sequence,
                                &pair.left_quality,
                                config.backbone_k,
                                repeated,
                                pair.fragment_id,
                                &mut local,
                            );
                            if let (Some(sequence), Some(quality)) = (
                                pair.right_sequence.as_deref(),
                                pair.right_quality.as_deref(),
                            ) {
                                scan_exact_u128(
                                    sequence,
                                    quality,
                                    config.backbone_k,
                                    repeated,
                                    pair.fragment_id,
                                    &mut local,
                                );
                            }
                        }
                        local
                    })
                    .reduce(FxHashMap::default, merge_local_maps)
            });
            for (key, value) in local {
                global.entry(key).or_default().add(value.exact);
            }
            Ok(())
        },
    )?;
    let mut keys: Vec<u128> = global
        .into_iter()
        .filter_map(|(key, evidence)| {
            (evidence.count >= config.min_count
                && evidence.fragment_count >= config.min_fragment_support
                && evidence.mean_quality(config.backbone_k) >= config.min_mean_quality)
                .then_some(key)
        })
        .collect();
    keys.sort_unstable();
    Ok(keys)
}

#[inline]
fn lane_support(byte: u8, base: u8) -> u8 {
    (byte >> (2 * (base & 0b11))) & 0b11
}

#[inline]
fn bump_lane(byte: &mut u8, base: u8) {
    let shift = 2 * (base & 0b11);
    let value = ((*byte >> shift) & 0b11).saturating_add(1).min(2);
    *byte = (*byte & !(0b11 << shift)) | (value << shift);
}

fn scan_backbone_edges(sequence: &[u8], graph: &BackboneGraph, support: &mut [u8]) {
    let mut roller = Roller128::new(graph.k);
    let mut previous: Option<u32> = None;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous = None;
            continue;
        };
        let Some((key, reverse)) = roller.push(bits) else {
            continue;
        };
        let current = graph
            .find_node(key)
            .map(|node| node * 2 + u32::from(reverse));
        if let (Some(source), Some(target)) = (previous, current) {
            bump_lane(&mut support[source as usize], bits);
            let reverse_source = target ^ 1;
            let reverse_base = 3 - graph.state_first_base_bits(source);
            bump_lane(&mut support[reverse_source as usize], reverse_base);
        }
        previous = current;
    }
}

fn count_backbone_edges(
    config: &HybridConfig,
    graph: &BackboneGraph,
    pool: &ThreadPool,
) -> Result<Vec<u8>> {
    let workers = pool.current_num_threads().max(1);
    let mut worker_support = vec![vec![0_u8; graph.state_count()]; workers];
    process_pair_batches(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |batch| {
            let batch_len = batch.len();
            pool.install(|| {
                worker_support
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(worker, support)| {
                        let start = batch_len * worker / workers;
                        let end = batch_len * (worker + 1) / workers;
                        for pair in &batch[start..end] {
                            scan_backbone_edges(&pair.left_sequence, graph, support);
                            if let Some(sequence) = pair.right_sequence.as_deref() {
                                scan_backbone_edges(sequence, graph, support);
                            }
                        }
                    });
            });
            Ok(())
        },
    )?;
    let mut merged = worker_support.remove(0);
    for source in worker_support {
        pool.install(|| {
            merged
                .par_iter_mut()
                .zip(source.par_iter())
                .for_each(|(left, &right)| {
                    let mut byte = 0_u8;
                    for base in 0..4_u8 {
                        let support = lane_support(*left, base)
                            .saturating_add(lane_support(right, base))
                            .min(2);
                        byte |= support << (2 * base);
                    }
                    *left = byte;
                });
        });
    }
    Ok(merged)
}

fn materialize_backbone(mut graph: BackboneGraph, support: &[u8]) -> Result<BackboneGraph> {
    let state_count = graph.state_count();
    let mut offsets = vec![0_u32; state_count + 1];
    let mut targets = Vec::new();
    let mut indegree = vec![0_u8; state_count];
    for source in 0..state_count as u32 {
        let mut local = [u32::MAX; 4];
        let mut count = 0_usize;
        for base in 0..4_u8 {
            if lane_support(support[source as usize], base) == 0 {
                continue;
            }
            if let Some(target) = graph.target_for_base(source, base) {
                local[count] = target;
                count += 1;
            }
        }
        local[..count].sort_unstable();
        let mut previous = u32::MAX;
        for &target in &local[..count] {
            if target == previous {
                continue;
            }
            previous = target;
            targets.push(target);
            indegree[target as usize] = indegree[target as usize].saturating_add(1);
        }
        if targets.len() > u32::MAX as usize {
            bail!("hybrid backbone has too many directed edges");
        }
        offsets[source as usize + 1] = targets.len() as u32;
    }
    graph.out_offsets = offsets;
    graph.out_targets = targets;
    graph.indegree = indegree;
    Ok(graph)
}

fn build_backbone(config: &HybridConfig, pool: &ThreadPool) -> Result<(BackboneGraph, usize)> {
    let (mut repeated, read_pairs) = discover_backbone(config)?;
    let keys = exact_backbone(config, &repeated, pool)?;
    repeated.release();
    if keys.len() > (u32::MAX as usize) / 2 {
        bail!("hybrid backbone has too many retained k-mers for u32 states");
    }
    let state_count = keys.len() * 2;
    let prefix = PrefixIndex::build(&keys, config.backbone_k);
    let skeleton = BackboneGraph {
        k: config.backbone_k,
        keys,
        prefix,
        out_offsets: vec![0; state_count + 1],
        out_targets: Vec::new(),
        indegree: vec![0; state_count],
    };
    let support = count_backbone_edges(config, &skeleton, pool)?;
    Ok((materialize_backbone(skeleton, &support)?, read_pairs))
}

#[derive(Debug)]
struct Neighborhood {
    id: u32,
    seed_state: u32,
    states: Vec<u32>,
    ambiguous_states: usize,
}

fn state_ambiguous(graph: &BackboneGraph, state: u32) -> bool {
    let indegree = graph.indegree[state as usize] as usize;
    let outdegree = graph.outdegree(state);
    indegree > 1 || outdegree > 1
}

fn graph_ambiguous_states(graph: &BackboneGraph) -> usize {
    (0..graph.state_count() as u32)
        .filter(|&state| state_ambiguous(graph, state))
        .count()
}

fn find_neighborhoods(config: &HybridConfig, graph: &BackboneGraph) -> Vec<Neighborhood> {
    let mut branch_seeds = Vec::new();
    let mut dead_seeds = Vec::new();
    for state in 0..graph.state_count() as u32 {
        let indegree = graph.indegree[state as usize] as usize;
        let outdegree = graph.outdegree(state);
        if indegree > 1 || outdegree > 1 {
            branch_seeds.push(state);
        } else if indegree + outdegree > 0 && (indegree == 0 || outdegree == 0) {
            dead_seeds.push(state);
        }
    }
    branch_seeds.sort_unstable_by_key(|&state| {
        let indegree = graph.indegree[state as usize] as usize;
        let outdegree = graph.outdegree(state);
        std::cmp::Reverse(indegree.saturating_sub(1) + outdegree.saturating_sub(1))
    });
    branch_seeds.extend(dead_seeds);

    let mut claimed = FxHashSet::default();
    let mut neighborhoods = Vec::new();
    for seed in branch_seeds {
        if neighborhoods.len() >= config.max_neighborhoods.min(64) {
            break;
        }
        if claimed.contains(&seed) {
            continue;
        }
        let mut seen = FxHashSet::default();
        let mut queue = VecDeque::new();
        seen.insert(seed);
        queue.push_back((seed, 0_usize));
        while let Some((state, depth)) = queue.pop_front() {
            if depth >= config.neighborhood_radius || seen.len() >= config.max_neighborhood_states {
                continue;
            }
            for edge in graph.out_range(state) {
                let next = graph.out_targets[edge];
                if seen.insert(next) {
                    queue.push_back((next, depth + 1));
                }
            }
            let reverse = state ^ 1;
            for edge in graph.out_range(reverse) {
                let previous = graph.out_targets[edge] ^ 1;
                if seen.insert(previous) {
                    queue.push_back((previous, depth + 1));
                }
            }
        }
        let mut states: Vec<u32> = seen.into_iter().collect();
        states.sort_unstable();
        let ambiguous_states = states
            .iter()
            .filter(|&&state| state_ambiguous(graph, state))
            .count();
        claimed.extend(states.iter().copied());
        neighborhoods.push(Neighborhood {
            id: neighborhoods.len() as u32,
            seed_state: seed,
            states,
            ambiguous_states,
        });
    }
    neighborhoods
}

fn build_route_index(
    graph: &BackboneGraph,
    neighborhoods: &[Neighborhood],
) -> FxHashMap<u128, u64> {
    let mut index = FxHashMap::default();
    for neighborhood in neighborhoods {
        let mask = 1_u64 << neighborhood.id;
        for &state in &neighborhood.states {
            let key = graph.keys[(state / 2) as usize];
            *index.entry(key).or_insert(NO_NEIGHBORHOOD) |= mask;
        }
    }
    index
}

fn route_mask(sequence: &[u8], k: usize, index: &FxHashMap<u128, u64>) -> u64 {
    let mut roller = Roller128::new(k);
    let mut mask = 0_u64;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            continue;
        };
        if let Some((key, _)) = roller.push(bits) {
            mask |= index.get(&key).copied().unwrap_or(NO_NEIGHBORHOOD);
        }
    }
    mask
}

fn longest_valid_run(sequence: &[u8]) -> usize {
    let mut best = 0_usize;
    let mut current = 0_usize;
    for &base in sequence {
        if base_bits(base).is_some() {
            current += 1;
            best = best.max(current);
        } else {
            current = 0;
        }
    }
    best
}

#[derive(Debug, Default, Clone)]
struct RouteStats {
    fragments: usize,
    top_spans: Vec<usize>,
}

impl RouteStats {
    fn observe(&mut self, span: usize, support: usize) {
        self.fragments += 1;
        self.top_spans.push(span);
        self.top_spans
            .sort_unstable_by(|left, right| right.cmp(left));
        self.top_spans.truncate(support.max(1));
    }

    fn supported_max_k(&self, support: usize, cap: usize) -> usize {
        if support == 0 || self.top_spans.len() < support {
            return 0;
        }
        odd_floor(self.top_spans[support - 1].min(cap))
    }
}

fn route_reads(
    config: &HybridConfig,
    graph: &BackboneGraph,
    neighborhoods: &[Neighborhood],
    route_index: &FxHashMap<u128, u64>,
) -> Result<(Vec<PathBuf>, Vec<RouteStats>)> {
    let route_dir = config.output_dir.join("hybrid_routes");
    fs::create_dir_all(&route_dir)?;
    let paths: Vec<PathBuf> = neighborhoods
        .iter()
        .map(|item| route_dir.join(format!("neighborhood_{:03}.tsv", item.id + 1)))
        .collect();
    let mut writers: Vec<BufWriter<File>> = paths
        .iter()
        .map(|path| File::create(path).map(BufWriter::new))
        .collect::<std::io::Result<_>>()?;
    let mut stats = vec![RouteStats::default(); neighborhoods.len()];
    for_each_pair(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |index, left, right| {
            let mut mask = route_mask(&left.sequence, graph.k, route_index);
            if let Some(right) = right.as_ref() {
                mask |= route_mask(&right.sequence, graph.k, route_index);
            }
            if mask == 0 {
                return Ok(());
            }
            let right_sequence = right
                .as_ref()
                .map(|record| record.sequence.as_slice())
                .unwrap_or_default();
            let right_quality = right
                .as_ref()
                .map(|record| record.quality.as_slice())
                .unwrap_or_default();
            let span = longest_valid_run(&left.sequence).max(longest_valid_run(right_sequence));
            let mut remaining = mask;
            while remaining != 0 {
                let bit = remaining.trailing_zeros() as usize;
                remaining &= remaining - 1;
                if bit >= writers.len() {
                    continue;
                }
                writeln!(
                    writers[bit],
                    "{}\t{}\t{}\t{}\t{}",
                    index + 1,
                    String::from_utf8_lossy(&left.sequence),
                    String::from_utf8_lossy(&left.quality),
                    String::from_utf8_lossy(right_sequence),
                    String::from_utf8_lossy(right_quality)
                )?;
                stats[bit].observe(span, config.min_fragment_support as usize);
            }
            Ok(())
        },
    )?;
    for writer in &mut writers {
        writer.flush()?;
    }
    Ok((paths, stats))
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash, Ord, PartialOrd)]
struct WideKey {
    words: [u64; WIDE_WORDS],
}

impl WideKey {
    fn shift_left_append(&mut self, bits: u8, k: usize) {
        let mut carry = 0_u64;
        for word in &mut self.words {
            let next = *word >> 62;
            *word = (*word << 2) | carry;
            carry = next;
        }
        self.words[0] |= u64::from(bits & 0b11);
        self.mask_to_k(k);
    }

    fn shift_right_two(&mut self) {
        let mut carry = 0_u64;
        for word in self.words.iter_mut().rev() {
            let next = *word & 0b11;
            *word = (*word >> 2) | (carry << 62);
            carry = next;
        }
    }

    fn get_two_bits(self, offset: usize) -> u8 {
        let word = offset / 64;
        let shift = offset % 64;
        if shift <= 62 {
            ((self.words[word] >> shift) & 0b11) as u8
        } else {
            let low = (self.words[word] >> 63) & 1;
            let high = if word + 1 < WIDE_WORDS {
                (self.words[word + 1] & 1) << 1
            } else {
                0
            };
            (low | high) as u8
        }
    }

    fn set_two_bits(&mut self, offset: usize, bits: u8) {
        let word = offset / 64;
        let shift = offset % 64;
        let bits = u64::from(bits & 0b11);
        if shift <= 62 {
            self.words[word] &= !(0b11_u64 << shift);
            self.words[word] |= bits << shift;
        } else {
            self.words[word] &= !(1_u64 << 63);
            self.words[word] |= (bits & 1) << 63;
            if word + 1 < WIDE_WORDS {
                self.words[word + 1] &= !1_u64;
                self.words[word + 1] |= (bits >> 1) & 1;
            }
        }
    }

    fn mask_to_k(&mut self, k: usize) {
        let bits = 2 * k;
        let full_words = bits / 64;
        let remainder = bits % 64;
        let keep_words = full_words + usize::from(remainder > 0);
        for word in self.words.iter_mut().skip(keep_words) {
            *word = 0;
        }
        if remainder > 0 {
            self.words[full_words] &= (1_u64 << remainder) - 1;
        }
    }

    fn reverse_complement(self, k: usize) -> Self {
        let mut output = Self::default();
        for position in 0..k {
            let bits = self.get_two_bits(2 * position);
            output.shift_left_append(3 - bits, k);
        }
        output
    }

    fn canonical(self, k: usize) -> (Self, bool) {
        let reverse = self.reverse_complement(k);
        if reverse < self {
            (reverse, true)
        } else {
            (self, false)
        }
    }

    fn to_sequence(self, k: usize) -> Vec<u8> {
        let mut sequence = vec![b'A'; k];
        for (position, base) in sequence.iter_mut().enumerate() {
            *base = match self.get_two_bits(2 * (k - 1 - position)) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
        }
        sequence
    }
}

#[derive(Clone, Copy, Debug)]
struct WideRoller {
    k: usize,
    forward: WideKey,
    reverse: WideKey,
    valid: usize,
}

impl WideRoller {
    fn new(k: usize) -> Self {
        Self {
            k,
            forward: WideKey::default(),
            reverse: WideKey::default(),
            valid: 0,
        }
    }

    fn reset(&mut self) {
        self.forward = WideKey::default();
        self.reverse = WideKey::default();
        self.valid = 0;
    }

    fn push(&mut self, bits: u8) -> Option<(WideKey, bool)> {
        self.forward.shift_left_append(bits, self.k);
        self.reverse.shift_right_two();
        self.reverse.set_two_bits(2 * (self.k - 1), 3 - bits);
        self.valid += 1;
        if self.valid < self.k {
            None
        } else if self.reverse < self.forward {
            Some((self.reverse, true))
        } else {
            Some((self.forward, false))
        }
    }
}

fn read_routed<F>(path: &Path, mut callback: F) -> Result<()>
where
    F: FnMut(u64, &[u8], &[u8], &[u8], &[u8]) -> Result<()>,
{
    let reader = BufReader::new(File::open(path)?);
    for line in reader.lines() {
        let line = line?;
        let mut fields = line.splitn(5, '\t');
        let fragment_id = fields
            .next()
            .context("missing routed fragment id")?
            .parse::<u64>()?;
        let left_sequence = fields.next().context("missing routed left sequence")?;
        let left_quality = fields.next().context("missing routed left quality")?;
        let right_sequence = fields.next().context("missing routed right sequence")?;
        let right_quality = fields.next().context("missing routed right quality")?;
        callback(
            fragment_id,
            left_sequence.as_bytes(),
            left_quality.as_bytes(),
            right_sequence.as_bytes(),
            right_quality.as_bytes(),
        )?;
    }
    Ok(())
}

fn scan_wide_exact(
    sequence: &[u8],
    quality: &[u8],
    k: usize,
    fragment_id: u64,
    evidence: &mut FxHashMap<WideKey, LocalEvidence>,
) {
    let mut roller = WideRoller::new(k);
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
        let entry = evidence.entry(key).or_default();
        entry.exact.count = entry.exact.count.saturating_add(1);
        entry.exact.quality_sum = entry.exact.quality_sum.saturating_add(quality_sum);
        if entry.last_fragment != fragment_id {
            entry.last_fragment = fragment_id;
            entry.exact.fragment_count = entry.exact.fragment_count.saturating_add(1);
        }
    }
}

fn wide_state_first_base_bits(key: WideKey, reverse: bool, k: usize) -> u8 {
    if reverse {
        3 - key.get_two_bits(0)
    } else {
        key.get_two_bits(2 * (k - 1))
    }
}

fn wide_target_for_base(
    keys: &[WideKey],
    index: &FxHashMap<WideKey, u32>,
    k: usize,
    source: u32,
    base: u8,
) -> Option<u32> {
    let key = keys[(source / 2) as usize];
    let mut oriented = if source & 1 == 0 {
        key
    } else {
        key.reverse_complement(k)
    };
    oriented.shift_left_append(base, k);
    let (canonical, reverse) = oriented.canonical(k);
    index
        .get(&canonical)
        .copied()
        .map(|node| node * 2 + u32::from(reverse))
}

fn scan_wide_edges(
    sequence: &[u8],
    k: usize,
    keys: &[WideKey],
    index: &FxHashMap<WideKey, u32>,
    support: &mut [u8],
) {
    let mut roller = WideRoller::new(k);
    let mut previous: Option<u32> = None;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous = None;
            continue;
        };
        let Some((key, reverse)) = roller.push(bits) else {
            continue;
        };
        let current = index
            .get(&key)
            .copied()
            .map(|node| node * 2 + u32::from(reverse));
        if let (Some(source), Some(target)) = (previous, current) {
            bump_lane(&mut support[source as usize], bits);
            let source_key = keys[(source / 2) as usize];
            let reverse_base = 3 - wide_state_first_base_bits(source_key, source & 1 != 0, k);
            bump_lane(&mut support[(target ^ 1) as usize], reverse_base);
        }
        previous = current;
    }
}

fn project_wide_key(key: WideKey, k: usize, backbone: &BackboneGraph) -> (usize, usize) {
    let sequence = key.to_sequence(k);
    let mut roller = Roller128::new(backbone.k);
    let mut previous = None;
    for &base in &sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous = None;
            continue;
        };
        let Some((low_key, reverse)) = roller.push(bits) else {
            continue;
        };
        let Some(node) = backbone.find_node(low_key) else {
            return (1, 0);
        };
        let state = node * 2 + u32::from(reverse);
        if let Some(previous_state) = previous {
            if !backbone.has_edge(previous_state, state) {
                return (0, 1);
            }
        }
        previous = Some(state);
    }
    (0, 0)
}

fn build_local_k(
    config: &HybridConfig,
    path: &Path,
    k: usize,
    backbone: &BackboneGraph,
) -> Result<LocalKSummary> {
    let mut evidence: FxHashMap<WideKey, LocalEvidence> = FxHashMap::default();
    read_routed(path, |fragment_id, left, left_q, right, right_q| {
        scan_wide_exact(left, left_q, k, fragment_id, &mut evidence);
        if !right.is_empty() {
            scan_wide_exact(right, right_q, k, fragment_id, &mut evidence);
        }
        Ok(())
    })?;
    let mut keys: Vec<WideKey> = evidence
        .into_iter()
        .filter_map(|(key, value)| {
            (value.exact.count >= config.min_count
                && value.exact.fragment_count >= config.min_fragment_support
                && value.exact.mean_quality(k) >= config.min_mean_quality)
                .then_some(key)
        })
        .collect();
    keys.sort_unstable();
    let index: FxHashMap<WideKey, u32> = keys
        .iter()
        .copied()
        .enumerate()
        .map(|(id, key)| (key, id as u32))
        .collect();
    let mut support = vec![0_u8; keys.len() * 2];
    read_routed(path, |_fragment_id, left, _left_q, right, _right_q| {
        scan_wide_edges(left, k, &keys, &index, &mut support);
        if !right.is_empty() {
            scan_wide_edges(right, k, &keys, &index, &mut support);
        }
        Ok(())
    })?;

    let mut indegree = vec![0_u8; keys.len() * 2];
    let mut outdegree = vec![0_u8; keys.len() * 2];
    let mut directed_edges = 0_usize;
    for source in 0..support.len() as u32 {
        let mut targets = [u32::MAX; 4];
        let mut count = 0_usize;
        for base in 0..4_u8 {
            if lane_support(support[source as usize], base) == 0 {
                continue;
            }
            if let Some(target) = wide_target_for_base(&keys, &index, k, source, base) {
                targets[count] = target;
                count += 1;
            }
        }
        targets[..count].sort_unstable();
        let mut previous = u32::MAX;
        for &target in &targets[..count] {
            if target == previous {
                continue;
            }
            previous = target;
            directed_edges += 1;
            outdegree[source as usize] = outdegree[source as usize].saturating_add(1);
            indegree[target as usize] = indegree[target as usize].saturating_add(1);
        }
    }
    let ambiguous_states = indegree
        .iter()
        .zip(outdegree.iter())
        .filter(|&(&incoming, &outgoing)| incoming > 1 || outgoing > 1)
        .count();
    let dead_end_states = indegree
        .iter()
        .zip(outdegree.iter())
        .filter(|&(&incoming, &outgoing)| {
            incoming + outgoing > 0 && (incoming == 0 || outgoing == 0)
        })
        .count();
    let (projection_missing_nodes, projection_missing_edges) = keys
        .par_iter()
        .map(|&key| project_wide_key(key, k, backbone))
        .reduce(
            || (0_usize, 0_usize),
            |left, right| (left.0 + right.0, left.1 + right.1),
        );
    let resolved = !keys.is_empty()
        && ambiguous_states == 0
        && projection_missing_nodes == 0
        && projection_missing_edges == 0;
    Ok(LocalKSummary {
        k,
        retained_kmers: keys.len(),
        directed_edges,
        ambiguous_states,
        dead_end_states,
        projection_missing_nodes,
        projection_missing_edges,
        resolved,
    })
}

fn odd_floor(value: usize) -> usize {
    if value.is_multiple_of(2) {
        value.saturating_sub(1)
    } else {
        value
    }
}

fn odd_ceil(value: usize) -> usize {
    if value.is_multiple_of(2) {
        value.saturating_add(1)
    } else {
        value
    }
}

fn adaptive_k_schedule(start: usize, max_k: usize) -> Vec<usize> {
    let max_k = odd_floor(max_k.min(WIDE_MAX_K));
    let mut current = odd_ceil(start).min(max_k);
    if current == 0 || current > max_k {
        return Vec::new();
    }
    let mut result = Vec::new();
    loop {
        result.push(current);
        if current >= max_k {
            break;
        }
        let growth = (current / 3).max(24);
        let mut next = odd_ceil((current + growth).min(max_k));
        if next > max_k {
            next = max_k;
        }
        if next <= current {
            break;
        }
        current = next;
    }
    if result.last().copied() != Some(max_k) && max_k >= start {
        result.push(max_k);
    }
    result
}

fn refine_one(
    config: &HybridConfig,
    neighborhood: &Neighborhood,
    path: &Path,
    stats: &RouteStats,
    backbone: &BackboneGraph,
) -> Result<HybridNeighborhoodSummary> {
    let supported_max_k = stats.supported_max_k(
        config.min_fragment_support as usize,
        config.max_local_k.min(WIDE_MAX_K),
    );
    let start = config.local_start_k.max(config.backbone_k + 2);
    let mut candidates = Vec::new();
    let mut selected_k = None;
    let mut resolved = false;
    if stats.fragments >= config.min_fragment_support as usize && supported_max_k >= start {
        for k in adaptive_k_schedule(start, supported_max_k) {
            let result = build_local_k(config, path, k, backbone)?;
            selected_k = Some(k);
            resolved = result.resolved;
            candidates.push(result);
            if resolved {
                break;
            }
        }
    }
    Ok(HybridNeighborhoodSummary {
        id: neighborhood.id + 1,
        seed_state: neighborhood.seed_state,
        states: neighborhood.states.len(),
        ambiguous_states: neighborhood.ambiguous_states,
        routed_fragments: stats.fragments,
        supported_max_k,
        selected_k,
        resolved,
        candidates,
    })
}

pub fn run_multik_hybrid(config: &HybridConfig, threads: usize) -> Result<HybridRunStats> {
    if config.backbone_k == 0 || config.backbone_k > BACKBONE_MAX_K {
        bail!("hybrid backbone k must be in 1..={BACKBONE_MAX_K}");
    }
    if config.min_count < 2 || config.min_fragment_support < 2 {
        bail!("hybrid Stage34 requires min-count>=2 and min-fragment-support>=2");
    }
    if config.max_local_k > WIDE_MAX_K {
        bail!("hybrid local k cannot exceed {WIDE_MAX_K}");
    }
    if config.max_neighborhoods == 0 || config.max_neighborhoods > 64 {
        bail!("hybrid max-neighborhoods must be in 1..=64");
    }
    if config.max_neighborhood_states == 0 {
        bail!("hybrid max-neighborhood-states must be positive");
    }
    if !(0.0..=60.0).contains(&config.min_mean_quality) {
        bail!("hybrid minimum mean quality must be in 0..=60");
    }
    fs::create_dir_all(&config.output_dir)?;
    let threads = threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let total_started = Instant::now();

    let backbone_started = Instant::now();
    let (backbone, read_pairs) = build_backbone(config, &pool)?;
    let backbone_seconds = backbone_started.elapsed().as_secs_f64();
    let backbone_ambiguous = graph_ambiguous_states(&backbone);
    let neighborhoods = find_neighborhoods(config, &backbone);
    eprintln!(
        "stage34 v5 hybrid backbone k={} complete: {} retained nodes, {} directed edges, {} ambiguous states, {} bounded neighborhoods in {:.3}s",
        config.backbone_k,
        backbone.keys.len(),
        backbone.out_targets.len(),
        backbone_ambiguous,
        neighborhoods.len(),
        backbone_seconds
    );

    let routing_started = Instant::now();
    let route_index = build_route_index(&backbone, &neighborhoods);
    let (route_paths, route_stats) = route_reads(config, &backbone, &neighborhoods, &route_index)?;
    let routing_seconds = routing_started.elapsed().as_secs_f64();
    let routed_fragment_copies: usize = route_stats.iter().map(|stats| stats.fragments).sum();
    eprintln!(
        "stage34 v5 conservative routing complete: {} routed fragment copies in {:.3}s",
        routed_fragment_copies, routing_seconds
    );

    let local_started = Instant::now();
    let local_results: Vec<Result<HybridNeighborhoodSummary>> = pool.install(|| {
        neighborhoods
            .par_iter()
            .zip(route_paths.par_iter())
            .zip(route_stats.par_iter())
            .map(|((neighborhood, path), stats)| {
                refine_one(config, neighborhood, path, stats, &backbone)
            })
            .collect()
    });
    let mut local = Vec::with_capacity(local_results.len());
    for result in local_results {
        local.push(result?);
    }
    local.sort_unstable_by_key(|item| item.id);
    let local_refinement_seconds = local_started.elapsed().as_secs_f64();
    let resolved = local.iter().filter(|item| item.resolved).count();
    let above_55 = local
        .iter()
        .filter_map(|item| item.selected_k)
        .filter(|&k| k > 55)
        .count();
    eprintln!(
        "stage34 v5 local refinement complete: {resolved}/{} neighborhoods resolved, {above_55} selected k>55 in {:.3}s",
        local.len(), local_refinement_seconds
    );

    for path in &route_paths {
        let _ = fs::remove_file(path);
    }
    let _ = fs::remove_dir(config.output_dir.join("hybrid_routes"));

    let summary = HybridSummary {
        version: "stage34-hybrid-adaptive-local-v5".to_string(),
        read_pairs,
        threads,
        backbone_k: config.backbone_k,
        backbone_retained_kmers: backbone.keys.len(),
        backbone_directed_edges: backbone.out_targets.len(),
        backbone_ambiguous_states: backbone_ambiguous,
        neighborhoods: neighborhoods.len(),
        routed_fragment_copies,
        local_start_k: config.local_start_k,
        max_local_k: config.max_local_k.min(WIDE_MAX_K),
        backbone_seconds,
        routing_seconds,
        local_refinement_seconds,
        total_seconds: total_started.elapsed().as_secs_f64(),
        local,
    };
    Ok(HybridRunStats { summary })
}

pub fn write_hybrid_outputs(run: &HybridRunStats, output: &Path) -> Result<()> {
    fs::create_dir_all(output)?;
    fs::write(
        output.join("hybrid_summary.json"),
        serde_json::to_vec_pretty(&run.summary)?,
    )?;
    let mut table = BufWriter::new(File::create(output.join("hybrid_neighborhoods.tsv"))?);
    writeln!(
        table,
        "id\tseed_state\tstates\tambiguous_states\trouted_fragments\tsupported_max_k\tselected_k\tresolved"
    )?;
    for item in &run.summary.local {
        writeln!(
            table,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            item.id,
            item.seed_state,
            item.states,
            item.ambiguous_states,
            item.routed_fragments,
            item.supported_max_k,
            item.selected_k
                .map_or_else(|| "NA".to_string(), |k| k.to_string()),
            item.resolved
        )?;
    }
    table.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wide_from_sequence(sequence: &[u8]) -> WideKey {
        let mut key = WideKey::default();
        for &base in sequence {
            key.shift_left_append(base_bits(base).unwrap(), sequence.len());
        }
        key
    }

    #[test]
    fn wide_key_roundtrips_past_u128_limit() {
        let sequence: Vec<u8> = (0..191)
            .map(|index| match index % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let key = wide_from_sequence(&sequence);
        assert_eq!(key.to_sequence(sequence.len()), sequence);
        let rc = key.reverse_complement(sequence.len());
        assert_eq!(rc.reverse_complement(sequence.len()), key);
    }

    #[test]
    fn supported_max_k_uses_independent_fragment_rank() {
        let mut stats = RouteStats::default();
        for span in [250_usize, 203, 202, 180] {
            stats.observe(span, 2);
        }
        assert_eq!(stats.supported_max_k(2, 255), 203);
        assert_eq!(stats.supported_max_k(3, 255), 0);
    }

    #[test]
    fn adaptive_schedule_reaches_read_supported_high_k() {
        let schedule = adaptive_k_schedule(55, 203);
        assert_eq!(schedule.first().copied(), Some(55));
        assert_eq!(schedule.last().copied(), Some(203));
        assert!(schedule.iter().any(|&k| k > 55));
        assert!(schedule.iter().all(|k| k % 2 == 1));
    }

    #[test]
    fn route_mask_is_conservative_for_any_matching_backbone_kmer() {
        let k = 5;
        let sequence = b"ACGTTGCA";
        let mut keys = Vec::new();
        collect_u128_keys(sequence, k, &mut keys);
        let mut index = FxHashMap::default();
        index.insert(keys[1], 1_u64 << 3);
        assert_eq!(route_mask(sequence, k, &index), 1_u64 << 3);
    }

    #[test]
    fn wide_roller_matches_direct_canonicalization() {
        let sequence: Vec<u8> = (0..220)
            .map(|index| match (index * 7) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let k = 191;
        let mut roller = WideRoller::new(k);
        for (end, &base) in sequence.iter().enumerate() {
            let observed = roller.push(base_bits(base).unwrap());
            if end + 1 < k {
                assert!(observed.is_none());
                continue;
            }
            let direct = wide_from_sequence(&sequence[end + 1 - k..=end]);
            assert_eq!(observed.unwrap(), direct.canonical(k));
        }
    }
}
