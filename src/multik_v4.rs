use crate::dna::base_bits;
use crate::fastq::{for_each_pair, FastqRecord};
use crate::graph::GraphSummary;
use crate::multik::{
    AdaptiveRescueCandidate, AdaptiveRescueSummary, MultiKConfig, MultiKLayerSummary,
    MultiKSummary, MultiKTimingSummary, ProjectionSummary,
};
use crate::multik_stream::Stage34RunStats;
use anyhow::{bail, Result};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;
use std::time::{Duration, Instant};

const MAX_COMPACT_K: usize = 63;
const BLOOM_BITS: usize = 1 << 30;
const BLOOM_WORDS: usize = BLOOM_BITS / 64;
const HLL_P: u32 = 16;
const HLL_REGISTERS: usize = 1 << HLL_P;
const BATCH_PAIRS: usize = 8192;
const PROGRESS_PAIRS: usize = 1_000_000;
const PREFIX_BITS: usize = 23;
const EXACT_GROUP_SIZE: usize = 2;
const NO_UNITIG: u32 = u32::MAX;
const PAIR_EVEN_MASK: u128 = 0x5555_5555_5555_5555_5555_5555_5555_5555;
const PAIR_ODD_MASK: u128 = 0xaaaa_aaaa_aaaa_aaaa_aaaa_aaaa_aaaa_aaaa;

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
    last_fragment: u32,
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

#[derive(Debug)]
struct PairLite {
    fragment_id: u32,
    left_sequence: Vec<u8>,
    left_quality: Vec<u8>,
    right_sequence: Option<Vec<u8>>,
    right_quality: Option<Vec<u8>>,
}

impl PairLite {
    fn new(index: usize, left: FastqRecord, right: Option<FastqRecord>) -> Result<Self> {
        let fragment_id = u32::try_from(index + 1)
            .map_err(|_| anyhow::anyhow!("too many read pairs for u32 fragment ids"))?;
        let (right_sequence, right_quality) = match right {
            Some(record) => (Some(record.sequence), Some(record.quality)),
            None => (None, None),
        };
        Ok(Self {
            fragment_id,
            left_sequence: left.sequence,
            left_quality: left.quality,
            right_sequence,
            right_quality,
        })
    }
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
struct PreparedLayer {
    k: usize,
    keys: Vec<u128>,
    prefix_index: PrefixIndex,
    node_hll: Hll,
    edge_hll: Hll,
    observations: u64,
    edge_observations: u64,
}

impl PreparedLayer {
    #[inline]
    fn state_count(&self) -> usize {
        self.keys.len() * 2
    }

    #[inline]
    fn find_node(&self, key: u128) -> Option<u32> {
        self.prefix_index.find(&self.keys, key)
    }

    #[inline]
    fn state_first_base_bits(&self, state: u32) -> u8 {
        let key = self.keys[(state / 2) as usize];
        if state & 1 == 0 {
            ((key >> (2 * (self.k - 1))) & 0b11) as u8
        } else {
            3 - (key & 0b11) as u8
        }
    }

    #[inline]
    fn target_for_base(&self, source: u32, base: u8) -> Option<u32> {
        let key = self.keys[(source / 2) as usize];
        let key_rc = reverse_complement_packed(key, self.k);
        let (oriented, reverse_oriented) = if source & 1 == 0 {
            (key, key_rc)
        } else {
            (key_rc, key)
        };
        let next = ((oriented << 2) | u128::from(base)) & packed_mask(self.k);
        let next_rc =
            (u128::from(3 - (base & 0b11)) << (2 * (self.k - 1))) | (reverse_oriented >> 2);
        let (canonical, reverse) = if next_rc < next {
            (next_rc, true)
        } else {
            (next, false)
        };
        self.find_node(canonical)
            .map(|node| node * 2 + u32::from(reverse))
    }
}

#[derive(Debug)]
struct CompactRawGraph {
    k: usize,
    keys: Vec<u128>,
    prefix_index: PrefixIndex,
    out_offsets: Vec<u32>,
    out_targets: Vec<u32>,
    indegree: Vec<u8>,
    singleton_solid_edges: usize,
}

impl CompactRawGraph {
    #[inline]
    fn state_count(&self) -> usize {
        self.keys.len() * 2
    }

    #[inline]
    fn out_range(&self, state: u32) -> std::ops::Range<usize> {
        let state = state as usize;
        self.out_offsets[state] as usize..self.out_offsets[state + 1] as usize
    }

    #[inline]
    fn outdegree(&self, state: u32) -> usize {
        self.out_range(state).len()
    }

    #[inline]
    fn find_node(&self, key: u128) -> Option<u32> {
        self.prefix_index.find(&self.keys, key)
    }

    #[inline]
    fn reverse_state(state: u32) -> u32 {
        state ^ 1
    }

    #[inline]
    fn state_last_base(&self, state: u32) -> u8 {
        let key = self.keys[(state / 2) as usize];
        let bits = if state & 1 == 0 {
            (key & 0b11) as u8
        } else {
            3 - ((key >> (2 * (self.k - 1))) & 0b11) as u8
        };
        bits_to_base(bits)
    }

    fn state_sequence(&self, state: u32) -> Vec<u8> {
        let node = (state / 2) as usize;
        let key = if state & 1 == 0 {
            self.keys[node]
        } else {
            reverse_complement_packed(self.keys[node], self.k)
        };
        packed_to_sequence(key, self.k)
    }

    fn edge_index(&self, source: u32, target: u32) -> Option<usize> {
        let range = self.out_range(source);
        self.out_targets[range.clone()]
            .binary_search(&target)
            .ok()
            .map(|local| range.start + local)
    }
}

#[derive(Debug)]
struct CompactUnitig {
    id: u32,
    states: Vec<u32>,
    sequence: Vec<u8>,
    start_state: u32,
    end_state: u32,
    length: usize,
}

#[derive(Debug)]
struct CompactUnitigGraph {
    unitigs: Vec<CompactUnitig>,
    edge_unitig: Vec<u32>,
    out_offsets: Vec<u32>,
    out_targets: Vec<u32>,
    indegree: Vec<u32>,
}

impl CompactUnitigGraph {
    #[inline]
    fn out_range(&self, unitig: u32) -> std::ops::Range<usize> {
        let unitig = unitig as usize;
        self.out_offsets[unitig] as usize..self.out_offsets[unitig + 1] as usize
    }

    #[inline]
    fn outdegree(&self, unitig: u32) -> usize {
        self.out_range(unitig).len()
    }

    fn unitig_for_edge(&self, raw: &CompactRawGraph, source: u32, target: u32) -> Option<u32> {
        let edge_index = raw.edge_index(source, target)?;
        let unitig = self.edge_unitig[edge_index];
        (unitig != NO_UNITIG).then_some(unitig)
    }
}

#[derive(Debug)]
struct CompactLayer {
    k: usize,
    raw: CompactRawGraph,
    unitigs: CompactUnitigGraph,
    summary: MultiKLayerSummary,
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

#[inline]
fn reverse_complement_packed(value: u128, k: usize) -> u128 {
    let reversed = value.reverse_bits();
    let swapped = ((reversed & PAIR_ODD_MASK) >> 1) | ((reversed & PAIR_EVEN_MASK) << 1);
    let aligned = swapped >> (128 - 2 * k);
    aligned ^ packed_mask(k)
}

#[inline]
fn bits_to_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

fn packed_to_sequence(mut value: u128, k: usize) -> Vec<u8> {
    let mut sequence = vec![b'A'; k];
    for base in sequence.iter_mut().rev() {
        *base = bits_to_base((value & 0b11) as u8);
        value >>= 2;
    }
    sequence
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
        batch.push(PairLite::new(index, left, right)?);
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

fn scan_discovery_sequence(
    sequence: &[u8],
    layer: &mut DiscoveryLayer,
    fragment_keys: &mut Vec<u128>,
) {
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

fn discover_repeated_parallel(
    config: &MultiKConfig,
    ks: &[usize],
    pool: &ThreadPool,
) -> Result<(Vec<DiscoveryLayer>, usize)> {
    let mut layers: Vec<DiscoveryLayer> = ks.iter().copied().map(DiscoveryLayer::new).collect();
    let mut processed = 0_usize;
    let read_pairs = process_pair_batches(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |batch| {
            pool.install(|| {
                layers.par_iter_mut().for_each(|layer| {
                    let mut fragment_keys = Vec::with_capacity(512);
                    for pair in batch {
                        fragment_keys.clear();
                        scan_discovery_sequence(&pair.left_sequence, layer, &mut fragment_keys);
                        if let Some(sequence) = pair.right_sequence.as_deref() {
                            scan_discovery_sequence(sequence, layer, &mut fragment_keys);
                        }
                        fragment_keys.sort_unstable();
                        fragment_keys.dedup();
                        for &key in &fragment_keys {
                            if layer.seen.contains(key) {
                                layer.repeated.insert(key);
                            } else {
                                layer.seen.insert(key);
                            }
                        }
                    }
                });
            });
            processed += batch.len();
            if processed / PROGRESS_PAIRS != (processed - batch.len()) / PROGRESS_PAIRS {
                eprintln!("stage34 v4 discovery progress: {processed} read pairs");
            }
            Ok(())
        },
    )?;
    for layer in &mut layers {
        layer.seen.release();
    }
    Ok((layers, read_pairs))
}

fn scan_exact_sequence(
    sequence: &[u8],
    quality: &[u8],
    k: usize,
    repeated: &Bloom,
    fragment_id: u32,
    evidence: &mut FxHashMap<u128, LocalEvidence>,
) {
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
        let entry = left.entry(key).or_default();
        entry.exact.add(value.exact);
    }
    left
}

fn exact_group(
    config: &MultiKConfig,
    layers: &[DiscoveryLayer],
    pool: &ThreadPool,
) -> Result<Vec<FxHashMap<u128, ExactEvidence>>> {
    let mut globals: Vec<FxHashMap<u128, ExactEvidence>> =
        layers.iter().map(|_| FxHashMap::default()).collect();
    let mut processed = 0_usize;
    process_pair_batches(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |batch| {
            let chunk = batch
                .len()
                .div_ceil(pool.current_num_threads().max(1))
                .max(1);
            let batch_maps: Vec<FxHashMap<u128, LocalEvidence>> = pool.install(|| {
                layers
                    .par_iter()
                    .map(|layer| {
                        batch
                            .par_chunks(chunk)
                            .map(|pairs| {
                                let mut local = FxHashMap::default();
                                for pair in pairs {
                                    scan_exact_sequence(
                                        &pair.left_sequence,
                                        &pair.left_quality,
                                        layer.k,
                                        &layer.repeated,
                                        pair.fragment_id,
                                        &mut local,
                                    );
                                    if let (Some(sequence), Some(quality)) = (
                                        pair.right_sequence.as_deref(),
                                        pair.right_quality.as_deref(),
                                    ) {
                                        scan_exact_sequence(
                                            sequence,
                                            quality,
                                            layer.k,
                                            &layer.repeated,
                                            pair.fragment_id,
                                            &mut local,
                                        );
                                    }
                                }
                                local
                            })
                            .reduce(FxHashMap::default, merge_local_maps)
                    })
                    .collect()
            });
            for (global, local) in globals.iter_mut().zip(batch_maps) {
                for (key, value) in local {
                    global.entry(key).or_default().add(value.exact);
                }
            }
            processed += batch.len();
            if processed / PROGRESS_PAIRS != (processed - batch.len()) / PROGRESS_PAIRS {
                let ks: Vec<_> = layers.iter().map(|layer| layer.k).collect();
                eprintln!("stage34 v4 exact k={ks:?} progress: {processed} read pairs");
            }
            Ok(())
        },
    )?;
    Ok(globals)
}

fn prepare_layer(
    mut discovery: DiscoveryLayer,
    evidence: FxHashMap<u128, ExactEvidence>,
    config: &MultiKConfig,
) -> Result<PreparedLayer> {
    let k = discovery.k;
    discovery.repeated.release();
    let mut keys: Vec<u128> = evidence
        .iter()
        .filter_map(|(key, value)| {
            (value.count >= config.min_count
                && value.fragment_count >= config.min_fragment_support
                && value.mean_quality(k) >= config.min_mean_quality)
                .then_some(*key)
        })
        .collect();
    drop(evidence);
    keys.sort_unstable();
    if keys.len() > (u32::MAX as usize) / 2 {
        bail!("k={k} graph has too many canonical nodes");
    }
    let prefix_index = PrefixIndex::build(&keys, k);
    eprintln!("stage34 v4 k={k}: retained {} nodes", keys.len());
    Ok(PreparedLayer {
        k,
        keys,
        prefix_index,
        node_hll: discovery.node_hll,
        edge_hll: discovery.edge_hll,
        observations: discovery.observations,
        edge_observations: discovery.edge_observations,
    })
}

#[inline]
fn lane_support(byte: u8, base: u8) -> u8 {
    (byte >> (2 * (base & 0b11))) & 0b11
}

#[inline]
fn bump_lane(byte: &mut u8, base: u8) {
    let shift = 2 * (base & 0b11);
    let value = (*byte >> shift) & 0b11;
    if value < 2 {
        *byte = (*byte & !(0b11 << shift)) | ((value + 1) << shift);
    }
}

#[inline]
fn merge_support_byte(left: u8, right: u8) -> u8 {
    let mut merged = 0_u8;
    for base in 0..4_u8 {
        let value = lane_support(left, base)
            .saturating_add(lane_support(right, base))
            .min(2);
        merged |= value << (2 * base);
    }
    merged
}

fn scan_edge_sequence(sequence: &[u8], layer: &PreparedLayer, support: &mut [u8]) {
    let mut roller = Roller::new(layer.k);
    let mut previous_state: Option<u32> = None;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous_state = None;
            continue;
        };
        let Some((key, reverse)) = roller.push(bits) else {
            continue;
        };
        let current_state = layer
            .find_node(key)
            .map(|node| node * 2 + u32::from(reverse));
        if let (Some(source), Some(target)) = (previous_state, current_state) {
            bump_lane(&mut support[source as usize], bits);
            let reverse_source = CompactRawGraph::reverse_state(target);
            let reverse_base = 3 - layer.state_first_base_bits(source);
            bump_lane(&mut support[reverse_source as usize], reverse_base);
        }
        previous_state = current_state;
    }
}

#[derive(Debug)]
struct EdgeWorker {
    supports: Vec<Vec<u8>>,
}

fn count_edges_all_k(
    config: &MultiKConfig,
    layers: &[PreparedLayer],
    pool: &ThreadPool,
) -> Result<Vec<Vec<u8>>> {
    let workers = pool.current_num_threads().max(1);
    let mut worker_state: Vec<EdgeWorker> = (0..workers)
        .map(|_| EdgeWorker {
            supports: layers
                .iter()
                .map(|layer| vec![0_u8; layer.state_count()])
                .collect(),
        })
        .collect();
    let mut processed = 0_usize;
    process_pair_batches(
        &config.read1,
        config.read2.as_deref(),
        config.max_pairs,
        |batch| {
            let batch_len = batch.len();
            pool.install(|| {
                worker_state
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(worker_index, worker)| {
                        let start = batch_len * worker_index / workers;
                        let end = batch_len * (worker_index + 1) / workers;
                        for pair in &batch[start..end] {
                            for (layer, support) in layers.iter().zip(worker.supports.iter_mut()) {
                                scan_edge_sequence(&pair.left_sequence, layer, support);
                                if let Some(sequence) = pair.right_sequence.as_deref() {
                                    scan_edge_sequence(sequence, layer, support);
                                }
                            }
                        }
                    });
            });
            processed += batch.len();
            if processed / PROGRESS_PAIRS != (processed - batch.len()) / PROGRESS_PAIRS {
                eprintln!("stage34 v4 all-k edge progress: {processed} read pairs");
            }
            Ok(())
        },
    )?;

    let mut workers_iter = worker_state.into_iter();
    let mut merged = workers_iter
        .next()
        .map(|worker| worker.supports)
        .unwrap_or_else(|| layers.iter().map(|layer| vec![0; layer.state_count()]).collect());
    for worker in workers_iter {
        for (target, source) in merged.iter_mut().zip(worker.supports) {
            pool.install(|| {
                target
                    .par_iter_mut()
                    .zip(source.par_iter())
                    .for_each(|(left, &right)| *left = merge_support_byte(*left, right));
            });
        }
    }
    Ok(merged)
}

fn build_raw_from_support(prepared: PreparedLayer, support: Vec<u8>) -> Result<CompactRawGraph> {
    let PreparedLayer {
        k,
        keys,
        prefix_index,
        ..
    } = prepared;
    if support.len() != keys.len() * 2 {
        bail!("k={k} packed edge support length mismatch");
    }
    let lookup = PreparedLayer {
        k,
        keys,
        prefix_index,
        node_hll: Hll::new(),
        edge_hll: Hll::new(),
        observations: 0,
        edge_observations: 0,
    };
    let state_count = lookup.state_count();
    let mut out_offsets = vec![0_u32; state_count + 1];
    let mut out_targets = Vec::with_capacity(state_count);
    let mut indegree = vec![0_u8; state_count];
    let mut singleton_solid_edges = 0_usize;

    for source in 0..state_count as u32 {
        let mut entries = [(0_u32, 0_u8); 4];
        let mut entry_count = 0_usize;
        let packed = support[source as usize];
        for base in 0..4_u8 {
            let lane = lane_support(packed, base);
            if lane == 0 {
                continue;
            }
            let target = lookup.target_for_base(source, base).ok_or_else(|| {
                anyhow::anyhow!("k={k} retained edge target missing for state {source} base {base}")
            })?;
            entries[entry_count] = (target, lane);
            entry_count += 1;
        }
        entries[..entry_count].sort_unstable_by_key(|&(target, _)| target);
        let mut index = 0_usize;
        while index < entry_count {
            let target = entries[index].0;
            let mut lane = entries[index].1;
            index += 1;
            while index < entry_count && entries[index].0 == target {
                lane = lane.saturating_add(entries[index].1).min(2);
                index += 1;
            }
            if lane < 2 {
                singleton_solid_edges += 1;
            }
            out_targets.push(target);
            indegree[target as usize] = indegree[target as usize].saturating_add(1);
        }
        if out_targets.len() > u32::MAX as usize {
            bail!("k={k} graph has too many directed edges for compact CSR");
        }
        out_offsets[source as usize + 1] = out_targets.len() as u32;
    }

    let PreparedLayer {
        keys,
        prefix_index,
        ..
    } = lookup;
    Ok(CompactRawGraph {
        k,
        keys,
        prefix_index,
        out_offsets,
        out_targets,
        indegree,
        singleton_solid_edges,
    })
}

fn compact_unitigs(raw: &CompactRawGraph) -> CompactUnitigGraph {
    let mut visited = vec![false; raw.out_targets.len()];
    let mut edge_unitig = vec![NO_UNITIG; raw.out_targets.len()];
    let mut unitigs = Vec::new();

    for source in 0..raw.state_count() as u32 {
        if raw.outdegree(source) == 0 {
            continue;
        }
        if raw.indegree[source as usize] == 1 && raw.outdegree(source) == 1 {
            continue;
        }
        for edge_index in raw.out_range(source) {
            if visited[edge_index] {
                continue;
            }
            let id = unitigs.len() as u32;
            let target = raw.out_targets[edge_index];
            let states = walk_unitig(
                raw,
                source,
                target,
                edge_index,
                id,
                &mut visited,
                &mut edge_unitig,
            );
            push_unitig(raw, states, id, &mut unitigs);
        }
    }

    for source in 0..raw.state_count() as u32 {
        for edge_index in raw.out_range(source) {
            if visited[edge_index] {
                continue;
            }
            let id = unitigs.len() as u32;
            let target = raw.out_targets[edge_index];
            let states = walk_unitig(
                raw,
                source,
                target,
                edge_index,
                id,
                &mut visited,
                &mut edge_unitig,
            );
            push_unitig(raw, states, id, &mut unitigs);
        }
    }

    for node in 0..raw.keys.len() as u32 {
        let state = node * 2;
        if raw.indegree[state as usize] == 0 && raw.outdegree(state) == 0 {
            let id = unitigs.len() as u32;
            push_unitig(raw, vec![state], id, &mut unitigs);
        }
    }

    let (out_offsets, out_targets, indegree) =
        build_unitig_adjacency(raw, &unitigs, &edge_unitig);
    CompactUnitigGraph {
        unitigs,
        edge_unitig,
        out_offsets,
        out_targets,
        indegree,
    }
}

fn walk_unitig(
    raw: &CompactRawGraph,
    source: u32,
    target: u32,
    first_edge: usize,
    unitig_id: u32,
    visited: &mut [bool],
    edge_unitig: &mut [u32],
) -> Vec<u32> {
    let mut states = vec![source, target];
    visited[first_edge] = true;
    edge_unitig[first_edge] = unitig_id;
    let mut current = target;
    let max_steps = raw.out_targets.len().saturating_add(1);
    for _ in 0..max_steps {
        if raw.indegree[current as usize] != 1 || raw.outdegree(current) != 1 {
            break;
        }
        let edge_index = raw.out_range(current).start;
        if visited[edge_index] {
            break;
        }
        let next = raw.out_targets[edge_index];
        visited[edge_index] = true;
        edge_unitig[edge_index] = unitig_id;
        states.push(next);
        current = next;
    }
    states
}

fn push_unitig(raw: &CompactRawGraph, states: Vec<u32>, id: u32, unitigs: &mut Vec<CompactUnitig>) {
    if states.is_empty() {
        return;
    }
    let start_state = states[0];
    let end_state = *states.last().expect("nonempty unitig states");
    let mut sequence = raw.state_sequence(start_state);
    sequence.reserve(states.len().saturating_sub(1));
    for &state in states.iter().skip(1) {
        sequence.push(raw.state_last_base(state));
    }
    let length = sequence.len();
    unitigs.push(CompactUnitig {
        id,
        states,
        sequence,
        start_state,
        end_state,
        length,
    });
}

fn build_unitig_adjacency(
    raw: &CompactRawGraph,
    unitigs: &[CompactUnitig],
    edge_unitig: &[u32],
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut out_offsets = vec![0_u32; unitigs.len() + 1];
    let mut out_targets = Vec::new();
    let mut indegree = vec![0_u32; unitigs.len()];
    for unitig in unitigs {
        let mut targets = [NO_UNITIG; 4];
        let mut count = 0_usize;
        for edge_index in raw.out_range(unitig.end_state) {
            let target = edge_unitig[edge_index];
            if target == NO_UNITIG || target == unitig.id {
                continue;
            }
            targets[count] = target;
            count += 1;
        }
        targets[..count].sort_unstable();
        let mut previous = NO_UNITIG;
        for &target in &targets[..count] {
            if target == previous {
                continue;
            }
            previous = target;
            out_targets.push(target);
            indegree[target as usize] = indegree[target as usize].saturating_add(1);
        }
        out_offsets[unitig.id as usize + 1] = out_targets.len() as u32;
    }
    (out_offsets, out_targets, indegree)
}

fn summarize(raw: &CompactRawGraph, unitigs: &CompactUnitigGraph) -> GraphSummary {
    let mut lengths: Vec<usize> = unitigs.unitigs.iter().map(|unitig| unitig.length).collect();
    let unitig_bases = lengths.iter().sum();
    let largest_unitig = lengths.iter().copied().max().unwrap_or(0);
    let unitig_n50 = n50(&mut lengths);
    let branching_unitigs = unitigs
        .unitigs
        .iter()
        .filter(|unitig| {
            unitigs.indegree[unitig.id as usize] > 1 || unitigs.outdegree(unitig.id) > 1
        })
        .count();
    GraphSummary {
        canonical_nodes: raw.keys.len(),
        oriented_states: raw.state_count(),
        directed_edges: raw.out_targets.len(),
        singleton_solid_edges: raw.singleton_solid_edges,
        unitigs: unitigs.unitigs.len(),
        unitig_edges: unitigs.out_targets.len(),
        branching_unitigs,
        unitig_bases,
        unitig_n50,
        largest_unitig,
    }
}

fn n50(lengths: &mut [usize]) -> usize {
    if lengths.is_empty() {
        return 0;
    }
    lengths.sort_unstable_by(|left, right| right.cmp(left));
    let total: usize = lengths.iter().sum();
    let mut cumulative = 0_usize;
    for &length in lengths.iter() {
        cumulative += length;
        if cumulative * 2 >= total {
            return length;
        }
    }
    0
}

fn build_compact_layer(prepared: PreparedLayer, support: Vec<u8>) -> Result<CompactLayer> {
    let k = prepared.k;
    let node_hll_estimate = prepared.node_hll.estimate();
    let edge_hll_estimate = prepared.edge_hll.estimate();
    let observations = prepared.observations;
    let edge_observations = prepared.edge_observations;
    let raw = build_raw_from_support(prepared, support)?;
    let unitigs = compact_unitigs(&raw);
    let graph_summary = summarize(&raw, &unitigs);
    let summary = MultiKLayerSummary {
        k,
        observations,
        edge_observations,
        distinct_kmers: node_hll_estimate,
        retained_kmers: raw.keys.len(),
        distinct_kplus1: edge_hll_estimate,
        retained_directed_edges: raw.out_targets.len(),
        graph: graph_summary,
    };
    Ok(CompactLayer {
        k,
        raw,
        unitigs,
        summary,
    })
}

#[derive(Clone, Copy, Debug, Default)]
struct ProjectionAccumulator {
    exact: usize,
    missing_node: usize,
    missing_edge: usize,
    visits: usize,
}

fn project_one(high_sequence: &[u8], low: &CompactLayer) -> ProjectionAccumulator {
    let mut roller = Roller::new(low.k);
    let mut previous_state: Option<u32> = None;
    let mut last_unitig: Option<u32> = None;
    let mut visits = 0_usize;
    for &base in high_sequence {
        let Some(bits) = base_bits(base) else {
            roller.reset();
            previous_state = None;
            last_unitig = None;
            continue;
        };
        let Some((key, reverse)) = roller.push(bits) else {
            continue;
        };
        let Some(node) = low.raw.find_node(key) else {
            return ProjectionAccumulator {
                missing_node: 1,
                ..ProjectionAccumulator::default()
            };
        };
        let state = node * 2 + u32::from(reverse);
        if let Some(previous) = previous_state {
            let Some(unitig) = low.unitigs.unitig_for_edge(&low.raw, previous, state) else {
                return ProjectionAccumulator {
                    missing_edge: 1,
                    ..ProjectionAccumulator::default()
                };
            };
            if last_unitig != Some(unitig) {
                visits += 1;
                last_unitig = Some(unitig);
            }
        }
        previous_state = Some(state);
    }
    ProjectionAccumulator {
        exact: 1,
        visits,
        ..ProjectionAccumulator::default()
    }
}

fn build_projection_summary(
    high: &CompactLayer,
    low: &CompactLayer,
    pool: &ThreadPool,
) -> ProjectionSummary {
    let total = pool.install(|| {
        high.unitigs
            .unitigs
            .par_iter()
            .map(|unitig| project_one(&unitig.sequence, low))
            .reduce(ProjectionAccumulator::default, |left, right| {
                ProjectionAccumulator {
                    exact: left.exact + right.exact,
                    missing_node: left.missing_node + right.missing_node,
                    missing_edge: left.missing_edge + right.missing_edge,
                    visits: left.visits + right.visits,
                }
            })
    });
    ProjectionSummary {
        high_k: high.k,
        low_k: low.k,
        high_unitigs: high.unitigs.unitigs.len(),
        exact_unitigs: total.exact,
        missing_node_unitigs: total.missing_node,
        missing_edge_unitigs: total.missing_edge,
        projected_low_unitig_visits: total.visits,
    }
}

#[derive(Debug, Default)]
struct RescueAccumulator {
    dead_end_unitigs: usize,
    low_anchor_missing: usize,
    branch_stops: usize,
    dead_stops: usize,
    cycle_stops: usize,
    maxed_out: usize,
    reanchored: usize,
    total_extension_bases: usize,
    max_extension_bases: usize,
    candidates: Vec<AdaptiveRescueCandidate>,
}

impl RescueAccumulator {
    fn merge(mut self, mut other: Self) -> Self {
        self.dead_end_unitigs += other.dead_end_unitigs;
        self.low_anchor_missing += other.low_anchor_missing;
        self.branch_stops += other.branch_stops;
        self.dead_stops += other.dead_stops;
        self.cycle_stops += other.cycle_stops;
        self.maxed_out += other.maxed_out;
        self.reanchored += other.reanchored;
        self.total_extension_bases += other.total_extension_bases;
        self.max_extension_bases = self.max_extension_bases.max(other.max_extension_bases);
        self.candidates.append(&mut other.candidates);
        self
    }
}

fn rescue_one(
    unitig: &CompactUnitig,
    high: &CompactLayer,
    low: &CompactLayer,
    max_bases: usize,
) -> RescueAccumulator {
    let mut result = RescueAccumulator::default();
    if high.unitigs.outdegree(unitig.id) != 0 || unitig.sequence.len() < high.k {
        return result;
    }
    result.dead_end_unitigs = 1;

    let anchor_start = unitig.sequence.len().saturating_sub(low.k);
    let mut low_roller = Roller::new(low.k);
    let mut anchor = None;
    for &base in &unitig.sequence[anchor_start..] {
        let Some(bits) = base_bits(base) else {
            low_roller.reset();
            continue;
        };
        anchor = low_roller.push(bits);
    }
    let Some((anchor_key, anchor_reverse)) = anchor else {
        result.low_anchor_missing = 1;
        return result;
    };
    let Some(anchor_node) = low.raw.find_node(anchor_key) else {
        result.low_anchor_missing = 1;
        return result;
    };

    let mut current = anchor_node * 2 + u32::from(anchor_reverse);
    let mut visited = FxHashSet::default();
    visited.insert(current);
    let context_start = unitig.sequence.len().saturating_sub(high.k - 1);
    let mut high_roller = Roller::new(high.k);
    for &base in &unitig.sequence[context_start..] {
        if let Some(bits) = base_bits(base) {
            let _ = high_roller.push(bits);
        }
    }
    let mut extension_bases = 0_usize;
    while extension_bases < max_bases {
        let range = low.raw.out_range(current);
        if range.is_empty() {
            result.dead_stops = 1;
            return result;
        }
        if range.len() != 1 {
            result.branch_stops = 1;
            return result;
        }
        let next = low.raw.out_targets[range.start];
        if !visited.insert(next) {
            result.cycle_stops = 1;
            return result;
        }
        let base = low.raw.state_last_base(next);
        extension_bases += 1;
        current = next;
        let Some(bits) = base_bits(base) else {
            continue;
        };
        let Some((key, reverse)) = high_roller.push(bits) else {
            continue;
        };
        if let Some(node) = high.raw.find_node(key) {
            let high_state = node * 2 + u32::from(reverse);
            if !unitig.states.contains(&high_state) {
                result.reanchored = 1;
                result.total_extension_bases = extension_bases;
                result.max_extension_bases = extension_bases;
                result.candidates.push(AdaptiveRescueCandidate {
                    high_k: high.k,
                    low_k: low.k,
                    high_unitig: unitig.id,
                    reanchor_state: high_state,
                    extension_bases,
                });
                return result;
            }
        }
    }
    result.maxed_out = 1;
    result
}

fn probe_adaptive_rescue(
    high: &CompactLayer,
    low: &CompactLayer,
    max_bases: usize,
    pool: &ThreadPool,
) -> (AdaptiveRescueSummary, Vec<AdaptiveRescueCandidate>) {
    let result = pool.install(|| {
        high.unitigs
            .unitigs
            .par_iter()
            .map(|unitig| rescue_one(unitig, high, low, max_bases))
            .reduce(RescueAccumulator::default, RescueAccumulator::merge)
    });
    let summary = AdaptiveRescueSummary {
        high_k: high.k,
        low_k: low.k,
        dead_end_unitigs: result.dead_end_unitigs,
        low_anchor_missing: result.low_anchor_missing,
        branch_stops: result.branch_stops,
        dead_stops: result.dead_stops,
        cycle_stops: result.cycle_stops,
        maxed_out: result.maxed_out,
        reanchored: result.reanchored,
        total_extension_bases: result.total_extension_bases,
        max_extension_bases: result.max_extension_bases,
    };
    (summary, result.candidates)
}

pub fn run_multik_v4(config: &MultiKConfig, threads: usize) -> Result<Stage34RunStats> {
    if config.min_count <= 1 || config.min_fragment_support <= 1 {
        bail!("Stage34 v4 requires min-count>=2 and min-fragment-support>=2");
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

    let threads = threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
    eprintln!("stage34 v4 using {threads} worker threads");
    let total_started = Instant::now();

    let discovery_started = Instant::now();
    let (mut discovery, read_pairs) = discover_repeated_parallel(config, &ks, &pool)?;
    let discovery_seconds = discovery_started.elapsed().as_secs_f64();
    eprintln!(
        "stage34 v4 repeat discovery complete: {read_pairs} pairs in {discovery_seconds:.3}s"
    );

    let exact_started = Instant::now();
    let mut prepared = Vec::with_capacity(discovery.len());
    while !discovery.is_empty() {
        let group_len = EXACT_GROUP_SIZE.min(discovery.len());
        let group: Vec<DiscoveryLayer> = discovery.drain(..group_len).collect();
        let group_ks: Vec<_> = group.iter().map(|layer| layer.k).collect();
        let started = Instant::now();
        let evidences = exact_group(config, &group, &pool)?;
        eprintln!(
            "stage34 v4 exact group {group_ks:?} complete in {:.3}s",
            started.elapsed().as_secs_f64()
        );
        for (layer, evidence) in group.into_iter().zip(evidences) {
            prepared.push(prepare_layer(layer, evidence, config)?);
        }
    }
    let exact_seconds = exact_started.elapsed().as_secs_f64();

    let edge_started = Instant::now();
    let supports = count_edges_all_k(config, &prepared, &pool)?;
    let edge_seconds = edge_started.elapsed().as_secs_f64();
    eprintln!("stage34 v4 all-k edge pass complete in {edge_seconds:.3}s");

    let graph_started = Instant::now();
    let mut projection_duration = Duration::ZERO;
    let mut rescue_duration = Duration::ZERO;
    let mut layer_summaries = Vec::with_capacity(prepared.len());
    let mut projection_summaries = Vec::with_capacity(prepared.len().saturating_sub(1));
    let mut rescue_summaries = Vec::with_capacity(prepared.len().saturating_sub(1));
    let mut rescue_candidates = Vec::new();
    let mut previous: Option<CompactLayer> = None;

    for (prepared_layer, support) in prepared.into_iter().zip(supports) {
        let k = prepared_layer.k;
        let started = Instant::now();
        let current = build_compact_layer(prepared_layer, support)?;
        eprintln!(
            "stage34 v4 k={k} graph complete in {:.3}s; {} nodes, {} directed edges, {} unitigs",
            started.elapsed().as_secs_f64(),
            current.summary.retained_kmers,
            current.summary.retained_directed_edges,
            current.summary.graph.unitigs
        );
        layer_summaries.push(current.summary.clone());

        if let Some(low) = previous.take() {
            let started = Instant::now();
            projection_summaries.push(build_projection_summary(&current, &low, &pool));
            projection_duration += started.elapsed();

            let started = Instant::now();
            let (rescue, mut candidates) =
                probe_adaptive_rescue(&current, &low, config.max_rescue_bases, &pool);
            rescue_duration += started.elapsed();
            rescue_summaries.push(rescue);
            rescue_candidates.append(&mut candidates);
            drop(low);
        }
        previous = Some(current);
    }
    drop(previous);
    let graph_seconds = graph_started.elapsed().as_secs_f64()
        - projection_duration.as_secs_f64()
        - rescue_duration.as_secs_f64();
    let finalize_seconds = exact_seconds + edge_seconds + graph_seconds;
    eprintln!(
        "stage34 v4 timings: discovery={discovery_seconds:.3}s exact={exact_seconds:.3}s edge={edge_seconds:.3}s graph={graph_seconds:.3}s projection={:.3}s rescue={:.3}s",
        projection_duration.as_secs_f64(),
        rescue_duration.as_secs_f64()
    );

    let summary = MultiKSummary {
        version: "stage34-layered-multik-v4-four-pass-packed-edge".to_string(),
        read_pairs,
        ks,
        min_count: config.min_count,
        min_fragment_support: config.min_fragment_support,
        min_mean_quality: config.min_mean_quality,
        max_rescue_bases: config.max_rescue_bases,
        layers: layer_summaries,
        projections: projection_summaries,
        adaptive_rescues: rescue_summaries,
        timings_seconds: MultiKTimingSummary {
            count_seconds: discovery_seconds,
            finalize_seconds,
            projection_seconds: projection_duration.as_secs_f64(),
            rescue_seconds: rescue_duration.as_secs_f64(),
            total_seconds: total_started.elapsed().as_secs_f64(),
        },
    };
    Ok(Stage34RunStats {
        summary,
        rescue_candidates,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multik_stream::run_multik_streaming;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn reverse_complement_fast_matches_reference() {
        fn reference(mut value: u128, k: usize) -> u128 {
            let mut output = 0_u128;
            for _ in 0..k {
                let bits = (value & 0b11) as u8;
                value >>= 2;
                output = (output << 2) | u128::from(3 - bits);
            }
            output
        }
        for k in [1_usize, 5, 21, 31, 55, 62] {
            let mask = packed_mask(k);
            for seed in 0..128_u128 {
                let value = mix64(seed as u64) as u128 & mask;
                assert_eq!(reverse_complement_packed(value, k), reference(value, k));
            }
        }
    }

    #[test]
    fn support_lanes_saturate_independently() {
        let mut byte = 0_u8;
        bump_lane(&mut byte, 0);
        bump_lane(&mut byte, 0);
        bump_lane(&mut byte, 0);
        bump_lane(&mut byte, 3);
        assert_eq!(lane_support(byte, 0), 2);
        assert_eq!(lane_support(byte, 1), 0);
        assert_eq!(lane_support(byte, 3), 1);
        assert_eq!(merge_support_byte(byte, byte), 0b10_00_00_10);
    }

    #[test]
    fn v4_matches_v3_fixture_summaries() {
        let dir = tempfile::tempdir().unwrap();
        let reads = dir.path().join("reads.fastq");
        let sequence = "ACGTTGCAACGTCAGTACGATCGTAGCTAACGTTGCA";
        let mut handle = File::create(&reads).unwrap();
        for index in 0..8 {
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
        let v3 = run_multik_streaming(&config, 2).unwrap();
        let v4 = run_multik_v4(&config, 2).unwrap();
        assert_eq!(
            serde_json::to_value(&v4.summary.layers).unwrap(),
            serde_json::to_value(&v3.summary.layers).unwrap()
        );
        assert_eq!(
            serde_json::to_value(&v4.summary.projections).unwrap(),
            serde_json::to_value(&v3.summary.projections).unwrap()
        );
        assert_eq!(
            serde_json::to_value(&v4.summary.adaptive_rescues).unwrap(),
            serde_json::to_value(&v3.summary.adaptive_rescues).unwrap()
        );
        assert_eq!(
            v4.summary.version,
            "stage34-layered-multik-v4-four-pass-packed-edge"
        );
    }
}
