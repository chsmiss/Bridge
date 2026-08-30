use crate::dna::{canonical_kmers, reverse_complement, KmerKey, OrientedKmer};
use crate::fastq::for_each_pair;
use crate::kmer::KmerSet;
use anyhow::{bail, Result};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::path::Path;

#[derive(Debug)]
pub struct RawGraph {
    pub k: usize,
    pub keys: Vec<KmerKey>,
    pub counts: Vec<u32>,
    pub out_offsets: Vec<u64>,
    pub out_targets: Vec<u32>,
    pub indegree: Vec<u32>,
}

impl RawGraph {
    #[inline]
    pub fn state_count(&self) -> usize {
        self.keys.len() * 2
    }

    #[inline]
    pub fn out_range(&self, state: u32) -> std::ops::Range<usize> {
        let state = state as usize;
        self.out_offsets[state] as usize..self.out_offsets[state + 1] as usize
    }

    #[inline]
    pub fn outdegree(&self, state: u32) -> usize {
        self.out_range(state).len()
    }

    #[inline]
    pub fn reverse_state(state: u32) -> u32 {
        state ^ 1
    }

    pub fn state_sequence(&self, state: u32) -> Vec<u8> {
        let node = (state / 2) as usize;
        if state & 1 == 0 {
            self.keys[node].to_sequence(self.k)
        } else {
            self.keys[node]
                .reverse_complement(self.k)
                .to_sequence(self.k)
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Unitig {
    pub id: u32,
    #[serde(skip_serializing)]
    pub states: Vec<u32>,
    #[serde(skip_serializing)]
    pub sequence: Vec<u8>,
    pub start_state: u32,
    pub end_state: u32,
    pub length: usize,
    pub mean_coverage: f32,
    pub min_coverage: u32,
    pub max_coverage: u32,
    pub circular: bool,
}

#[derive(Debug)]
pub struct UnitigGraph {
    pub k: usize,
    pub unitigs: Vec<Unitig>,
    pub edge_to_unitig: FxHashMap<(u32, u32), u32>,
    pub reverse_unitig: Vec<u32>,
    pub out_offsets: Vec<u64>,
    pub out_targets: Vec<u32>,
    pub indegree: Vec<u32>,
}

impl UnitigGraph {
    #[inline]
    pub fn out_range(&self, unitig: u32) -> std::ops::Range<usize> {
        let unitig = unitig as usize;
        self.out_offsets[unitig] as usize..self.out_offsets[unitig + 1] as usize
    }

    #[inline]
    pub fn outdegree(&self, unitig: u32) -> usize {
        self.out_range(unitig).len()
    }

    pub fn incoming_at_state(&self, state: u32) -> Vec<u32> {
        self.unitigs
            .iter()
            .filter_map(|unitig| (unitig.end_state == state).then_some(unitig.id))
            .collect()
    }

    pub fn outgoing_at_state(&self, state: u32) -> Vec<u32> {
        self.unitigs
            .iter()
            .filter_map(|unitig| (unitig.start_state == state).then_some(unitig.id))
            .collect()
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct GraphSummary {
    pub canonical_nodes: usize,
    pub oriented_states: usize,
    pub directed_edges: usize,
    pub unitigs: usize,
    pub unitig_edges: usize,
    pub branching_unitigs: usize,
    pub unitig_bases: usize,
    pub unitig_n50: usize,
    pub largest_unitig: usize,
}

pub fn build_raw_graph(
    read1: &Path,
    read2: Option<&Path>,
    kmer_set: &KmerSet,
    max_pairs: Option<usize>,
) -> Result<RawGraph> {
    let mut keys: Vec<KmerKey> = kmer_set.retained.iter().copied().collect();
    keys.sort_unstable();
    if keys.len() > (u32::MAX as usize) / 2 {
        bail!("graph has too many canonical nodes for u32 oriented IDs");
    }
    let index: FxHashMap<KmerKey, u32> = keys
        .iter()
        .enumerate()
        .map(|(node_id, key)| (*key, node_id as u32))
        .collect();
    let counts: Vec<u32> = keys
        .iter()
        .map(|key| kmer_set.evidence.get(key).map_or(0, |value| value.count))
        .collect();

    // Count (k+1)-mer transitions rather than materializing every observed
    // adjacency. This follows the edge-centric logic used by succinct DBG
    // assemblers and removes many one-off error branches that survive node-only
    // abundance filtering. Mercy-rescued paths remain eligible at count one.
    let mut edge_counts: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    for_each_pair(read1, read2, max_pairs, |_pair_index, left, right| {
        add_record_edges(
            &left.sequence,
            kmer_set.summary.k,
            &index,
            &mut edge_counts,
        )?;
        if let Some(right) = right {
            add_record_edges(
                &right.sequence,
                kmer_set.summary.k,
                &index,
                &mut edge_counts,
            )?;
        }
        Ok(())
    })?;

    let min_edge_count = kmer_set.summary.min_count.max(1);
    let mut edges: Vec<(u32, u32)> = edge_counts
        .into_iter()
        .filter_map(|((source, target), support)| {
            let source_key = keys[(source / 2) as usize];
            let target_key = keys[(target / 2) as usize];
            let mercy_edge = kmer_set.rescued.contains(&source_key)
                || kmer_set.rescued.contains(&target_key);
            (support >= min_edge_count || mercy_edge).then_some((source, target))
        })
        .collect();
    edges.sort_unstable();

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

    Ok(RawGraph {
        k: kmer_set.summary.k,
        keys,
        counts,
        out_offsets,
        out_targets,
        indegree,
    })
}

fn add_record_edges(
    sequence: &[u8],
    k: usize,
    index: &FxHashMap<KmerKey, u32>,
    edge_counts: &mut FxHashMap<(u32, u32), u32>,
) -> Result<()> {
    let kmers = canonical_kmers(sequence, k)?;
    for pair in kmers.windows(2) {
        let left = pair[0];
        let right = pair[1];
        if right.position != left.position + 1 {
            continue;
        }
        let (Some(&left_id), Some(&right_id)) = (index.get(&left.key), index.get(&right.key))
        else {
            continue;
        };
        let left_state = oriented_state(left_id, left);
        let right_state = oriented_state(right_id, right);
        let forward = edge_counts.entry((left_state, right_state)).or_insert(0);
        *forward = forward.saturating_add(1);
        let reverse = edge_counts
            .entry((
                RawGraph::reverse_state(right_state),
                RawGraph::reverse_state(left_state),
            ))
            .or_insert(0);
        *reverse = reverse.saturating_add(1);
    }
    Ok(())
}

#[inline]
fn oriented_state(node_id: u32, kmer: OrientedKmer) -> u32 {
    node_id * 2 + u32::from(kmer.reverse)
}

pub fn compact_unitigs(graph: &RawGraph) -> UnitigGraph {
    let mut visited = vec![false; graph.out_targets.len()];
    let mut unitigs = Vec::new();
    let mut edge_to_unitig: FxHashMap<(u32, u32), u32> = FxHashMap::default();

    for source in 0..graph.state_count() as u32 {
        if graph.outdegree(source) == 0 {
            continue;
        }
        if graph.indegree[source as usize] == 1 && graph.outdegree(source) == 1 {
            continue;
        }
        for edge_index in graph.out_range(source) {
            if visited[edge_index] {
                continue;
            }
            let target = graph.out_targets[edge_index];
            let states = walk_unitig(graph, source, target, edge_index, &mut visited);
            push_unitig(graph, states, false, &mut unitigs, &mut edge_to_unitig);
        }
    }

    for source in 0..graph.state_count() as u32 {
        for edge_index in graph.out_range(source) {
            if visited[edge_index] {
                continue;
            }
            let target = graph.out_targets[edge_index];
            let states = walk_unitig(graph, source, target, edge_index, &mut visited);
            push_unitig(graph, states, true, &mut unitigs, &mut edge_to_unitig);
        }
    }

    for node in 0..graph.keys.len() as u32 {
        let state = node * 2;
        if graph.indegree[state as usize] == 0 && graph.outdegree(state) == 0 {
            push_unitig(graph, vec![state], false, &mut unitigs, &mut edge_to_unitig);
        }
    }

    let reverse_unitig = compute_reverse_unitigs(&unitigs);
    let (out_offsets, out_targets, indegree) = build_unitig_adjacency(&unitigs);
    UnitigGraph {
        k: graph.k,
        unitigs,
        edge_to_unitig,
        reverse_unitig,
        out_offsets,
        out_targets,
        indegree,
    }
}

fn walk_unitig(
    graph: &RawGraph,
    source: u32,
    target: u32,
    first_edge: usize,
    visited: &mut [bool],
) -> Vec<u32> {
    let mut states = vec![source, target];
    visited[first_edge] = true;
    let mut current = target;
    let max_steps = graph.out_targets.len().saturating_add(1);

    for _ in 0..max_steps {
        if graph.indegree[current as usize] != 1 || graph.outdegree(current) != 1 {
            break;
        }
        let edge_index = graph.out_range(current).start;
        if visited[edge_index] {
            break;
        }
        let next = graph.out_targets[edge_index];
        visited[edge_index] = true;
        states.push(next);
        current = next;
    }
    states
}

fn push_unitig(
    graph: &RawGraph,
    states: Vec<u32>,
    circular: bool,
    unitigs: &mut Vec<Unitig>,
    edge_to_unitig: &mut FxHashMap<(u32, u32), u32>,
) {
    if states.is_empty() {
        return;
    }
    let id = unitigs.len() as u32;
    let start_state = states[0];
    let end_state = *states.last().expect("nonempty unitig states");
    let mut sequence = graph.state_sequence(start_state);
    for &state in states.iter().skip(1) {
        let state_sequence = graph.state_sequence(state);
        if let Some(&base) = state_sequence.last() {
            sequence.push(base);
        }
    }

    let mut coverage_sum = 0_u64;
    let mut min_coverage = u32::MAX;
    let mut max_coverage = 0_u32;
    for &state in &states {
        let coverage = graph.counts[(state / 2) as usize];
        coverage_sum += u64::from(coverage);
        min_coverage = min_coverage.min(coverage);
        max_coverage = max_coverage.max(coverage);
    }
    if min_coverage == u32::MAX {
        min_coverage = 0;
    }

    for pair in states.windows(2) {
        edge_to_unitig.insert((pair[0], pair[1]), id);
    }
    let length = sequence.len();
    let state_count = states.len().max(1);
    unitigs.push(Unitig {
        id,
        states,
        sequence,
        start_state,
        end_state,
        length,
        mean_coverage: coverage_sum as f32 / state_count as f32,
        min_coverage,
        max_coverage,
        circular,
    });
}

fn build_unitig_adjacency(unitigs: &[Unitig]) -> (Vec<u64>, Vec<u32>, Vec<u32>) {
    let mut starts: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for unitig in unitigs {
        starts
            .entry(unitig.start_state)
            .or_default()
            .push(unitig.id);
    }

    let mut edges = Vec::new();
    for unitig in unitigs {
        if let Some(targets) = starts.get(&unitig.end_state) {
            for &target in targets {
                if target != unitig.id {
                    edges.push((unitig.id, target));
                }
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();

    let mut out_offsets = vec![0_u64; unitigs.len() + 1];
    let mut indegree = vec![0_u32; unitigs.len()];
    for &(source, target) in &edges {
        out_offsets[source as usize + 1] += 1;
        indegree[target as usize] = indegree[target as usize].saturating_add(1);
    }
    for index in 1..out_offsets.len() {
        out_offsets[index] += out_offsets[index - 1];
    }
    let out_targets = edges.into_iter().map(|(_, target)| target).collect();
    (out_offsets, out_targets, indegree)
}

fn compute_reverse_unitigs(unitigs: &[Unitig]) -> Vec<u32> {
    let mut index: FxHashMap<Vec<u8>, u32> = FxHashMap::default();
    for unitig in unitigs {
        index.insert(unitig.sequence.clone(), unitig.id);
    }
    unitigs
        .iter()
        .map(|unitig| {
            let reverse = reverse_complement(&unitig.sequence);
            index.get(&reverse).copied().unwrap_or(unitig.id)
        })
        .collect()
}

pub fn summarize(graph: &RawGraph, unitigs: &UnitigGraph) -> GraphSummary {
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
        canonical_nodes: graph.keys.len(),
        oriented_states: graph.state_count(),
        directed_edges: graph.out_targets.len(),
        unitigs: unitigs.unitigs.len(),
        unitig_edges: unitigs.out_targets.len(),
        branching_unitigs,
        unitig_bases,
        unitig_n50,
        largest_unitig,
    }
}

pub fn n50(lengths: &mut [usize]) -> usize {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calculates_n50() {
        let mut lengths = vec![10, 30, 20, 40];
        assert_eq!(n50(&mut lengths), 30);
    }
}
