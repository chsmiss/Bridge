use crate::dna::{canonical_kmers, reverse_complement};
use crate::fastq::for_each_pair;
use crate::graph::{
    build_raw_graph, compact_unitigs, summarize, GraphSummary, RawGraph, UnitigGraph,
};
use crate::kmer::{count_and_filter, KmerCountSummary};
use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct AssembleConfig {
    pub read1: PathBuf,
    pub read2: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub k: usize,
    pub min_count: u32,
    pub mercy_max_kmers: usize,
    pub mercy_min_support: u16,
    pub min_read_support: u32,
    pub min_pair_support: u32,
    pub min_contig_length: usize,
    pub max_pairs: Option<usize>,
    pub threads: usize,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub read_pairs: u32,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssemblyStats {
    pub version: String,
    pub kmer: KmerCountSummary,
    pub graph: GraphSummary,
    pub threaded_reads: u64,
    pub threaded_pairs: u64,
    pub direct_transitions: usize,
    pub pair_bridges: usize,
    pub primary_contigs: usize,
    pub primary_bases: usize,
    pub primary_n50: usize,
    pub largest_primary: usize,
    pub timings_seconds: FxHashMap<String, f64>,
}

#[derive(Debug)]
pub struct AssemblyProduct {
    pub raw_graph: RawGraph,
    pub unitig_graph: UnitigGraph,
    pub transitions: FxHashMap<(u32, u32), TransitionEvidence>,
    pub primary_paths: Vec<Vec<u32>>,
    pub primary_sequences: Vec<Vec<u8>>,
    pub stats: AssemblyStats,
}

#[derive(Debug)]
struct ThreadingResult {
    transitions: FxHashMap<(u32, u32), TransitionEvidence>,
    threaded_reads: u64,
    threaded_pairs: u64,
}

pub fn assemble(config: &AssembleConfig) -> Result<AssemblyProduct> {
    fs::create_dir_all(&config.output_dir).with_context(|| {
        format!(
            "failed to create output directory {}",
            config.output_dir.display()
        )
    })?;
    if config.threads == 0 {
        anyhow::bail!("threads must be positive");
    }
    let _ = ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global();

    let total_started = Instant::now();
    let mut timings = FxHashMap::default();

    let started = Instant::now();
    let kmer_set = count_and_filter(
        &config.read1,
        config.read2.as_deref(),
        config.k,
        config.min_count,
        config.mercy_max_kmers,
        config.mercy_min_support,
        config.max_pairs,
    )?;
    timings.insert(
        "kmer_count_filter".to_string(),
        started.elapsed().as_secs_f64(),
    );

    let started = Instant::now();
    let raw_graph = build_raw_graph(
        &config.read1,
        config.read2.as_deref(),
        &kmer_set,
        config.max_pairs,
    )?;
    timings.insert("graph_build".to_string(), started.elapsed().as_secs_f64());

    let started = Instant::now();
    let unitig_graph = compact_unitigs(&raw_graph);
    let graph_summary = summarize(&raw_graph, &unitig_graph);
    timings.insert(
        "unitig_compaction".to_string(),
        started.elapsed().as_secs_f64(),
    );

    let started = Instant::now();
    let ThreadingResult {
        transitions,
        threaded_reads,
        threaded_pairs,
    } = thread_reads(
        &config.read1,
        config.read2.as_deref(),
        &raw_graph,
        &unitig_graph,
        config.max_pairs,
    )?;
    timings.insert(
        "read_pair_threading".to_string(),
        started.elapsed().as_secs_f64(),
    );

    let started = Instant::now();
    let primary_paths = safe_primary_paths(&unitig_graph, &transitions, config.min_read_support);
    let mut primary_sequences =
        deduplicate_primary_sequences(&unitig_graph, &primary_paths, config.min_contig_length);
    primary_sequences
        .sort_unstable_by(|left, right| right.len().cmp(&left.len()).then(left.cmp(right)));
    timings.insert(
        "safe_walk_emission".to_string(),
        started.elapsed().as_secs_f64(),
    );

    let mut primary_lengths: Vec<usize> = primary_sequences.iter().map(Vec::len).collect();
    let primary_bases = primary_lengths.iter().sum();
    let largest_primary = primary_lengths.iter().copied().max().unwrap_or(0);
    let primary_n50 = crate::graph::n50(&mut primary_lengths);
    let direct_transitions = transitions
        .values()
        .filter(|evidence| evidence.direct_reads > 0)
        .count();
    let pair_bridges = transitions
        .values()
        .filter(|evidence| evidence.read_pairs >= config.min_pair_support)
        .count();
    timings.insert("total".to_string(), total_started.elapsed().as_secs_f64());

    let stats = AssemblyStats {
        version: env!("CARGO_PKG_VERSION").to_string(),
        kmer: kmer_set.summary,
        graph: graph_summary,
        threaded_reads,
        threaded_pairs,
        direct_transitions,
        pair_bridges,
        primary_contigs: primary_sequences.len(),
        primary_bases,
        primary_n50,
        largest_primary,
        timings_seconds: timings,
    };

    Ok(AssemblyProduct {
        raw_graph,
        unitig_graph,
        transitions,
        primary_paths,
        primary_sequences,
        stats,
    })
}

fn thread_reads(
    read1: &Path,
    read2: Option<&Path>,
    graph: &RawGraph,
    unitigs: &UnitigGraph,
    max_pairs: Option<usize>,
) -> Result<ThreadingResult> {
    let node_index: FxHashMap<_, _> = graph
        .keys
        .iter()
        .enumerate()
        .map(|(node_id, key)| (*key, node_id as u32))
        .collect();
    let mut transitions: FxHashMap<(u32, u32), TransitionEvidence> = FxHashMap::default();
    let mut threaded_reads = 0_u64;
    let mut threaded_pairs = 0_u64;

    for_each_pair(read1, read2, max_pairs, |_pair_index, left, right| {
        let left_path = thread_record(&left.sequence, graph.k, &node_index, unitigs)?;
        if !left_path.is_empty() {
            threaded_reads += 1;
            add_direct_transitions(&left_path, &mut transitions);
        }

        if let Some(right) = right {
            let right_path = thread_record(&right.sequence, graph.k, &node_index, unitigs)?;
            if !right_path.is_empty() {
                threaded_reads += 1;
                add_direct_transitions(&right_path, &mut transitions);
            }
            if !left_path.is_empty() && !right_path.is_empty() {
                let molecular_right = reverse_unitig_path(&right_path, unitigs);
                if let (Some(&left_end), Some(&right_start)) =
                    (left_path.last(), molecular_right.first())
                {
                    if left_end != right_start {
                        let evidence = transitions.entry((left_end, right_start)).or_default();
                        evidence.read_pairs = evidence.read_pairs.saturating_add(1);
                        threaded_pairs += 1;
                    }
                }
            }
        }
        Ok(())
    })?;

    Ok(ThreadingResult {
        transitions,
        threaded_reads,
        threaded_pairs,
    })
}

fn thread_record(
    sequence: &[u8],
    k: usize,
    node_index: &FxHashMap<crate::dna::KmerKey, u32>,
    unitigs: &UnitigGraph,
) -> Result<Vec<u32>> {
    let kmers = canonical_kmers(sequence, k)?;
    let mut path = Vec::new();
    for pair in kmers.windows(2) {
        if pair[1].position != pair[0].position + 1 {
            continue;
        }
        let (Some(&left_id), Some(&right_id)) =
            (node_index.get(&pair[0].key), node_index.get(&pair[1].key))
        else {
            continue;
        };
        let left_state = left_id * 2 + u32::from(pair[0].reverse);
        let right_state = right_id * 2 + u32::from(pair[1].reverse);
        let Some(&unitig_id) = unitigs.edge_to_unitig.get(&(left_state, right_state)) else {
            continue;
        };
        if path.last().copied() != Some(unitig_id) {
            path.push(unitig_id);
        }
    }
    Ok(path)
}

fn reverse_unitig_path(path: &[u32], unitigs: &UnitigGraph) -> Vec<u32> {
    path.iter()
        .rev()
        .map(|&unitig_id| unitigs.reverse_unitig[unitig_id as usize])
        .collect()
}

fn add_direct_transitions(
    path: &[u32],
    transitions: &mut FxHashMap<(u32, u32), TransitionEvidence>,
) {
    for pair in path.windows(2) {
        if pair[0] == pair[1] {
            continue;
        }
        let evidence = transitions.entry((pair[0], pair[1])).or_default();
        evidence.direct_reads = evidence.direct_reads.saturating_add(1);
    }
}

fn safe_primary_paths(
    unitigs: &UnitigGraph,
    transitions: &FxHashMap<(u32, u32), TransitionEvidence>,
    min_read_support: u32,
) -> Vec<Vec<u32>> {
    let unitig_count = unitigs.unitigs.len();
    let mut outgoing: Vec<Vec<u32>> = vec![Vec::new(); unitig_count];
    let mut incoming: Vec<Vec<u32>> = vec![Vec::new(); unitig_count];

    for (&(source, target), evidence) in transitions {
        if evidence.direct_reads < min_read_support || source == target {
            continue;
        }
        outgoing[source as usize].push(target);
        incoming[target as usize].push(source);
    }
    for neighbors in outgoing.iter_mut().chain(incoming.iter_mut()) {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    let mut used = vec![false; unitig_count];
    let mut paths = Vec::new();

    for start in 0..unitig_count as u32 {
        if used[start as usize] {
            continue;
        }
        let is_internal =
            incoming[start as usize].len() == 1 && outgoing[start as usize].len() == 1;
        if is_internal {
            continue;
        }
        let path = extend_safe_path(start, &outgoing, &incoming, &mut used);
        if !path.is_empty() {
            paths.push(path);
        }
    }

    for start in 0..unitig_count as u32 {
        if used[start as usize] {
            continue;
        }
        let path = extend_safe_path(start, &outgoing, &incoming, &mut used);
        if !path.is_empty() {
            paths.push(path);
        }
    }
    paths
}

fn extend_safe_path(
    start: u32,
    outgoing: &[Vec<u32>],
    incoming: &[Vec<u32>],
    used: &mut [bool],
) -> Vec<u32> {
    let mut path = vec![start];
    used[start as usize] = true;
    let mut current = start;

    for _ in 0..used.len() {
        if outgoing[current as usize].len() != 1 {
            break;
        }
        let next = outgoing[current as usize][0];
        if incoming[next as usize].len() != 1 || used[next as usize] {
            break;
        }
        path.push(next);
        used[next as usize] = true;
        current = next;
    }
    path
}

fn deduplicate_primary_sequences(
    unitigs: &UnitigGraph,
    paths: &[Vec<u32>],
    min_length: usize,
) -> Vec<Vec<u8>> {
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut sequences = Vec::new();
    for path in paths {
        let sequence = assemble_unitig_path(unitigs, path);
        if sequence.len() < min_length {
            continue;
        }
        let reverse = reverse_complement(&sequence);
        let canonical = if reverse < sequence {
            reverse
        } else {
            sequence
        };
        if seen.insert(canonical.clone()) {
            sequences.push(canonical);
        }
    }
    sequences
}

pub fn assemble_unitig_path(unitigs: &UnitigGraph, path: &[u32]) -> Vec<u8> {
    let Some((&first, rest)) = path.split_first() else {
        return Vec::new();
    };
    let mut sequence = unitigs.unitigs[first as usize].sequence.clone();
    for &unitig_id in rest {
        let next = &unitigs.unitigs[unitig_id as usize].sequence;
        let skip = unitigs.k.min(next.len());
        sequence.extend_from_slice(&next[skip..]);
    }
    sequence
}
