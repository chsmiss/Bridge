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
    pub min_primary_support: u32,
    pub primary_dominance: f32,
    pub min_contig_length: usize,
    pub max_pairs: Option<usize>,
    pub threads: usize,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub gapped_reads: u32,
    pub read_pairs: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct TransitionCandidate {
    node: u32,
    direct_reads: u32,
    read_pairs: u32,
    solid_topology: bool,
}

#[derive(Clone, Copy, Debug)]
struct PathSelectionConfig {
    min_read_support: u32,
    min_pair_support: u32,
    min_primary_support: u32,
    min_count: u32,
    dominance: f32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ThreadedSegment {
    unitigs: Vec<u32>,
    start_edge_position: usize,
    end_edge_position: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct BubbleAllele {
    pub bubble_id: u32,
    pub allele_index: u32,
    pub unitig_id: u32,
    pub path: Vec<u32>,
    pub length: usize,
    pub mean_coverage: f32,
    pub left_support: u32,
    pub right_support: u32,
    pub physically_flanked: bool,
    #[serde(skip_serializing)]
    pub allele_sequence: Vec<u8>,
    #[serde(skip_serializing)]
    pub haplotig_sequence: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssemblyStats {
    pub version: String,
    pub kmer: KmerCountSummary,
    pub graph: GraphSummary,
    pub threaded_reads: u64,
    pub threaded_pairs: u64,
    pub direct_transitions: usize,
    pub gapped_transitions: usize,
    pub pair_bridges: usize,
    pub dominant_transitions: usize,
    pub simple_bubbles: usize,
    pub variant_alleles: usize,
    pub haplotigs: usize,
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
    pub bubble_alleles: Vec<BubbleAllele>,
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
    if !(0.5..=1.0).contains(&config.primary_dominance) {
        anyhow::bail!("primary dominance must be in 0.5..=1.0");
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
    let bubble_alleles = detect_simple_bubbles(
        &unitig_graph,
        &transitions,
        config.min_read_support,
        config.min_contig_length,
    );
    let (primary_paths, dominant_transitions) = primary_paths(
        &unitig_graph,
        &transitions,
        &bubble_alleles,
        PathSelectionConfig {
            min_read_support: config.min_read_support,
            min_pair_support: config.min_pair_support,
            min_primary_support: config.min_primary_support,
            min_count: config.min_count,
            dominance: config.primary_dominance,
        },
    );
    let mut primary_sequences =
        deduplicate_primary_sequences(&unitig_graph, &primary_paths, config.min_contig_length);
    primary_sequences
        .sort_unstable_by(|left, right| right.len().cmp(&left.len()).then(left.cmp(right)));
    timings.insert(
        "branch_aware_emission".to_string(),
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
    let gapped_transitions = transitions
        .values()
        .filter(|evidence| evidence.gapped_reads > 0)
        .count();
    let pair_bridges = transitions
        .values()
        .filter(|evidence| evidence.read_pairs >= config.min_pair_support)
        .count();
    let simple_bubbles = bubble_alleles
        .iter()
        .map(|allele| allele.bubble_id)
        .collect::<FxHashSet<_>>()
        .len();
    let haplotigs = bubble_alleles
        .iter()
        .filter(|allele| allele.haplotig_sequence.is_some())
        .count();
    timings.insert("total".to_string(), total_started.elapsed().as_secs_f64());

    let stats = AssemblyStats {
        version: env!("CARGO_PKG_VERSION").to_string(),
        kmer: kmer_set.summary,
        graph: graph_summary,
        threaded_reads,
        threaded_pairs,
        direct_transitions,
        gapped_transitions,
        pair_bridges,
        dominant_transitions,
        simple_bubbles,
        variant_alleles: bubble_alleles.len(),
        haplotigs,
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
        bubble_alleles,
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
        let left_segments = thread_record(&left.sequence, graph.k, &node_index, unitigs)?;
        if !left_segments.is_empty() {
            threaded_reads += 1;
            add_direct_transitions(&left_segments, unitigs, &mut transitions);
        }

        if let Some(right) = right {
            let right_segments = thread_record(&right.sequence, graph.k, &node_index, unitigs)?;
            if !right_segments.is_empty() {
                threaded_reads += 1;
                add_direct_transitions(&right_segments, unitigs, &mut transitions);
            }
            if let (Some(left_end), Some(right_start)) = (
                last_threaded_unitig(&left_segments),
                first_molecular_unitig(&right_segments, unitigs),
            ) {
                if left_end != right_start {
                    let evidence = transitions.entry((left_end, right_start)).or_default();
                    evidence.read_pairs = evidence.read_pairs.saturating_add(1);
                    threaded_pairs += 1;
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
) -> Result<Vec<ThreadedSegment>> {
    let kmers = canonical_kmers(sequence, k)?;
    let mut segments = Vec::new();
    let mut current = Vec::new();
    let mut segment_start = None;
    let mut segment_end = 0_usize;

    for pair in kmers.windows(2) {
        let edge_position = pair[0].position;
        let unitig_id = if pair[1].position == edge_position + 1 {
            let ids = (node_index.get(&pair[0].key), node_index.get(&pair[1].key));
            match ids {
                (Some(&left_id), Some(&right_id)) => {
                    let left_state = left_id * 2 + u32::from(pair[0].reverse);
                    let right_state = right_id * 2 + u32::from(pair[1].reverse);
                    unitigs
                        .edge_to_unitig
                        .get(&(left_state, right_state))
                        .copied()
                }
                _ => None,
            }
        } else {
            None
        };

        match unitig_id {
            Some(unitig_id) => {
                segment_start.get_or_insert(edge_position);
                segment_end = edge_position;
                if current.last().copied() != Some(unitig_id) {
                    current.push(unitig_id);
                }
            }
            None => {
                flush_threaded_segment(&mut segments, &mut current, &mut segment_start, segment_end)
            }
        }
    }
    flush_threaded_segment(&mut segments, &mut current, &mut segment_start, segment_end);
    Ok(segments)
}

fn flush_threaded_segment(
    segments: &mut Vec<ThreadedSegment>,
    current: &mut Vec<u32>,
    segment_start: &mut Option<usize>,
    segment_end: usize,
) {
    if current.is_empty() {
        *segment_start = None;
        return;
    }
    let start_edge_position = segment_start.take().unwrap_or(segment_end);
    segments.push(ThreadedSegment {
        unitigs: std::mem::take(current),
        start_edge_position,
        end_edge_position: segment_end,
    });
}

fn last_threaded_unitig(segments: &[ThreadedSegment]) -> Option<u32> {
    segments
        .last()
        .and_then(|segment| segment.unitigs.last())
        .copied()
}

fn first_molecular_unitig(segments: &[ThreadedSegment], unitigs: &UnitigGraph) -> Option<u32> {
    segments
        .last()
        .and_then(|segment| segment.unitigs.last())
        .map(|&unitig_id| unitigs.reverse_unitig[unitig_id as usize])
}

fn add_direct_transitions(
    segments: &[ThreadedSegment],
    unitigs: &UnitigGraph,
    transitions: &mut FxHashMap<(u32, u32), TransitionEvidence>,
) {
    for segment in segments {
        for pair in segment.unitigs.windows(2) {
            if pair[0] == pair[1] {
                continue;
            }
            let evidence = transitions.entry((pair[0], pair[1])).or_default();
            evidence.direct_reads = evidence.direct_reads.saturating_add(1);
        }
    }

    // Exact k-mer threading can be interrupted by a short filtered or
    // low-quality window even when the same read anchors the two adjacent
    // unitigs. Recover only a single existing UnitigGraph edge and only across
    // a very short gap; this adds physical evidence without inventing graph
    // topology or searching arbitrary alternative paths.
    const MAX_GAPPED_EDGE_KMERS: usize = 4;
    for pair in segments.windows(2) {
        let left = &pair[0];
        let right = &pair[1];
        let Some(&source) = left.unitigs.last() else {
            continue;
        };
        let Some(&target) = right.unitigs.first() else {
            continue;
        };
        if source == target {
            continue;
        }
        let missing_edges = right
            .start_edge_position
            .saturating_sub(left.end_edge_position.saturating_add(1));
        if missing_edges > MAX_GAPPED_EDGE_KMERS || !unitigs.has_edge(source, target) {
            continue;
        }
        let evidence = transitions.entry((source, target)).or_default();
        evidence.direct_reads = evidence.direct_reads.saturating_add(1);
        evidence.gapped_reads = evidence.gapped_reads.saturating_add(1);
    }
}

fn primary_paths(
    unitigs: &UnitigGraph,
    transitions: &FxHashMap<(u32, u32), TransitionEvidence>,
    bubble_alleles: &[BubbleAllele],
    selection: PathSelectionConfig,
) -> (Vec<Vec<u32>>, usize) {
    let unitig_count = unitigs.unitigs.len();
    let excluded = non_primary_bubble_alleles(bubble_alleles);
    let mut outgoing_candidates: Vec<Vec<TransitionCandidate>> = vec![Vec::new(); unitig_count];
    let mut incoming_candidates: Vec<Vec<TransitionCandidate>> = vec![Vec::new(); unitig_count];

    // Use graph adjacency as the candidate set. Every retained unitig edge was
    // observed at least min_count times, or was explicitly mercy-rescued.
    // Direct-read counts rank ambiguous exits, while excluding non-primary
    // simple-bubble alleles collapses local reconvergent variation in the
    // primary backbone without deleting it from variants/haplotigs.
    for source in 0..unitig_count as u32 {
        if excluded.contains(&source) {
            continue;
        }
        for edge_index in unitigs.out_range(source) {
            let target = unitigs.out_targets[edge_index];
            if excluded.contains(&target) || source == target {
                continue;
            }
            let evidence = transitions
                .get(&(source, target))
                .copied()
                .unwrap_or_default();
            let solid_topology = unitigs.unitigs[source as usize].min_coverage
                >= selection.min_count
                && unitigs.unitigs[target as usize].min_coverage >= selection.min_count;
            outgoing_candidates[source as usize].push(TransitionCandidate {
                node: target,
                direct_reads: evidence.direct_reads,
                read_pairs: evidence.read_pairs,
                solid_topology,
            });
            incoming_candidates[target as usize].push(TransitionCandidate {
                node: source,
                direct_reads: evidence.direct_reads,
                read_pairs: evidence.read_pairs,
                solid_topology,
            });
        }
    }

    let selected_out: Vec<Option<u32>> = outgoing_candidates
        .iter()
        .map(|candidates| {
            choose_transition(
                candidates,
                selection.min_read_support,
                selection.min_pair_support,
                selection.min_primary_support,
                selection.dominance,
            )
        })
        .collect();
    let selected_in: Vec<Option<u32>> = incoming_candidates
        .iter()
        .map(|candidates| {
            choose_transition(
                candidates,
                selection.min_read_support,
                selection.min_pair_support,
                selection.min_primary_support,
                selection.dominance,
            )
        })
        .collect();

    let mut successor = vec![None; unitig_count];
    let mut predecessor = vec![None; unitig_count];
    let mut dominant_transitions = 0_usize;
    for source in 0..unitig_count as u32 {
        if excluded.contains(&source) {
            continue;
        }
        let Some(target) = selected_out[source as usize] else {
            continue;
        };
        if selected_in[target as usize] == Some(source) {
            successor[source as usize] = Some(target);
            predecessor[target as usize] = Some(source);
            if outgoing_candidates[source as usize].len() > 1
                || incoming_candidates[target as usize].len() > 1
            {
                dominant_transitions += 1;
            }
        }
    }

    let mut used = vec![false; unitig_count];
    for &unitig_id in &excluded {
        used[unitig_id as usize] = true;
    }
    let mut paths = Vec::new();
    for start in 0..unitig_count as u32 {
        if used[start as usize] || predecessor[start as usize].is_some() {
            continue;
        }
        paths.push(extend_selected_path(start, &successor, &mut used));
    }
    for start in 0..unitig_count as u32 {
        if !used[start as usize] {
            paths.push(extend_selected_path(start, &successor, &mut used));
        }
    }
    (paths, dominant_transitions)
}

fn non_primary_bubble_alleles(bubble_alleles: &[BubbleAllele]) -> FxHashSet<u32> {
    let mut groups: FxHashMap<u32, Vec<&BubbleAllele>> = FxHashMap::default();
    for allele in bubble_alleles {
        groups.entry(allele.bubble_id).or_default().push(allele);
    }

    let mut excluded = FxHashSet::default();
    for alleles in groups.values() {
        // Sharing graph boundaries is not enough to prove a biological
        // bubble. Collapse alternatives only when at least two alleles have
        // independent direct-read support on both flanks. This prevents
        // repeat/orientation artefacts from truncating a linear primary path.
        let supported: Vec<&BubbleAllele> = alleles
            .iter()
            .copied()
            .filter(|allele| allele.physically_flanked)
            .collect();
        if supported.len() < 2 {
            continue;
        }
        let Some(primary) = supported.iter().copied().max_by(|left, right| {
            left.mean_coverage
                .total_cmp(&right.mean_coverage)
                .then_with(|| {
                    (left.left_support + left.right_support)
                        .cmp(&(right.left_support + right.right_support))
                })
                .then_with(|| right.unitig_id.cmp(&left.unitig_id))
        }) else {
            continue;
        };
        for allele in supported {
            if allele.unitig_id != primary.unitig_id {
                excluded.insert(allele.unitig_id);
            }
        }
    }
    excluded
}

fn choose_transition(
    candidates: &[TransitionCandidate],
    min_read_support: u32,
    min_pair_support: u32,
    min_primary_support: u32,
    dominance: f32,
) -> Option<u32> {
    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        let candidate = candidates[0];
        // A unique retained edge has already passed k-mer/edge filtering. It
        // may be extended without a separately materialized read transition
        // when both incident unitigs are solid. Mercy-only links still need
        // direct-read or paired-fragment support.
        return (candidate.direct_reads >= min_read_support
            || candidate.read_pairs >= min_pair_support
            || candidate.solid_topology)
            .then_some(candidate.node);
    }

    let max_direct = candidates
        .iter()
        .map(|candidate| candidate.direct_reads)
        .max()
        .unwrap_or(0);
    let use_direct = max_direct >= min_read_support;
    let mut ranked = candidates.to_vec();
    if use_direct {
        ranked.sort_unstable_by(|left, right| {
            right
                .direct_reads
                .cmp(&left.direct_reads)
                .then_with(|| right.read_pairs.cmp(&left.read_pairs))
                .then_with(|| right.solid_topology.cmp(&left.solid_topology))
                .then_with(|| left.node.cmp(&right.node))
        });
    } else {
        ranked.sort_unstable_by(|left, right| {
            right
                .read_pairs
                .cmp(&left.read_pairs)
                .then_with(|| right.direct_reads.cmp(&left.direct_reads))
                .then_with(|| right.solid_topology.cmp(&left.solid_topology))
                .then_with(|| left.node.cmp(&right.node))
        });
    }

    let total: u64 = ranked
        .iter()
        .map(|candidate| {
            u64::from(if use_direct {
                candidate.direct_reads
            } else {
                candidate.read_pairs
            })
        })
        .sum();
    let best = ranked[0];
    let best_support = if use_direct {
        best.direct_reads
    } else {
        best.read_pairs
    };
    let minimum = if use_direct {
        min_primary_support
    } else {
        min_pair_support
    };
    let fraction = best_support as f32 / total.max(1) as f32;
    (best_support >= minimum && fraction >= dominance).then_some(best.node)
}

fn extend_selected_path(start: u32, successor: &[Option<u32>], used: &mut [bool]) -> Vec<u32> {
    let mut path = Vec::new();
    let mut current = start;
    for _ in 0..used.len() {
        if used[current as usize] {
            break;
        }
        used[current as usize] = true;
        path.push(current);
        let Some(next) = successor[current as usize] else {
            break;
        };
        current = next;
    }
    path
}

fn detect_simple_bubbles(
    unitigs: &UnitigGraph,
    transitions: &FxHashMap<(u32, u32), TransitionEvidence>,
    min_read_support: u32,
    min_contig_length: usize,
) -> Vec<BubbleAllele> {
    let mut groups: FxHashMap<(u32, u32), Vec<u32>> = FxHashMap::default();
    let mut starts: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut ends: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for unitig in &unitigs.unitigs {
        groups
            .entry((unitig.start_state, unitig.end_state))
            .or_default()
            .push(unitig.id);
        starts
            .entry(unitig.start_state)
            .or_default()
            .push(unitig.id);
        ends.entry(unitig.end_state).or_default().push(unitig.id);
    }

    let mut group_entries: Vec<_> = groups.into_iter().collect();
    group_entries.sort_unstable_by_key(|(boundary, _)| *boundary);
    let mut output = Vec::new();
    let mut bubble_id = 0_u32;

    for ((start_state, end_state), mut alleles) in group_entries {
        alleles.sort_unstable();
        alleles.dedup();
        if alleles.len() < 2 || alleles.len() > 8 || start_state == end_state {
            continue;
        }
        let allele_set: FxHashSet<u32> = alleles.iter().copied().collect();
        let left_candidates: Vec<u32> = ends
            .get(&start_state)
            .into_iter()
            .flatten()
            .copied()
            .filter(|unitig| !allele_set.contains(unitig))
            .collect();
        let right_candidates: Vec<u32> = starts
            .get(&end_state)
            .into_iter()
            .flatten()
            .copied()
            .filter(|unitig| !allele_set.contains(unitig))
            .collect();
        let unique_left = (left_candidates.len() == 1).then(|| left_candidates[0]);
        let unique_right = (right_candidates.len() == 1).then(|| right_candidates[0]);

        for (allele_index, unitig_id) in alleles.into_iter().enumerate() {
            let unitig = &unitigs.unitigs[unitig_id as usize];
            let left_support = unique_left
                .and_then(|left| transitions.get(&(left, unitig_id)))
                .map_or(0, |evidence| evidence.direct_reads);
            let right_support = unique_right
                .and_then(|right| transitions.get(&(unitig_id, right)))
                .map_or(0, |evidence| evidence.direct_reads);
            let physically_flanked = unique_left.is_some()
                && unique_right.is_some()
                && left_support >= min_read_support
                && right_support >= min_read_support;
            let mut path = Vec::new();
            if let Some(left) = unique_left.filter(|_| physically_flanked) {
                path.push(left);
            }
            path.push(unitig_id);
            if let Some(right) = unique_right.filter(|_| physically_flanked) {
                path.push(right);
            }
            let haplotig_sequence = if physically_flanked {
                let sequence = assemble_unitig_path(unitigs, &path);
                (sequence.len() >= min_contig_length).then_some(sequence)
            } else {
                None
            };
            output.push(BubbleAllele {
                bubble_id,
                allele_index: allele_index as u32,
                unitig_id,
                path,
                length: unitig.length,
                mean_coverage: unitig.mean_coverage,
                left_support,
                right_support,
                physically_flanked,
                allele_sequence: unitig.sequence.clone(),
                haplotig_sequence,
            });
        }
        bubble_id += 1;
    }
    output
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

#[cfg(test)]
mod transition_selection_tests {
    use super::*;

    fn candidate(
        node: u32,
        direct_reads: u32,
        read_pairs: u32,
        solid_topology: bool,
    ) -> TransitionCandidate {
        TransitionCandidate {
            node,
            direct_reads,
            read_pairs,
            solid_topology,
        }
    }

    #[test]
    fn unique_solid_topology_is_safe_without_materialized_transition() {
        let candidates = [candidate(7, 0, 0, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(7));
    }

    #[test]
    fn unique_mercy_like_edge_requires_physical_support() {
        let candidates = [candidate(7, 0, 0, false)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), None);
    }

    #[test]
    fn pair_support_can_resolve_an_adjacent_edge_when_reads_do_not() {
        let candidates = [candidate(3, 0, 7, false), candidate(4, 0, 1, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(3));
    }

    #[test]
    fn direct_reads_take_precedence_when_they_clear_the_read_gate() {
        let candidates = [candidate(3, 6, 0, true), candidate(4, 1, 20, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), Some(3));
    }

    #[test]
    fn ambiguous_physical_evidence_remains_unresolved() {
        let candidates = [candidate(3, 0, 5, true), candidate(4, 0, 4, true)];
        assert_eq!(choose_transition(&candidates, 2, 2, 5, 0.75), None);
    }
}

#[cfg(test)]
mod gapped_threading_tests {
    use super::*;
    use crate::graph::{Unitig, UnitigGraph};

    fn unitig(id: u32) -> Unitig {
        Unitig {
            id,
            states: Vec::new(),
            sequence: b"ACGT".to_vec(),
            start_state: id,
            end_state: id + 1,
            length: 4,
            mean_coverage: 3.0,
            min_coverage: 3,
            max_coverage: 3,
            circular: false,
        }
    }

    fn two_unitig_graph(linked: bool) -> UnitigGraph {
        UnitigGraph {
            k: 3,
            unitigs: vec![unitig(0), unitig(1)],
            edge_to_unitig: FxHashMap::default(),
            reverse_unitig: vec![0, 1],
            out_offsets: if linked { vec![0, 1, 1] } else { vec![0, 0, 0] },
            out_targets: if linked { vec![1] } else { Vec::new() },
            indegree: if linked { vec![0, 1] } else { vec![0, 0] },
        }
    }

    #[test]
    fn short_gap_on_existing_edge_adds_physical_transition() {
        let segments = vec![
            ThreadedSegment {
                unitigs: vec![0],
                start_edge_position: 0,
                end_edge_position: 2,
            },
            ThreadedSegment {
                unitigs: vec![1],
                start_edge_position: 5,
                end_edge_position: 6,
            },
        ];
        let mut transitions = FxHashMap::default();
        add_direct_transitions(&segments, &two_unitig_graph(true), &mut transitions);
        let evidence = transitions.get(&(0, 1)).copied().unwrap_or_default();
        assert_eq!(evidence.direct_reads, 1);
        assert_eq!(evidence.gapped_reads, 1);
    }

    #[test]
    fn gap_does_not_create_missing_graph_edge() {
        let segments = vec![
            ThreadedSegment {
                unitigs: vec![0],
                start_edge_position: 0,
                end_edge_position: 2,
            },
            ThreadedSegment {
                unitigs: vec![1],
                start_edge_position: 5,
                end_edge_position: 6,
            },
        ];
        let mut transitions = FxHashMap::default();
        add_direct_transitions(&segments, &two_unitig_graph(false), &mut transitions);
        assert!(!transitions.contains_key(&(0, 1)));
    }
}
