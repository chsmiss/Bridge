from pathlib import Path
import re

assembler = Path("src/assembler.rs")
text = assembler.read_text()

text = text.replace(
'''pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub read_pairs: u32,
}
''',
'''pub struct TransitionEvidence {
    pub direct_reads: u32,
    pub gapped_reads: u32,
    pub read_pairs: u32,
}
''',
1,
)

marker = '''struct PathSelectionConfig {
    min_read_support: u32,
    min_pair_support: u32,
    min_primary_support: u32,
    min_count: u32,
    dominance: f32,
}
'''
insert = marker + '''
#[derive(Clone, Debug, Eq, PartialEq)]
struct ThreadedSegment {
    unitigs: Vec<u32>,
    start_edge_position: usize,
    end_edge_position: usize,
}
'''
if marker not in text:
    raise SystemExit("PathSelectionConfig marker missing")
text = text.replace(marker, insert, 1)

text = text.replace(
'''    pub direct_transitions: usize,
    pub pair_bridges: usize,
''',
'''    pub direct_transitions: usize,
    pub gapped_transitions: usize,
    pub pair_bridges: usize,
''',
1,
)

text = text.replace(
'''        if !left_segments.is_empty() {
            threaded_reads += 1;
            add_direct_transitions(&left_segments, &mut transitions);
        }
''',
'''        if !left_segments.is_empty() {
            threaded_reads += 1;
            add_direct_transitions(&left_segments, unitigs, &mut transitions);
        }
''',
1,
)
text = text.replace(
'''            if !right_segments.is_empty() {
                threaded_reads += 1;
                add_direct_transitions(&right_segments, &mut transitions);
            }
''',
'''            if !right_segments.is_empty() {
                threaded_reads += 1;
                add_direct_transitions(&right_segments, unitigs, &mut transitions);
            }
''',
1,
)

pattern = re.compile(
    r"fn thread_record\(\n.*?\n\}\n\nfn last_threaded_unitig\(.*?\n\}\n\nfn first_molecular_unitig\(.*?\n\}\n\nfn add_direct_transitions\(.*?\n\}\n",
    re.S,
)
replacement = '''fn thread_record(
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
            None => flush_threaded_segment(
                &mut segments,
                &mut current,
                &mut segment_start,
                segment_end,
            ),
        }
    }
    flush_threaded_segment(
        &mut segments,
        &mut current,
        &mut segment_start,
        segment_end,
    );
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
'''
text, count = pattern.subn(replacement, text, count=1)
if count != 1:
    raise SystemExit(f"threading block replacements: {count}")

text = text.replace(
'''    let pair_bridges = transitions
        .values()
        .filter(|evidence| evidence.read_pairs >= config.min_pair_support)
        .count();
''',
'''    let gapped_transitions = transitions
        .values()
        .filter(|evidence| evidence.gapped_reads > 0)
        .count();
    let pair_bridges = transitions
        .values()
        .filter(|evidence| evidence.read_pairs >= config.min_pair_support)
        .count();
''',
1,
)
text = text.replace(
'''        direct_transitions,
        pair_bridges,
''',
'''        direct_transitions,
        gapped_transitions,
        pair_bridges,
''',
1,
)

text += '''

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
'''

assembler.write_text(text)

graph = Path("src/graph.rs")
text = graph.read_text()
marker = '''    #[inline]
    pub fn outdegree(&self, unitig: u32) -> usize {
        self.out_range(unitig).len()
    }
'''
insert = marker + '''
    #[inline]
    pub fn has_edge(&self, source: u32, target: u32) -> bool {
        self.out_targets[self.out_range(source)].binary_search(&target).is_ok()
    }
'''
if marker not in text:
    raise SystemExit("UnitigGraph outdegree marker missing")
graph.write_text(text.replace(marker, insert, 1))

output = Path("src/output.rs")
text = output.read_text()
old = '''                "L\\tu{}\\t+\\tu{}\\t+\\t{}M\\tDR:i:{}\\tPE:i:{}",
                source, target, product.unitig_graph.k, evidence.direct_reads, evidence.read_pairs
'''
new = '''                "L\\tu{}\\t+\\tu{}\\t+\\t{}M\\tDR:i:{}\\tGR:i:{}\\tPE:i:{}",
                source,
                target,
                product.unitig_graph.k,
                evidence.direct_reads,
                evidence.gapped_reads,
                evidence.read_pairs
'''
if old not in text:
    raise SystemExit("GFA link output block missing")
output.write_text(text.replace(old, new, 1))
