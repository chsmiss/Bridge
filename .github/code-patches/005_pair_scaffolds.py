from pathlib import Path

scaffold = Path("src/scaffold.rs")
scaffold.write_text(r'''use crate::assembler::TransitionEvidence;
use crate::dna::reverse_complement;
use crate::graph::UnitigGraph;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct ScaffoldLink {
    pub source_component: u32,
    pub target_component: u32,
    pub pair_support: u32,
    pub gap_bases: usize,
}

#[derive(Debug)]
pub struct ScaffoldResult {
    pub sequences: Vec<Vec<u8>>,
    pub links: Vec<ScaffoldLink>,
}

#[derive(Clone, Copy, Debug)]
struct Candidate {
    component: u32,
    support: u32,
}

pub fn build_scaffolds(
    unitigs: &UnitigGraph,
    primary_paths: &[Vec<u32>],
    transitions: &FxHashMap<(u32, u32), TransitionEvidence>,
    min_pair_support: u32,
    dominance: f32,
    gap_bases: usize,
) -> ScaffoldResult {
    if primary_paths.is_empty() {
        return ScaffoldResult {
            sequences: Vec::new(),
            links: Vec::new(),
        };
    }

    let path_sequences: Vec<Vec<u8>> = primary_paths
        .iter()
        .map(|path| assemble_path(unitigs, path))
        .collect();
    let canonical_sequences: Vec<Vec<u8>> = path_sequences
        .iter()
        .map(|sequence| canonical_sequence(sequence))
        .collect();

    let mut start_to_component = FxHashMap::default();
    let mut end_to_component = FxHashMap::default();
    for (component, path) in primary_paths.iter().enumerate() {
        let Some((&start, rest)) = path.split_first() else {
            continue;
        };
        let end = rest.last().copied().unwrap_or(start);
        start_to_component.insert(start, component as u32);
        end_to_component.insert(end, component as u32);
    }

    let mut pair_support: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    for (&(source, target), evidence) in transitions {
        if evidence.read_pairs < min_pair_support {
            continue;
        }
        let (Some(&source_component), Some(&target_component)) = (
            end_to_component.get(&source),
            start_to_component.get(&target),
        ) else {
            continue;
        };
        if source_component == target_component
            || canonical_sequences[source_component as usize]
                == canonical_sequences[target_component as usize]
        {
            continue;
        }
        let entry = pair_support
            .entry((source_component, target_component))
            .or_insert(0);
        *entry = entry.saturating_add(evidence.read_pairs);
    }

    let component_count = primary_paths.len();
    let mut outgoing = vec![Vec::new(); component_count];
    let mut incoming = vec![Vec::new(); component_count];
    for (&(source, target), &support) in &pair_support {
        outgoing[source as usize].push(Candidate {
            component: target,
            support,
        });
        incoming[target as usize].push(Candidate {
            component: source,
            support,
        });
    }

    let selected_out: Vec<Option<u32>> = outgoing
        .iter()
        .map(|candidates| choose_candidate(candidates, min_pair_support, dominance))
        .collect();
    let selected_in: Vec<Option<u32>> = incoming
        .iter()
        .map(|candidates| choose_candidate(candidates, min_pair_support, dominance))
        .collect();

    let mut successor = vec![None; component_count];
    let mut predecessor = vec![None; component_count];
    for source in 0..component_count as u32 {
        let Some(target) = selected_out[source as usize] else {
            continue;
        };
        if selected_in[target as usize] == Some(source) {
            successor[source as usize] = Some(target);
            predecessor[target as usize] = Some(source);
        }
    }

    let mut used = vec![false; component_count];
    let mut scaffold_components = Vec::new();
    for start in 0..component_count as u32 {
        if used[start as usize] || predecessor[start as usize].is_some() {
            continue;
        }
        scaffold_components.push(extend_component_path(start, &successor, &mut used));
    }
    // Reciprocal links can still form a cycle. Do not turn a cycle into an
    // arbitrary linear scaffold: leave every remaining component independent.
    for component in 0..component_count as u32 {
        if !used[component as usize] {
            used[component as usize] = true;
            scaffold_components.push(vec![component]);
        }
    }

    let mut links = Vec::new();
    let mut seen_sequences: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut sequences = Vec::new();
    for components in scaffold_components {
        let Some((&first, rest)) = components.split_first() else {
            continue;
        };
        let mut sequence = path_sequences[first as usize].clone();
        for &target in rest {
            let source = components
                .iter()
                .position(|&component| component == target)
                .and_then(|index| index.checked_sub(1))
                .map(|index| components[index])
                .unwrap_or(first);
            let support = pair_support.get(&(source, target)).copied().unwrap_or(0);
            sequence.extend(std::iter::repeat_n(b'N', gap_bases));
            sequence.extend_from_slice(&path_sequences[target as usize]);
            links.push(ScaffoldLink {
                source_component: source,
                target_component: target,
                pair_support: support,
                gap_bases,
            });
        }
        let canonical = canonical_sequence(&sequence);
        if seen_sequences.insert(canonical.clone()) {
            sequences.push(canonical);
        }
    }

    ScaffoldResult { sequences, links }
}

fn choose_candidate(candidates: &[Candidate], min_support: u32, dominance: f32) -> Option<u32> {
    if candidates.is_empty() {
        return None;
    }
    let total: u64 = candidates
        .iter()
        .map(|candidate| u64::from(candidate.support))
        .sum();
    let best = candidates.iter().copied().max_by(|left, right| {
        left.support
            .cmp(&right.support)
            .then_with(|| right.component.cmp(&left.component))
    })?;
    let fraction = best.support as f32 / total.max(1) as f32;
    (best.support >= min_support && fraction >= dominance).then_some(best.component)
}

fn extend_component_path(start: u32, successor: &[Option<u32>], used: &mut [bool]) -> Vec<u32> {
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

fn assemble_path(unitigs: &UnitigGraph, path: &[u32]) -> Vec<u8> {
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

fn canonical_sequence(sequence: &[u8]) -> Vec<u8> {
    let reverse = reverse_complement(sequence);
    if reverse.as_slice() < sequence {
        reverse
    } else {
        sequence.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Unitig;

    fn unitig(id: u32, sequence: &[u8]) -> Unitig {
        Unitig {
            id,
            states: Vec::new(),
            sequence: sequence.to_vec(),
            start_state: id,
            end_state: id + 1,
            length: sequence.len(),
            mean_coverage: 4.0,
            min_coverage: 4,
            max_coverage: 4,
            circular: false,
        }
    }

    fn graph() -> UnitigGraph {
        UnitigGraph {
            k: 3,
            unitigs: vec![
                unitig(0, b"AAAC"),
                unitig(1, b"CCCG"),
                unitig(2, b"GGGT"),
            ],
            edge_to_unitig: FxHashMap::default(),
            reverse_unitig: vec![0, 1, 2],
            out_offsets: vec![0, 0, 0, 0],
            out_targets: Vec::new(),
            indegree: vec![0, 0, 0],
        }
    }

    #[test]
    fn reciprocal_pair_link_builds_an_n_gap_scaffold() {
        let mut transitions = FxHashMap::default();
        transitions.insert(
            (0, 1),
            TransitionEvidence {
                read_pairs: 8,
                ..TransitionEvidence::default()
            },
        );
        let result = build_scaffolds(
            &graph(),
            &[vec![0], vec![1]],
            &transitions,
            3,
            0.75,
            10,
        );
        assert_eq!(result.links.len(), 1);
        assert_eq!(result.sequences.len(), 1);
        assert_eq!(result.sequences[0].len(), 18);
        assert!(result.sequences[0].windows(10).any(|window| window == b"NNNNNNNNNN"));
    }

    #[test]
    fn ambiguous_pair_exits_remain_separate() {
        let mut transitions = FxHashMap::default();
        for (target, support) in [(1, 5), (2, 4)] {
            transitions.insert(
                (0, target),
                TransitionEvidence {
                    read_pairs: support,
                    ..TransitionEvidence::default()
                },
            );
        }
        let result = build_scaffolds(
            &graph(),
            &[vec![0], vec![1], vec![2]],
            &transitions,
            3,
            0.75,
            10,
        );
        assert!(result.links.is_empty());
        assert_eq!(result.sequences.len(), 3);
    }
}
''')

lib = Path("src/lib.rs")
text = lib.read_text()
if "pub mod scaffold;\n" not in text:
    text += "pub mod scaffold;\n"
lib.write_text(text)

main = Path("src/main.rs")
text = main.read_text()
text = text.replace(
'''        #[arg(long, default_value_t = 200)]
        min_contig_length: usize,
''',
'''        #[arg(long, default_value_t = 200)]
        min_contig_length: usize,
        #[arg(long, default_value_t = 100)]
        scaffold_gap_bases: usize,
''',
1,
)
text = text.replace(
'''            min_contig_length,
            max_pairs,
''',
'''            min_contig_length,
            scaffold_gap_bases,
            max_pairs,
''',
1,
)
text = text.replace(
'''                min_contig_length,
                max_pairs,
''',
'''                min_contig_length,
                scaffold_gap_bases,
                max_pairs,
''',
1,
)
main.write_text(text)

assembler = Path("src/assembler.rs")
text = assembler.read_text()
text = text.replace(
'''use crate::kmer::{count_and_filter, KmerCountSummary, KmerFilterConfig};
''',
'''use crate::kmer::{count_and_filter, KmerCountSummary, KmerFilterConfig};
use crate::scaffold::{build_scaffolds, ScaffoldLink};
''',
1,
)
text = text.replace(
'''    pub min_contig_length: usize,
    pub max_pairs: Option<usize>,
''',
'''    pub min_contig_length: usize,
    pub scaffold_gap_bases: usize,
    pub max_pairs: Option<usize>,
''',
1,
)
text = text.replace(
'''    pub largest_primary: usize,
    pub timings_seconds: FxHashMap<String, f64>,
''',
'''    pub largest_primary: usize,
    pub scaffolds: usize,
    pub scaffold_links: usize,
    pub scaffold_bases: usize,
    pub scaffold_n50: usize,
    pub largest_scaffold: usize,
    pub scaffold_gap_bases: usize,
    pub timings_seconds: FxHashMap<String, f64>,
''',
1,
)
text = text.replace(
'''    pub primary_sequences: Vec<Vec<u8>>,
    pub bubble_alleles: Vec<BubbleAllele>,
''',
'''    pub primary_sequences: Vec<Vec<u8>>,
    pub scaffold_sequences: Vec<Vec<u8>>,
    pub scaffold_links: Vec<ScaffoldLink>,
    pub bubble_alleles: Vec<BubbleAllele>,
''',
1,
)
marker = '''    primary_sequences
        .sort_unstable_by(|left, right| right.len().cmp(&left.len()).then(left.cmp(right)));
    timings.insert(
        "branch_aware_emission".to_string(),
        started.elapsed().as_secs_f64(),
    );

'''
insert = marker + '''    let started = Instant::now();
    let scaffold_result = build_scaffolds(
        &unitig_graph,
        &primary_paths,
        &transitions,
        config.min_pair_support,
        config.primary_dominance,
        config.scaffold_gap_bases,
    );
    let mut scaffold_sequences = scaffold_result.sequences;
    scaffold_sequences
        .sort_unstable_by(|left, right| right.len().cmp(&left.len()).then(left.cmp(right)));
    let scaffold_links = scaffold_result.links;
    timings.insert(
        "paired_scaffolding".to_string(),
        started.elapsed().as_secs_f64(),
    );

'''
if marker not in text:
    raise SystemExit("branch-aware emission marker not found")
text = text.replace(marker, insert, 1)
marker = '''    let primary_n50 = crate::graph::n50(&mut primary_lengths);
'''
insert = marker + '''    let mut scaffold_lengths: Vec<usize> = scaffold_sequences.iter().map(Vec::len).collect();
    let scaffold_bases = scaffold_lengths.iter().sum();
    let largest_scaffold = scaffold_lengths.iter().copied().max().unwrap_or(0);
    let scaffold_n50 = crate::graph::n50(&mut scaffold_lengths);
'''
text = text.replace(marker, insert, 1)
text = text.replace(
'''        largest_primary,
        timings_seconds: timings,
''',
'''        largest_primary,
        scaffolds: scaffold_sequences.len(),
        scaffold_links: scaffold_links.len(),
        scaffold_bases,
        scaffold_n50,
        largest_scaffold,
        scaffold_gap_bases: scaffold_links
            .iter()
            .map(|link| link.gap_bases)
            .sum(),
        timings_seconds: timings,
''',
1,
)
text = text.replace(
'''        primary_sequences,
        bubble_alleles,
''',
'''        primary_sequences,
        scaffold_sequences,
        scaffold_links,
        bubble_alleles,
''',
1,
)
assembler.write_text(text)

output = Path("src/output.rs")
text = output.read_text()
text = text.replace(
'''use crate::assembler::{AssemblyProduct, BubbleAllele};
''',
'''use crate::assembler::{AssemblyProduct, BubbleAllele};
use crate::scaffold::ScaffoldLink;
''',
1,
)
marker = '''    write_fasta(
        &output_dir.join("unitigs.fasta"),
'''
insert = '''    write_fasta(
        &output_dir.join("primary_scaffolds.fasta"),
        product
            .scaffold_sequences
            .iter()
            .enumerate()
            .map(|(index, sequence)| {
                (
                    format!("scaffold_{:06} len={}", index + 1, sequence.len()),
                    sequence.as_slice(),
                )
            }),
    )?;
    write_scaffold_links(
        &output_dir.join("scaffold_links.tsv"),
        &product.scaffold_links,
    )?;

''' + marker
if marker not in text:
    raise SystemExit("unitig output marker missing")
text = text.replace(marker, insert, 1)
marker = '''fn write_bubble_table(path: &Path, alleles: &[BubbleAllele]) -> Result<()> {
'''
insert = '''fn write_scaffold_links(path: &Path, links: &[ScaffoldLink]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "source_component\\ttarget_component\\tpair_support\\tgap_bases"
    )?;
    for link in links {
        writeln!(
            writer,
            "{}\\t{}\\t{}\\t{}",
            link.source_component,
            link.target_component,
            link.pair_support,
            link.gap_bases
        )?;
    }
    writer.flush()?;
    Ok(())
}

''' + marker
text = text.replace(marker, insert, 1)
output.write_text(text)

tests = Path("tests/integration.rs")
text = tests.read_text()
text = text.replace(
'''        min_contig_length: 20,
        max_pairs: None,
''',
'''        min_contig_length: 20,
        scaffold_gap_bases: 100,
        max_pairs: None,
''',
1,
)
tests.write_text(text)
