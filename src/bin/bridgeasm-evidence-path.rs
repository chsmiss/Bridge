use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Conservative,
    Balanced,
    Exploratory,
}

#[derive(Parser, Debug)]
#[command(
    name = "bridgeasm-evidence-path",
    about = "Fuse nucleotide, protein-assembly, ESM, and DNA-LM evidence on existing GFA links"
)]
struct Cli {
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    edge_evidence: PathBuf,
    #[arg(long)]
    esm_scores: Option<PathBuf>,
    #[arg(long)]
    dna_lm_scores: Option<PathBuf>,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    report: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = Mode::Balanced)]
    mode: Mode,
    #[arg(long, default_value_t = 1)]
    min_direct: u32,
    #[arg(long, default_value_t = 1)]
    min_pair: u32,
    #[arg(long, default_value_t = 1)]
    min_gapped: u32,
    #[arg(long, default_value_t = 0.18)]
    min_coverage_ratio: f32,
    #[arg(long, default_value_t = 0.45)]
    min_protein_score: f32,
    #[arg(long, default_value_t = 0.45)]
    max_protein_ambiguity: f32,
    #[arg(long, default_value_t = 0.90)]
    min_frame_consistency: f32,
    #[arg(long, default_value_t = 6)]
    min_protein_kmers: u32,
    #[arg(long, default_value_t = 2400.0)]
    protein_weight: f64,
    #[arg(long, default_value_t = 300.0)]
    esm_weight: f64,
    #[arg(long, default_value_t = 180.0)]
    dna_lm_weight: f64,
    #[arg(long, default_value_t = 200)]
    min_length: usize,
}

#[derive(Clone, Debug)]
struct Segment {
    name: String,
    sequence: Vec<u8>,
    coverage: f32,
}

#[derive(Clone, Copy, Debug)]
struct Link {
    source: usize,
    target: usize,
    overlap: usize,
    direct: u32,
    gapped: u32,
    pairs: u32,
}

#[derive(Clone, Debug, Default)]
struct EdgeEvidence {
    protein_score: f32,
    protein_ambiguity: f32,
    frame_consistency: f32,
    unique_kmers: u32,
    breakpoint_class: String,
    protein_id: String,
    esm_delta: f32,
    dna_lm_delta: f32,
}

#[derive(Clone, Debug)]
struct RankedLink {
    link: Link,
    score: f64,
    coverage_ratio: f32,
    physical: bool,
    protein_ok: bool,
    evidence: EdgeEvidence,
}

#[derive(Clone, Debug)]
struct RawLink {
    source: String,
    target: String,
    source_orientation: String,
    target_orientation: String,
    overlap: String,
    tags: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate(&cli)?;
    let (segments, links, name_to_id, skipped_oriented) = read_gfa(&cli.gfa)?;
    if segments.is_empty() {
        bail!("GFA contains no segments");
    }
    let mut evidence = read_edge_evidence(&cli.edge_evidence, &name_to_id)?;
    if let Some(path) = cli.esm_scores.as_ref() {
        merge_scalar_scores(path, "esm_delta", &name_to_id, &mut evidence, |row, value| {
            row.esm_delta = value;
        })?;
    }
    if let Some(path) = cli.dna_lm_scores.as_ref() {
        merge_scalar_scores(
            path,
            "dna_lm_delta",
            &name_to_id,
            &mut evidence,
            |row, value| row.dna_lm_delta = value,
        )?;
    }

    let (selected, ranked) = select_links(&segments, &links, &evidence, &cli);
    let paths = build_paths(segments.len(), &selected);
    let mut sequences = Vec::new();
    let mut seen = FxHashSet::default();
    for path in &paths {
        let sequence = assemble_path(path, &segments, &selected);
        if sequence.len() < cli.min_length {
            continue;
        }
        let canonical = canonical_sequence(sequence);
        if seen.insert(canonical.clone()) {
            sequences.push((path.clone(), canonical));
        }
    }
    sequences.sort_unstable_by(|left, right| {
        right
            .1
            .len()
            .cmp(&left.1.len())
            .then_with(|| left.1.cmp(&right.1))
    });
    write_fasta(&cli.output, &sequences)?;
    if let Some(report) = cli.report.as_ref() {
        write_report(report, &segments, &ranked, &selected)?;
    }

    let mut lengths: Vec<usize> = sequences.iter().map(|(_, sequence)| sequence.len()).collect();
    let total: usize = lengths.iter().sum();
    let n50 = n50(&mut lengths);
    let largest = lengths.iter().copied().max().unwrap_or(0);
    let selected_protein = selected
        .values()
        .filter(|link| {
            evidence
                .get(&(link.source, link.target))
                .is_some_and(|row| protein_usable(row, &cli))
        })
        .count();
    eprintln!(
        "evidence-path mode={:?}: contigs={} bp={} N50={} largest={} selected_edges={} protein_supported_edges={} skipped_oriented_links={}",
        cli.mode,
        sequences.len(),
        total,
        n50,
        largest,
        selected.len(),
        selected_protein,
        skipped_oriented
    );
    Ok(())
}

fn validate(cli: &Cli) -> Result<()> {
    for (name, value) in [
        ("min-coverage-ratio", cli.min_coverage_ratio),
        ("min-protein-score", cli.min_protein_score),
        ("max-protein-ambiguity", cli.max_protein_ambiguity),
        ("min-frame-consistency", cli.min_frame_consistency),
    ] {
        if !(0.0..=1.0).contains(&value) {
            bail!("{name} must be in 0..=1");
        }
    }
    Ok(())
}

fn read_gfa(
    path: &Path,
) -> Result<(Vec<Segment>, Vec<Link>, FxHashMap<String, usize>, usize)> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut segments = Vec::new();
    let mut name_to_id = FxHashMap::default();
    let mut raw_links = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with("H\t") || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                if fields[2] == "*" {
                    bail!("GFA segment {} has no sequence", fields[1]);
                }
                let name = fields[1].to_string();
                let coverage = tag_f32(&fields[3..], "KC:f:").unwrap_or(1.0);
                let id = segments.len();
                name_to_id.insert(name.clone(), id);
                segments.push(Segment {
                    name,
                    sequence: fields[2].as_bytes().to_vec(),
                    coverage,
                });
            }
            Some("L") if fields.len() >= 6 => raw_links.push(RawLink {
                source: fields[1].to_string(),
                source_orientation: fields[2].to_string(),
                target: fields[3].to_string(),
                target_orientation: fields[4].to_string(),
                overlap: fields[5].to_string(),
                tags: fields[6..].iter().map(|value| (*value).to_string()).collect(),
            }),
            _ => {}
        }
    }

    let mut links = Vec::with_capacity(raw_links.len());
    let mut skipped_oriented = 0_usize;
    for raw in raw_links {
        if raw.source_orientation != "+" || raw.target_orientation != "+" {
            skipped_oriented += 1;
            continue;
        }
        let Some(&source) = name_to_id.get(&raw.source) else {
            continue;
        };
        let Some(&target) = name_to_id.get(&raw.target) else {
            continue;
        };
        if source == target {
            continue;
        }
        let tag_refs: Vec<&str> = raw.tags.iter().map(String::as_str).collect();
        links.push(Link {
            source,
            target,
            overlap: parse_overlap(&raw.overlap),
            direct: tag_u32(&tag_refs, "DR:i:").unwrap_or(0),
            gapped: tag_u32(&tag_refs, "GR:i:").unwrap_or(0),
            pairs: tag_u32(&tag_refs, "PE:i:").unwrap_or(0),
        });
    }
    Ok((segments, links, name_to_id, skipped_oriented))
}

fn parse_overlap(value: &str) -> usize {
    value
        .strip_suffix('M')
        .and_then(|number| number.parse::<usize>().ok())
        .unwrap_or(0)
}

fn tag_u32(fields: &[&str], prefix: &str) -> Option<u32> {
    fields
        .iter()
        .find_map(|field| field.strip_prefix(prefix)?.parse().ok())
}

fn tag_f32(fields: &[&str], prefix: &str) -> Option<f32> {
    fields
        .iter()
        .find_map(|field| field.strip_prefix(prefix)?.parse().ok())
}

fn header_map(header: &str) -> FxHashMap<String, usize> {
    header
        .split('\t')
        .enumerate()
        .map(|(index, name)| (name.trim().to_string(), index))
        .collect()
}

fn required_column(columns: &FxHashMap<String, usize>, name: &str) -> Result<usize> {
    columns
        .get(name)
        .copied()
        .with_context(|| format!("TSV is missing required column {name:?}"))
}

fn optional_field<'a>(fields: &'a [&str], columns: &FxHashMap<String, usize>, name: &str) -> &'a str {
    columns
        .get(name)
        .and_then(|index| fields.get(*index))
        .copied()
        .unwrap_or("")
}

fn parse_f32(value: &str, default: f32) -> f32 {
    value.parse::<f32>().unwrap_or(default)
}

fn parse_u32(value: &str, default: u32) -> u32 {
    value.parse::<u32>().unwrap_or(default)
}

fn read_edge_evidence(
    path: &Path,
    name_to_id: &FxHashMap<String, usize>,
) -> Result<FxHashMap<(usize, usize), EdgeEvidence>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .with_context(|| format!("{} is empty", path.display()))?;
    let columns = header_map(&header);
    let source_column = required_column(&columns, "source")?;
    let target_column = required_column(&columns, "target")?;
    let mut evidence = FxHashMap::default();

    for line in lines {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let Some(source_name) = fields.get(source_column) else {
            continue;
        };
        let Some(target_name) = fields.get(target_column) else {
            continue;
        };
        let Some(&source) = name_to_id.get(*source_name) else {
            continue;
        };
        let Some(&target) = name_to_id.get(*target_name) else {
            continue;
        };
        evidence.insert(
            (source, target),
            EdgeEvidence {
                protein_score: parse_f32(
                    optional_field(&fields, &columns, "protein_score"),
                    0.0,
                ),
                protein_ambiguity: parse_f32(
                    optional_field(&fields, &columns, "ambiguity"),
                    1.0,
                ),
                frame_consistency: parse_f32(
                    optional_field(&fields, &columns, "frame_consistency"),
                    0.0,
                ),
                unique_kmers: parse_u32(
                    optional_field(&fields, &columns, "unique_kmers"),
                    0,
                ),
                breakpoint_class: optional_field(&fields, &columns, "breakpoint_class")
                    .to_string(),
                protein_id: optional_field(&fields, &columns, "protein_id").to_string(),
                esm_delta: parse_f32(optional_field(&fields, &columns, "esm_delta"), 0.0),
                dna_lm_delta: parse_f32(
                    optional_field(&fields, &columns, "dna_lm_delta"),
                    0.0,
                ),
            },
        );
    }
    Ok(evidence)
}

fn merge_scalar_scores<F>(
    path: &Path,
    value_column_name: &str,
    name_to_id: &FxHashMap<String, usize>,
    evidence: &mut FxHashMap<(usize, usize), EdgeEvidence>,
    mut assign: F,
) -> Result<()>
where
    F: FnMut(&mut EdgeEvidence, f32),
{
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .with_context(|| format!("{} is empty", path.display()))?;
    let columns = header_map(&header);
    let source_column = required_column(&columns, "source")?;
    let target_column = required_column(&columns, "target")?;
    let value_column = required_column(&columns, value_column_name)?;
    for line in lines {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(source_name), Some(target_name), Some(value)) = (
            fields.get(source_column),
            fields.get(target_column),
            fields.get(value_column),
        ) else {
            continue;
        };
        let (Some(&source), Some(&target)) =
            (name_to_id.get(*source_name), name_to_id.get(*target_name))
        else {
            continue;
        };
        let parsed = value.parse::<f32>().unwrap_or(0.0);
        assign(evidence.entry((source, target)).or_default(), parsed);
    }
    Ok(())
}

fn protein_usable(row: &EdgeEvidence, cli: &Cli) -> bool {
    row.protein_score >= cli.min_protein_score
        && row.protein_ambiguity <= cli.max_protein_ambiguity
        && row.frame_consistency >= cli.min_frame_consistency
        && row.unique_kmers >= cli.min_protein_kmers
        && row.breakpoint_class == "same_orf_supported"
}

fn select_links(
    segments: &[Segment],
    links: &[Link],
    evidence: &FxHashMap<(usize, usize), EdgeEvidence>,
    cli: &Cli,
) -> (FxHashMap<usize, Link>, Vec<RankedLink>) {
    let node_count = segments.len();
    let mut outgoing: Vec<Vec<Link>> = vec![Vec::new(); node_count];
    let mut incoming: Vec<Vec<Link>> = vec![Vec::new(); node_count];
    for &link in links {
        outgoing[link.source].push(link);
        incoming[link.target].push(link);
    }

    let mut ranked = Vec::new();
    for &link in links {
        let row = evidence
            .get(&(link.source, link.target))
            .cloned()
            .unwrap_or_default();
        let source_cov = segments[link.source].coverage.max(0.001);
        let target_cov = segments[link.target].coverage.max(0.001);
        let coverage_ratio = source_cov.min(target_cov) / source_cov.max(target_cov);
        let physical = link.direct >= cli.min_direct
            || link.gapped >= cli.min_gapped
            || link.pairs >= cli.min_pair;
        let protein_ok = protein_usable(&row, cli);
        let topology_unique = outgoing[link.source].len() == 1 && incoming[link.target].len() == 1;

        let eligible = match cli.mode {
            Mode::Conservative => physical && coverage_ratio >= cli.min_coverage_ratio,
            Mode::Balanced => {
                coverage_ratio >= cli.min_coverage_ratio && (physical || protein_ok || topology_unique)
            }
            Mode::Exploratory => {
                physical
                    || protein_ok
                    || topology_unique
                    || row.dna_lm_delta > 0.25
            }
        };
        if !eligible {
            continue;
        }

        let physical_score = f64::from(link.direct) * 1000.0
            + f64::from(link.gapped) * 420.0
            + f64::from(link.pairs) * 220.0;
        let topology_bonus = if outgoing[link.source].len() == 1 { 220.0 } else { 0.0 }
            + if incoming[link.target].len() == 1 { 220.0 } else { 0.0 };
        let coverage_score = f64::from(coverage_ratio) * 650.0;
        let protein_score = if protein_ok {
            cli.protein_weight
                * f64::from(row.protein_score)
                * f64::from(1.0 - row.protein_ambiguity)
                * f64::from(row.frame_consistency)
        } else if row.protein_score > 0.0 {
            -500.0 * f64::from(row.protein_ambiguity.max(0.25))
        } else {
            0.0
        };
        // Model scores are weak modifiers.  They are never used to create an edge;
        // every candidate already exists in the nucleotide GFA.
        let esm_score = if protein_ok {
            cli.esm_weight * f64::from(row.esm_delta).tanh()
        } else {
            0.0
        };
        let dna_lm_score = if physical || protein_ok {
            cli.dna_lm_weight * f64::from(row.dna_lm_delta).tanh()
        } else {
            0.0
        };
        let gain = segments[link.target]
            .sequence
            .len()
            .saturating_sub(link.overlap) as f64;
        let score = physical_score
            + topology_bonus
            + coverage_score
            + protein_score
            + esm_score
            + dna_lm_score
            + gain.ln_1p() * 18.0;
        ranked.push(RankedLink {
            link,
            score,
            coverage_ratio,
            physical,
            protein_ok,
            evidence: row,
        });
    }

    ranked.sort_unstable_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| right.protein_ok.cmp(&left.protein_ok))
            .then_with(|| right.physical.cmp(&left.physical))
            .then_with(|| right.coverage_ratio.total_cmp(&left.coverage_ratio))
            .then_with(|| left.link.source.cmp(&right.link.source))
            .then_with(|| left.link.target.cmp(&right.link.target))
    });

    let mut successor: FxHashMap<usize, Link> = FxHashMap::default();
    let mut predecessor = vec![None; node_count];
    let mut parent: Vec<usize> = (0..node_count).collect();
    for candidate in &ranked {
        let link = candidate.link;
        if successor.contains_key(&link.source) || predecessor[link.target].is_some() {
            continue;
        }
        let source_root = find(&mut parent, link.source);
        let target_root = find(&mut parent, link.target);
        if source_root == target_root {
            continue;
        }
        successor.insert(link.source, link);
        predecessor[link.target] = Some(link.source);
        union(&mut parent, source_root, target_root);
    }
    (successor, ranked)
}

fn build_paths(node_count: usize, successor: &FxHashMap<usize, Link>) -> Vec<Vec<usize>> {
    let mut predecessor = vec![None; node_count];
    for link in successor.values() {
        predecessor[link.target] = Some(link.source);
    }
    let mut used = vec![false; node_count];
    let mut paths = Vec::new();
    for start in 0..node_count {
        if predecessor[start].is_some() || used[start] {
            continue;
        }
        let path = extend_path(start, successor, &mut used);
        if !path.is_empty() {
            paths.push(path);
        }
    }
    for start in 0..node_count {
        if !used[start] {
            let path = extend_path(start, successor, &mut used);
            if !path.is_empty() {
                paths.push(path);
            }
        }
    }
    paths
}

fn extend_path(
    start: usize,
    successor: &FxHashMap<usize, Link>,
    used: &mut [bool],
) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = start;
    for _ in 0..used.len() {
        if used[current] {
            break;
        }
        used[current] = true;
        path.push(current);
        let Some(link) = successor.get(&current) else {
            break;
        };
        current = link.target;
    }
    path
}

fn assemble_path(
    path: &[usize],
    segments: &[Segment],
    successor: &FxHashMap<usize, Link>,
) -> Vec<u8> {
    let Some((&first, rest)) = path.split_first() else {
        return Vec::new();
    };
    let mut sequence = segments[first].sequence.clone();
    let mut previous = first;
    for &node in rest {
        let overlap = successor
            .get(&previous)
            .filter(|link| link.target == node)
            .map_or(0, |link| link.overlap.min(segments[node].sequence.len()));
        sequence.extend_from_slice(&segments[node].sequence[overlap..]);
        previous = node;
    }
    sequence
}

fn canonical_sequence(sequence: Vec<u8>) -> Vec<u8> {
    let reverse = reverse_complement(&sequence);
    if reverse < sequence {
        reverse
    } else {
        sequence
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn write_fasta(path: &Path, sequences: &[(Vec<usize>, Vec<u8>)]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    for (index, (nodes, sequence)) in sequences.iter().enumerate() {
        writeln!(
            writer,
            ">evidence_path_{:06} len={} unitigs={}",
            index + 1,
            sequence.len(),
            nodes.len()
        )?;
        for chunk in sequence.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_report(
    path: &Path,
    segments: &[Segment],
    ranked: &[RankedLink],
    selected: &FxHashMap<usize, Link>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(
        writer,
        "source\ttarget\tselected\tscore\tdirect\tgapped\tpairs\tcoverage_ratio\tphysical\tprotein_ok\tprotein_score\tprotein_ambiguity\tframe_consistency\tunique_kmers\tprotein_id\tbreakpoint_class\tesm_delta\tdna_lm_delta"
    )?;
    for candidate in ranked {
        let link = candidate.link;
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{:.5}\t{:.5}",
            segments[link.source].name,
            segments[link.target].name,
            selected
                .get(&link.source)
                .is_some_and(|selected_link| selected_link.target == link.target),
            candidate.score,
            link.direct,
            link.gapped,
            link.pairs,
            candidate.coverage_ratio,
            candidate.physical,
            candidate.protein_ok,
            candidate.evidence.protein_score,
            candidate.evidence.protein_ambiguity,
            candidate.evidence.frame_consistency,
            candidate.evidence.unique_kmers,
            candidate.evidence.protein_id,
            candidate.evidence.breakpoint_class,
            candidate.evidence.esm_delta,
            candidate.evidence.dna_lm_delta
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn find(parent: &mut [usize], mut node: usize) -> usize {
    let mut root = node;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[node] != node {
        let next = parent[node];
        parent[node] = root;
        node = next;
    }
    root
}

fn union(parent: &mut [usize], left: usize, right: usize) {
    let left = find(parent, left);
    let right = find(parent, right);
    if left != right {
        parent[right] = left;
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
