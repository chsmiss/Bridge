use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Conservative,
    Megahit,
    Aggressive,
    Coverage,
}

#[derive(Parser, Debug)]
#[command(
    name = "bridgeasm-flowpath",
    about = "Coverage-flow path cover over an existing BridgeAsm UnitigGraph GFA"
)]
struct Cli {
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    report: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = Mode::Megahit)]
    mode: Mode,
    #[arg(long, default_value_t = 2)]
    min_direct: u32,
    #[arg(long, default_value_t = 2)]
    min_pair: u32,
    #[arg(long, default_value_t = 0.10)]
    disconnect_ratio: f32,
    #[arg(long, default_value_t = 0.20)]
    low_local_ratio: f32,
    #[arg(long, default_value_t = 0.10)]
    sibling_ratio: f32,
    #[arg(long, default_value_t = 200)]
    max_tip_len: usize,
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

#[derive(Clone, Copy, Debug)]
struct RankedLink {
    link: Link,
    score: f64,
    coverage_ratio: f32,
    source_fraction: f32,
    target_fraction: f32,
    physical: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate(&cli)?;
    let (segments, links) = read_gfa(&cli.gfa)?;
    if segments.is_empty() {
        bail!("GFA contains no segments");
    }

    let (selected, ranked) = select_links(&segments, &links, &cli);
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

    let mut lengths: Vec<usize> = sequences
        .iter()
        .map(|(_, sequence)| sequence.len())
        .collect();
    let total: usize = lengths.iter().sum();
    let n50 = n50(&mut lengths);
    let largest = lengths.iter().copied().max().unwrap_or(0);
    eprintln!(
        "flowpath mode={:?}: {} contigs, {} bp, N50 {}, largest {}, selected_edges={}",
        cli.mode,
        sequences.len(),
        total,
        n50,
        largest,
        selected.len()
    );
    Ok(())
}

fn validate(cli: &Cli) -> Result<()> {
    for (name, value) in [
        ("disconnect-ratio", cli.disconnect_ratio),
        ("low-local-ratio", cli.low_local_ratio),
        ("sibling-ratio", cli.sibling_ratio),
    ] {
        if !(0.0..=1.0).contains(&value) {
            bail!("{name} must be in 0..=1");
        }
    }
    Ok(())
}

fn read_gfa(path: &PathBuf) -> Result<(Vec<Segment>, Vec<Link>)> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut segments = Vec::new();
    let mut name_to_id = FxHashMap::default();
    let mut raw_links: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('H') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                let name = fields[1].to_string();
                let sequence = fields[2].as_bytes().to_vec();
                let coverage = tag_f32(&fields[3..], "KC:f:").unwrap_or(1.0);
                let id = segments.len();
                name_to_id.insert(name.clone(), id);
                segments.push(Segment {
                    name,
                    sequence,
                    coverage,
                });
            }
            Some("L") if fields.len() >= 6 => raw_links.push(line),
            _ => {}
        }
    }

    let mut links = Vec::with_capacity(raw_links.len());
    for line in raw_links {
        let fields: Vec<&str> = line.split('\t').collect();
        let Some(&source) = name_to_id.get(fields[1]) else {
            continue;
        };
        let Some(&target) = name_to_id.get(fields[3]) else {
            continue;
        };
        if source == target {
            continue;
        }
        let overlap = fields[5]
            .strip_suffix('M')
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(0);
        links.push(Link {
            source,
            target,
            overlap,
            direct: tag_u32(&fields[6..], "DR:i:").unwrap_or(0),
            gapped: tag_u32(&fields[6..], "GR:i:").unwrap_or(0),
            pairs: tag_u32(&fields[6..], "PE:i:").unwrap_or(0),
        });
    }
    Ok((segments, links))
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

fn select_links(
    segments: &[Segment],
    links: &[Link],
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
        let source_cov = segments[link.source].coverage.max(0.001);
        let target_cov = segments[link.target].coverage.max(0.001);
        let coverage_ratio = source_cov.min(target_cov) / source_cov.max(target_cov);
        let best_out_cov = outgoing[link.source]
            .iter()
            .map(|edge| segments[edge.target].coverage)
            .fold(0.0_f32, f32::max)
            .max(0.001);
        let best_in_cov = incoming[link.target]
            .iter()
            .map(|edge| segments[edge.source].coverage)
            .fold(0.0_f32, f32::max)
            .max(0.001);
        let source_fraction = target_cov / best_out_cov;
        let target_fraction = source_cov / best_in_cov;
        let physical = link.direct >= cli.min_direct || link.pairs >= cli.min_pair;

        let source_is_tip = incoming[link.source].is_empty()
            && segments[link.source].sequence.len() <= cli.max_tip_len;
        let target_is_tip = outgoing[link.target].is_empty()
            && segments[link.target].sequence.len() <= cli.max_tip_len;
        let weak_tip = (source_is_tip && source_cov < target_cov * cli.disconnect_ratio)
            || (target_is_tip && target_cov < source_cov * cli.disconnect_ratio);
        if weak_tip && !physical {
            continue;
        }

        let eligible = match cli.mode {
            Mode::Conservative => {
                physical
                    && coverage_ratio >= cli.low_local_ratio
                    && source_fraction >= cli.sibling_ratio
                    && target_fraction >= cli.sibling_ratio
            }
            Mode::Megahit => {
                physical
                    || (coverage_ratio >= cli.low_local_ratio
                        && source_fraction >= cli.disconnect_ratio
                        && target_fraction >= cli.disconnect_ratio)
            }
            Mode::Aggressive => {
                physical
                    || (coverage_ratio >= cli.disconnect_ratio
                        && source_fraction >= cli.sibling_ratio
                        && target_fraction >= cli.sibling_ratio)
            }
            Mode::Coverage => {
                coverage_ratio >= cli.disconnect_ratio
                    && source_fraction >= cli.sibling_ratio
                    && target_fraction >= cli.sibling_ratio
            }
        };
        if !eligible {
            continue;
        }

        let physical_score = f64::from(link.direct) * 1000.0
            + f64::from(link.gapped) * 350.0
            + f64::from(link.pairs) * 180.0;
        let topology_bonus = if outgoing[link.source].len() == 1 {
            250.0
        } else {
            0.0
        } + if incoming[link.target].len() == 1 {
            250.0
        } else {
            0.0
        };
        let coverage_score = f64::from(coverage_ratio) * 600.0
            + f64::from(source_fraction.min(1.0)) * 250.0
            + f64::from(target_fraction.min(1.0)) * 250.0;
        let gain = segments[link.target]
            .sequence
            .len()
            .saturating_sub(link.overlap) as f64;
        let score = physical_score + topology_bonus + coverage_score + gain.ln_1p() * 20.0;
        ranked.push(RankedLink {
            link,
            score,
            coverage_ratio,
            source_fraction,
            target_fraction,
            physical,
        });
    }

    ranked.sort_unstable_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(Ordering::Equal)
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

fn extend_path(start: usize, successor: &FxHashMap<usize, Link>, used: &mut [bool]) -> Vec<usize> {
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

fn write_fasta(path: &PathBuf, sequences: &[(Vec<usize>, Vec<u8>)]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    for (index, (nodes, sequence)) in sequences.iter().enumerate() {
        writeln!(
            writer,
            ">flowpath_{:06} len={} unitigs={}",
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
    path: &PathBuf,
    segments: &[Segment],
    ranked: &[RankedLink],
    selected: &FxHashMap<usize, Link>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(
        writer,
        "source\ttarget\tselected\tscore\tdirect\tgapped\tpairs\tcoverage_ratio\tsource_fraction\ttarget_fraction\tphysical"
    )?;
    for candidate in ranked {
        let link = candidate.link;
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}",
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
            candidate.source_fraction,
            candidate.target_fraction,
            candidate.physical
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
