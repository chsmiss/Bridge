use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Preset {
    Conservative,
    Megahit,
    Aggressive,
}

#[derive(Parser, Debug)]
#[command(
    name = "bridgeasm-cleanpath",
    about = "Iterative consensus graph cleaning and recompaction over BridgeAsm GFA"
)]
struct Cli {
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    removed: Option<PathBuf>,
    #[arg(long)]
    report: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = Preset::Megahit)]
    preset: Preset,
    #[arg(long, default_value_t = 5)]
    rounds: usize,
    #[arg(long, default_value_t = 2)]
    min_direct: u32,
    #[arg(long, default_value_t = 2)]
    min_pair: u32,
    #[arg(long, default_value_t = 200)]
    max_tip_len: usize,
    #[arg(long, default_value_t = 0.10)]
    disconnect_ratio: f32,
    #[arg(long, default_value_t = 0.20)]
    low_local_ratio: f32,
    #[arg(long, default_value_t = 0.60)]
    bubble_dominance: f32,
    #[arg(long, default_value_t = 200)]
    min_length: usize,
}

#[derive(Clone, Debug)]
struct Segment {
    name: String,
    sequence: Vec<u8>,
    coverage: f32,
    active: bool,
    removed_reason: Option<&'static str>,
}

#[derive(Clone, Copy, Debug)]
struct Link {
    source: usize,
    target: usize,
    overlap: usize,
    direct: u32,
    gapped: u32,
    pairs: u32,
    active: bool,
}

#[derive(Default)]
struct RoundStats {
    tips: usize,
    weak_edges: usize,
    bubble_nodes: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate(&cli)?;
    let (mut segments, mut links) = read_gfa(&cli.gfa)?;
    if segments.is_empty() {
        bail!("GFA contains no segments");
    }

    let mut round_stats = Vec::new();
    for _ in 0..cli.rounds {
        let mut stats = RoundStats::default();
        stats.tips += remove_weak_tips(&mut segments, &mut links, &cli);
        stats.weak_edges += remove_low_local_edges(&segments, &mut links, &cli);
        stats.bubble_nodes += pop_simple_bubbles(&mut segments, &mut links, &cli);
        deactivate_incident_links(&segments, &mut links);
        let changed = stats.tips + stats.weak_edges + stats.bubble_nodes;
        round_stats.push(stats);
        if changed == 0 {
            break;
        }
    }

    let paths = compact_clean_graph(&segments, &links);
    let mut records = paths
        .into_iter()
        .map(|path| {
            let sequence = assemble_path(&path, &segments, &links);
            (path, canonical_sequence(sequence))
        })
        .filter(|(_, sequence)| sequence.len() >= cli.min_length)
        .collect::<Vec<_>>();
    exact_dedup(&mut records);
    records.sort_unstable_by(|left, right| {
        right
            .1
            .len()
            .cmp(&left.1.len())
            .then_with(|| left.1.cmp(&right.1))
    });
    write_fasta(&cli.output, "cleanpath", &records)?;

    if let Some(path) = cli.removed.as_ref() {
        let removed_records = segments
            .iter()
            .enumerate()
            .filter(|(_, segment)| !segment.active && segment.sequence.len() >= cli.min_length)
            .map(|(index, segment)| (vec![index], canonical_sequence(segment.sequence.clone())))
            .collect::<Vec<_>>();
        write_fasta(path, "removed", &removed_records)?;
    }
    if let Some(path) = cli.report.as_ref() {
        write_report(path, &segments, &links, &round_stats)?;
    }

    let mut lengths = records.iter().map(|(_, sequence)| sequence.len()).collect::<Vec<_>>();
    let total: usize = lengths.iter().sum();
    let n50 = n50(&mut lengths);
    let largest = lengths.iter().copied().max().unwrap_or(0);
    let active_nodes = segments.iter().filter(|segment| segment.active).count();
    let active_edges = links.iter().filter(|link| link.active).count();
    eprintln!(
        "cleanpath preset={:?}: contigs={} total_bp={} N50={} largest={} active_nodes={} active_edges={} rounds={}",
        cli.preset,
        records.len(),
        total,
        n50,
        largest,
        active_nodes,
        active_edges,
        round_stats.len()
    );
    Ok(())
}

fn validate(cli: &Cli) -> Result<()> {
    if cli.rounds == 0 {
        bail!("rounds must be positive");
    }
    for (name, value) in [
        ("disconnect-ratio", cli.disconnect_ratio),
        ("low-local-ratio", cli.low_local_ratio),
        ("bubble-dominance", cli.bubble_dominance),
    ] {
        if !(0.0..=1.0).contains(&value) {
            bail!("{name} must be in 0..=1");
        }
    }
    Ok(())
}

fn effective_ratios(cli: &Cli) -> (f32, f32, f32) {
    match cli.preset {
        Preset::Conservative => (
            cli.disconnect_ratio * 0.5,
            cli.low_local_ratio * 0.5,
            cli.bubble_dominance.max(0.70),
        ),
        Preset::Megahit => (
            cli.disconnect_ratio,
            cli.low_local_ratio,
            cli.bubble_dominance,
        ),
        Preset::Aggressive => (
            (cli.disconnect_ratio * 2.0).min(0.35),
            (cli.low_local_ratio * 1.5).min(0.50),
            (cli.bubble_dominance * 0.85).max(0.50),
        ),
    }
}

fn read_gfa(path: &PathBuf) -> Result<(Vec<Segment>, Vec<Link>)> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut segments = Vec::new();
    let mut name_to_id = FxHashMap::default();
    let mut raw_links = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('H') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                let name = fields[1].to_string();
                let id = segments.len();
                name_to_id.insert(name.clone(), id);
                segments.push(Segment {
                    name,
                    sequence: fields[2].as_bytes().to_vec(),
                    coverage: tag_f32(&fields[3..], "KC:f:").unwrap_or(1.0),
                    active: true,
                    removed_reason: None,
                });
            }
            Some("L") if fields.len() >= 6 => raw_links.push(line),
            _ => {}
        }
    }

    let mut links = Vec::with_capacity(raw_links.len());
    for line in raw_links {
        let fields = line.split('\t').collect::<Vec<_>>();
        let (Some(&source), Some(&target)) =
            (name_to_id.get(fields[1]), name_to_id.get(fields[3]))
        else {
            continue;
        };
        if source == target {
            continue;
        }
        links.push(Link {
            source,
            target,
            overlap: fields[5]
                .strip_suffix('M')
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0),
            direct: tag_u32(&fields[6..], "DR:i:").unwrap_or(0),
            gapped: tag_u32(&fields[6..], "GR:i:").unwrap_or(0),
            pairs: tag_u32(&fields[6..], "PE:i:").unwrap_or(0),
            active: true,
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

fn adjacency(segments: &[Segment], links: &[Link]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut outgoing = vec![Vec::new(); segments.len()];
    let mut incoming = vec![Vec::new(); segments.len()];
    for (edge_id, link) in links.iter().enumerate() {
        if link.active && segments[link.source].active && segments[link.target].active {
            outgoing[link.source].push(edge_id);
            incoming[link.target].push(edge_id);
        }
    }
    (outgoing, incoming)
}

fn physical(link: &Link, cli: &Cli) -> bool {
    link.direct >= cli.min_direct || link.pairs >= cli.min_pair
}

fn remove_weak_tips(segments: &mut [Segment], links: &mut [Link], cli: &Cli) -> usize {
    let (disconnect, _, _) = effective_ratios(cli);
    let (outgoing, incoming) = adjacency(segments, links);
    let mut remove = Vec::new();
    for node in 0..segments.len() {
        if !segments[node].active || segments[node].sequence.len() > cli.max_tip_len {
            continue;
        }
        let terminal = incoming[node].is_empty() || outgoing[node].is_empty();
        if !terminal {
            continue;
        }
        let incident = incoming[node]
            .iter()
            .chain(outgoing[node].iter())
            .copied()
            .collect::<Vec<_>>();
        if incident.iter().any(|&edge_id| physical(&links[edge_id], cli)) {
            continue;
        }
        let neighbor_cov = incident
            .iter()
            .map(|&edge_id| {
                let link = links[edge_id];
                let other = if link.source == node { link.target } else { link.source };
                segments[other].coverage
            })
            .fold(0.0_f32, f32::max);
        if neighbor_cov > 0.0 && segments[node].coverage <= neighbor_cov * disconnect {
            remove.push(node);
        }
    }
    for node in &remove {
        segments[*node].active = false;
        segments[*node].removed_reason = Some("tip");
    }
    deactivate_incident_links(segments, links);
    remove.len()
}

fn remove_low_local_edges(segments: &[Segment], links: &mut [Link], cli: &Cli) -> usize {
    let (disconnect, low_local, _) = effective_ratios(cli);
    let (outgoing, incoming) = adjacency(segments, links);
    let mut remove = Vec::new();
    for (edge_id, link) in links.iter().enumerate() {
        if !link.active || physical(link, cli) {
            continue;
        }
        let source_branch = outgoing[link.source].len() > 1;
        let target_branch = incoming[link.target].len() > 1;
        if !source_branch && !target_branch {
            continue;
        }
        let source_cov = segments[link.source].coverage.max(0.001);
        let target_cov = segments[link.target].coverage.max(0.001);
        let local_ratio = source_cov.min(target_cov) / source_cov.max(target_cov);
        let best_out = outgoing[link.source]
            .iter()
            .map(|&id| segments[links[id].target].coverage)
            .fold(0.0_f32, f32::max)
            .max(0.001);
        let best_in = incoming[link.target]
            .iter()
            .map(|&id| segments[links[id].source].coverage)
            .fold(0.0_f32, f32::max)
            .max(0.001);
        let outgoing_fraction = target_cov / best_out;
        let incoming_fraction = source_cov / best_in;
        if local_ratio < low_local
            || (source_branch && outgoing_fraction < disconnect)
            || (target_branch && incoming_fraction < disconnect)
        {
            remove.push(edge_id);
        }
    }
    for edge_id in &remove {
        links[*edge_id].active = false;
    }
    remove.len()
}

fn pop_simple_bubbles(segments: &mut [Segment], links: &mut [Link], cli: &Cli) -> usize {
    let (_, _, dominance) = effective_ratios(cli);
    let (outgoing, incoming) = adjacency(segments, links);
    let mut remove_nodes = FxHashSet::default();

    for source in 0..segments.len() {
        if !segments[source].active || outgoing[source].len() < 2 || outgoing[source].len() > 8 {
            continue;
        }
        let mut by_sink: FxHashMap<usize, Vec<(usize, usize, usize)>> = FxHashMap::default();
        for &first_edge_id in &outgoing[source] {
            let first = links[first_edge_id];
            let middle = first.target;
            if !segments[middle].active || incoming[middle].len() != 1 || outgoing[middle].len() != 1 {
                continue;
            }
            let second_edge_id = outgoing[middle][0];
            let second = links[second_edge_id];
            if second.target == source {
                continue;
            }
            by_sink
                .entry(second.target)
                .or_default()
                .push((middle, first_edge_id, second_edge_id));
        }
        for candidates in by_sink.values() {
            if candidates.len() < 2 {
                continue;
            }
            let mut scored = candidates
                .iter()
                .map(|&(middle, left_edge, right_edge)| {
                    let left = links[left_edge];
                    let right = links[right_edge];
                    let support = left.direct.min(right.direct) as f32
                        + 0.5 * left.pairs.min(right.pairs) as f32
                        + 0.2 * left.gapped.min(right.gapped) as f32;
                    let score = segments[middle].coverage.max(0.001) * (1.0 + support);
                    (middle, score, physical(&left, cli) && physical(&right, cli))
                })
                .collect::<Vec<_>>();
            scored.sort_unstable_by(|left, right| right.1.total_cmp(&left.1));
            let best = scored[0];
            let second_score = scored.get(1).map_or(0.0, |entry| entry.1);
            if best.1 <= 0.0 || (best.1 - second_score) / best.1 < dominance {
                continue;
            }
            for &(middle, score, both_physical) in scored.iter().skip(1) {
                if both_physical {
                    continue;
                }
                if score <= best.1 * (1.0 - dominance) {
                    remove_nodes.insert(middle);
                }
            }
        }
    }

    for &node in &remove_nodes {
        segments[node].active = false;
        segments[node].removed_reason = Some("bubble");
    }
    deactivate_incident_links(segments, links);
    remove_nodes.len()
}

fn deactivate_incident_links(segments: &[Segment], links: &mut [Link]) {
    for link in links {
        if !segments[link.source].active || !segments[link.target].active {
            link.active = false;
        }
    }
}

fn compact_clean_graph(segments: &[Segment], links: &[Link]) -> Vec<Vec<usize>> {
    let (outgoing, incoming) = adjacency(segments, links);
    let mut used = vec![false; segments.len()];
    let mut paths = Vec::new();

    for start in 0..segments.len() {
        if !segments[start].active || used[start] {
            continue;
        }
        let indeg = incoming[start].len();
        let outdeg = outgoing[start].len();
        if indeg == 1 && outdeg == 1 {
            continue;
        }
        if outdeg == 0 {
            used[start] = true;
            paths.push(vec![start]);
            continue;
        }
        for &edge_id in &outgoing[start] {
            if used[start] && outgoing[start].len() == 1 {
                break;
            }
            let mut path = vec![start];
            let mut current = links[edge_id].target;
            let mut local_seen = FxHashSet::default();
            local_seen.insert(start);
            while segments[current].active && !local_seen.contains(&current) {
                path.push(current);
                local_seen.insert(current);
                if incoming[current].len() != 1 || outgoing[current].len() != 1 {
                    break;
                }
                current = links[outgoing[current][0]].target;
            }
            for &node in &path {
                if incoming[node].len() <= 1 && outgoing[node].len() <= 1 {
                    used[node] = true;
                }
            }
            paths.push(path);
        }
    }

    for node in 0..segments.len() {
        if segments[node].active && !used[node] {
            paths.push(vec![node]);
        }
    }
    paths
}

fn assemble_path(path: &[usize], segments: &[Segment], links: &[Link]) -> Vec<u8> {
    let Some((&first, rest)) = path.split_first() else {
        return Vec::new();
    };
    let mut sequence = segments[first].sequence.clone();
    let mut previous = first;
    for &node in rest {
        let overlap = links
            .iter()
            .find(|link| link.active && link.source == previous && link.target == node)
            .map_or(0, |link| link.overlap.min(segments[node].sequence.len()));
        sequence.extend_from_slice(&segments[node].sequence[overlap..]);
        previous = node;
    }
    sequence
}

fn canonical_sequence(sequence: Vec<u8>) -> Vec<u8> {
    let reverse = reverse_complement(&sequence);
    if reverse < sequence { reverse } else { sequence }
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

fn exact_dedup(records: &mut Vec<(Vec<usize>, Vec<u8>)>) {
    let mut seen = FxHashSet::default();
    records.retain(|(_, sequence)| seen.insert(sequence.clone()));
}

fn write_fasta(path: &PathBuf, prefix: &str, records: &[(Vec<usize>, Vec<u8>)]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    for (index, (nodes, sequence)) in records.iter().enumerate() {
        writeln!(
            writer,
            ">{prefix}_{:06} len={} unitigs={}",
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
    links: &[Link],
    rounds: &[RoundStats],
) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(writer, "round\ttips\tweak_edges\tbubble_nodes")?;
    for (index, stats) in rounds.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            index + 1,
            stats.tips,
            stats.weak_edges,
            stats.bubble_nodes
        )?;
    }
    writeln!(writer, "\nremoved_segment\treason\tlength\tcoverage")?;
    for segment in segments.iter().filter(|segment| !segment.active) {
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.3}",
            segment.name,
            segment.removed_reason.unwrap_or("incident"),
            segment.sequence.len(),
            segment.coverage
        )?;
    }
    writeln!(
        writer,
        "\nactive_nodes\t{}\nactive_edges\t{}",
        segments.iter().filter(|segment| segment.active).count(),
        links.iter().filter(|link| link.active).count()
    )?;
    writer.flush()?;
    Ok(())
}

fn n50(lengths: &mut [usize]) -> usize {
    if lengths.is_empty() {
        return 0;
    }
    lengths.sort_unstable_by(|left, right| right.cmp(left));
    let total: usize = lengths.iter().sum();
    let mut cumulative = 0;
    for &length in lengths.iter() {
        cumulative += length;
        if cumulative * 2 >= total {
            return length;
        }
    }
    0
}
