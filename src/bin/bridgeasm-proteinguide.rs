use anyhow::{bail, Context, Result};
use bridgeasm::dna::reverse_complement;
use clap::{Parser, ValueEnum};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Mode {
    /// Only join anchors that have a validated suffix-prefix overlap.
    OverlapOnly,
    /// Permit only short, nearly full-length, very high-identity guide gaps.
    Conservative,
    /// Permit longer protein-guided gaps, while retaining all provenance.
    Protein,
}

#[derive(Parser, Debug)]
#[command(
    name = "bridgeasm-proteinguide",
    about = "Build an auditable path cover from PenguiN/PLASS-guided nucleotide contigs"
)]
struct Cli {
    /// BridgeAsm UnitigGraph GFA.
    #[arg(long)]
    gfa: PathBuf,
    /// PenguiN guided_nuclassemble nucleotide FASTA (or another protein-guided FASTA).
    #[arg(long)]
    guide: PathBuf,
    /// PAF produced by mapping GFA segments (queries) to guide contigs (targets).
    #[arg(long)]
    paf: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    report: Option<PathBuf>,
    #[arg(long)]
    summary: Option<PathBuf>,
    /// Optional model evidence TSV: source, target, score, decision, scorer.
    /// Oriented IDs use forms such as u12+ and u57-; decision can be neutral or veto.
    #[arg(long)]
    model_scores: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = Mode::Conservative)]
    mode: Mode,
    #[arg(long, default_value_t = 30)]
    min_mapq: u8,
    #[arg(long, default_value_t = 90)]
    min_alignment: usize,
    #[arg(long, default_value_t = 0.80)]
    min_query_fraction: f32,
    #[arg(long, default_value_t = 0.97)]
    min_identity: f32,
    /// Reject a unitig when a competing alignment has at least this fraction of the best score.
    #[arg(long, default_value_t = 0.90)]
    secondary_ratio: f32,
    #[arg(long, default_value_t = 0.15)]
    min_coverage_ratio: f32,
    #[arg(long, default_value_t = 300)]
    max_gap: usize,
    #[arg(long, default_value_t = 90)]
    conservative_max_gap: usize,
    #[arg(long, default_value_t = 2_000)]
    max_overlap: usize,
    #[arg(long, default_value_t = 15)]
    min_overlap: usize,
    #[arg(long, default_value_t = 0.97)]
    min_overlap_identity: f32,
    #[arg(long, default_value_t = 32)]
    overlap_slack: usize,
    #[arg(long, default_value_t = 30)]
    max_gap_disagreement: usize,
    #[arg(long, default_value_t = 0.05)]
    max_bridge_n_fraction: f32,
    #[arg(long, default_value_t = 1)]
    min_guide_support: usize,
    /// Weight for a normalized external model score. Models rank or veto; they never create bases.
    #[arg(long, default_value_t = 400.0)]
    model_weight: f64,
    /// When set, reject model-scored candidates below this value. Unscored candidates remain eligible.
    #[arg(long)]
    min_model_score: Option<f64>,
    #[arg(long, default_value_t = 200)]
    min_length: usize,
}

#[derive(Clone, Debug)]
struct Segment {
    name: String,
    sequence: Vec<u8>,
    coverage: f32,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct OrientedNode {
    id: usize,
    reverse: bool,
}

#[derive(Debug)]
struct GfaGraph {
    segments: Vec<Segment>,
    name_to_id: FxHashMap<String, usize>,
    links: FxHashMap<(OrientedNode, OrientedNode), usize>,
}

#[derive(Clone, Debug)]
struct Anchor {
    segment_id: usize,
    guide_name: String,
    reverse: bool,
    full_start: i64,
    full_end: i64,
    identity: f32,
    query_fraction: f32,
    mapq: u8,
    aligned_length: usize,
    rank_score: f64,
    primary: bool,
}

#[derive(Default, Debug, Serialize)]
struct AnchorStats {
    paf_records: usize,
    malformed_or_unknown: usize,
    threshold_filtered: usize,
    low_mapq_or_non_primary: usize,
    ambiguous: usize,
    accepted: usize,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct EdgeKey {
    source: OrientedNode,
    target: OrientedNode,
}

#[derive(Clone, Debug)]
struct Candidate {
    key: EdgeKey,
    guide_name: String,
    projected_gap: i64,
    overlap: usize,
    bridge: Vec<u8>,
    direct_graph: bool,
    identity: f32,
    query_fraction: f32,
    mapq: u8,
    aligned_length: usize,
    coverage_ratio: f32,
    base_score: f64,
}

#[derive(Debug)]
struct CandidateAggregate {
    best: Candidate,
    guides: FxHashSet<String>,
    min_gap: i64,
    max_gap: i64,
    conflicting_gap: bool,
}

#[derive(Clone, Debug, Default)]
struct ModelEvidence {
    score_sum: f64,
    votes: usize,
    veto: bool,
    scorers: FxHashSet<String>,
}

#[derive(Clone, Debug)]
struct RankedCandidate {
    best: Candidate,
    guide_support: usize,
    model_score: Option<f64>,
    model_votes: usize,
    model_scorers: usize,
    score: f64,
    eligible: bool,
    reason: String,
    selected: bool,
}

#[derive(Clone, Debug)]
struct OutputRecord {
    sequence: Vec<u8>,
    unitigs: usize,
    guide_edges: usize,
    guide_bases: usize,
}

#[derive(Debug, Serialize)]
struct RunSummary {
    mode: String,
    segments: usize,
    guide_contigs: usize,
    anchor_stats: AnchorStats,
    candidate_edges: usize,
    eligible_edges: usize,
    selected_edges: usize,
    model_scored_candidates: usize,
    model_vetoes: usize,
    output_contigs: usize,
    output_bases: usize,
    output_n50: usize,
    largest_contig: usize,
    inserted_guide_bases: usize,
    rejection_counts: BTreeMap<String, usize>,
}

#[derive(Clone, Copy, Debug)]
struct SelectionConfig {
    mode: Mode,
    min_coverage_ratio: f32,
    max_gap: usize,
    conservative_max_gap: usize,
    max_overlap: usize,
    min_overlap: usize,
    min_overlap_identity: f32,
    overlap_slack: usize,
    max_gap_disagreement: usize,
    max_bridge_n_fraction: f32,
    min_guide_support: usize,
    model_weight: f64,
    min_model_score: Option<f64>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate_cli(&cli)?;

    let graph = read_gfa(&cli.gfa)?;
    if graph.segments.is_empty() {
        bail!("GFA contains no segments");
    }
    let guides = read_fasta(&cli.guide)?;
    if guides.is_empty() {
        bail!("guide FASTA contains no sequences");
    }

    let (anchors_by_guide, anchor_stats) = read_unique_anchors(&cli.paf, &graph, &guides, &cli)?;
    let model_evidence = if let Some(path) = cli.model_scores.as_deref() {
        read_model_scores(path, &graph.name_to_id)?
    } else {
        FxHashMap::default()
    };
    let config = SelectionConfig::from(&cli);
    let (mut ranked, rejection_counts) = build_candidates(
        &graph,
        &guides,
        &anchors_by_guide,
        &model_evidence,
        config,
    );
    let (successor, predecessor, orientation) =
        select_path_cover(&mut ranked, graph.segments.len());
    let records = build_output_records(
        &graph.segments,
        &successor,
        &predecessor,
        &orientation,
        cli.min_length,
    );
    write_fasta(&cli.output, &records)?;
    if let Some(path) = cli.report.as_deref() {
        write_report(path, &ranked, &graph.segments)?;
    }

    let mut lengths: Vec<usize> = records.iter().map(|record| record.sequence.len()).collect();
    let output_bases: usize = lengths.iter().sum();
    let output_n50 = n50(&mut lengths);
    let largest_contig = lengths.iter().copied().max().unwrap_or(0);
    let selected_edges = ranked.iter().filter(|candidate| candidate.selected).count();
    let inserted_guide_bases: usize = records.iter().map(|record| record.guide_bases).sum();
    let summary = RunSummary {
        mode: format!("{:?}", cli.mode),
        segments: graph.segments.len(),
        guide_contigs: guides.len(),
        anchor_stats,
        candidate_edges: ranked.len(),
        eligible_edges: ranked.iter().filter(|candidate| candidate.eligible).count(),
        selected_edges,
        model_scored_candidates: ranked
            .iter()
            .filter(|candidate| candidate.model_score.is_some())
            .count(),
        model_vetoes: ranked
            .iter()
            .filter(|candidate| candidate.reason == "model_veto")
            .count(),
        output_contigs: records.len(),
        output_bases,
        output_n50,
        largest_contig,
        inserted_guide_bases,
        rejection_counts: rejection_counts
            .into_iter()
            .map(|(reason, count)| (reason.to_string(), count))
            .collect(),
    };
    if let Some(path) = cli.summary.as_deref() {
        let file = File::create(path)
            .with_context(|| format!("failed to create {}", path.display()))?;
        serde_json::to_writer_pretty(BufWriter::new(file), &summary)
            .context("failed to write protein-guide summary")?;
    }

    eprintln!(
        "proteinguide mode={:?}: anchors={}, candidates={} (eligible={}), selected_edges={}, contigs={}, bases={}, N50={}, largest={}, inserted_guide_bases={}",
        cli.mode,
        summary.anchor_stats.accepted,
        summary.candidate_edges,
        summary.eligible_edges,
        summary.selected_edges,
        summary.output_contigs,
        summary.output_bases,
        summary.output_n50,
        summary.largest_contig,
        summary.inserted_guide_bases
    );
    Ok(())
}

impl From<&Cli> for SelectionConfig {
    fn from(cli: &Cli) -> Self {
        Self {
            mode: cli.mode,
            min_coverage_ratio: cli.min_coverage_ratio,
            max_gap: cli.max_gap,
            conservative_max_gap: cli.conservative_max_gap,
            max_overlap: cli.max_overlap,
            min_overlap: cli.min_overlap,
            min_overlap_identity: cli.min_overlap_identity,
            overlap_slack: cli.overlap_slack,
            max_gap_disagreement: cli.max_gap_disagreement,
            max_bridge_n_fraction: cli.max_bridge_n_fraction,
            min_guide_support: cli.min_guide_support,
            model_weight: cli.model_weight,
            min_model_score: cli.min_model_score,
        }
    }
}

fn validate_cli(cli: &Cli) -> Result<()> {
    for (name, value) in [
        ("min-query-fraction", cli.min_query_fraction),
        ("min-identity", cli.min_identity),
        ("secondary-ratio", cli.secondary_ratio),
        ("min-coverage-ratio", cli.min_coverage_ratio),
        ("min-overlap-identity", cli.min_overlap_identity),
        ("max-bridge-n-fraction", cli.max_bridge_n_fraction),
    ] {
        if !(0.0..=1.0).contains(&value) {
            bail!("{name} must be in 0..=1");
        }
    }
    if cli.min_overlap > cli.max_overlap {
        bail!("min-overlap cannot exceed max-overlap");
    }
    if cli.min_guide_support == 0 {
        bail!("min-guide-support must be at least one");
    }
    if let Some(score) = cli.min_model_score {
        if !score.is_finite() {
            bail!("min-model-score must be finite");
        }
    }
    Ok(())
}

fn read_gfa(path: &Path) -> Result<GfaGraph> {
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
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                if fields[2] == "*" {
                    bail!("GFA segment {} has no sequence", fields[1]);
                }
                let name = fields[1].to_string();
                if name_to_id.contains_key(&name) {
                    bail!("duplicate GFA segment name: {name}");
                }
                let coverage = tag_f32(&fields[3..], "KC:f:")
                    .or_else(|| tag_f32(&fields[3..], "dp:f:"))
                    .unwrap_or(1.0);
                let id = segments.len();
                name_to_id.insert(name.clone(), id);
                segments.push(Segment {
                    name,
                    sequence: fields[2].as_bytes().to_ascii_uppercase(),
                    coverage,
                });
            }
            Some("L") if fields.len() >= 6 => raw_links.push((
                fields[1].to_string(),
                fields[2] == "-",
                fields[3].to_string(),
                fields[4] == "-",
                parse_overlap(fields[5]),
            )),
            _ => {}
        }
    }

    let mut links = FxHashMap::default();
    for (source_name, source_reverse, target_name, target_reverse, overlap) in raw_links {
        let Some(&source_id) = name_to_id.get(&source_name) else {
            continue;
        };
        let Some(&target_id) = name_to_id.get(&target_name) else {
            continue;
        };
        let source = OrientedNode {
            id: source_id,
            reverse: source_reverse,
        };
        let target = OrientedNode {
            id: target_id,
            reverse: target_reverse,
        };
        links.insert((source, target), overlap);
        links.insert(
            (
                OrientedNode {
                    id: target_id,
                    reverse: !target_reverse,
                },
                OrientedNode {
                    id: source_id,
                    reverse: !source_reverse,
                },
            ),
            overlap,
        );
    }

    Ok(GfaGraph {
        segments,
        name_to_id,
        links,
    })
}

fn parse_overlap(text: &str) -> usize {
    text.bytes()
        .take_while(u8::is_ascii_digit)
        .fold(0_usize, |value, digit| {
            value
                .saturating_mul(10)
                .saturating_add(usize::from(digit - b'0'))
        })
}

fn tag_f32(fields: &[&str], prefix: &str) -> Option<f32> {
    fields
        .iter()
        .find_map(|field| field.strip_prefix(prefix)?.parse().ok())
}

fn read_fasta(path: &Path) -> Result<FxHashMap<String, Vec<u8>>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut records = FxHashMap::default();
    let mut name: Option<String> = None;
    let mut sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = name.take() {
                if records.insert(previous.clone(), sequence).is_some() {
                    bail!("duplicate FASTA record: {previous}");
                }
                sequence = Vec::new();
            }
            let record_name = header
                .split_whitespace()
                .next()
                .filter(|value| !value.is_empty())
                .context("empty FASTA header")?;
            name = Some(record_name.to_string());
        } else if !line.is_empty() {
            sequence.extend(line.bytes().map(|base| base.to_ascii_uppercase()));
        }
    }
    if let Some(previous) = name {
        if records.insert(previous.clone(), sequence).is_some() {
            bail!("duplicate FASTA record: {previous}");
        }
    }
    Ok(records)
}

fn read_unique_anchors(
    path: &Path,
    graph: &GfaGraph,
    guides: &FxHashMap<String, Vec<u8>>,
    cli: &Cli,
) -> Result<(FxHashMap<String, Vec<Anchor>>, AnchorStats)> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut stats = AnchorStats::default();
    let mut by_segment: Vec<Vec<Anchor>> = vec![Vec::new(); graph.segments.len()];

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        stats.paf_records += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            stats.malformed_or_unknown += 1;
            continue;
        }
        let Some(&segment_id) = graph.name_to_id.get(fields[0]) else {
            stats.malformed_or_unknown += 1;
            continue;
        };
        let Some(guide_sequence) = guides.get(fields[5]) else {
            stats.malformed_or_unknown += 1;
            continue;
        };
        let parsed = (
            fields[1].parse::<usize>(),
            fields[2].parse::<usize>(),
            fields[3].parse::<usize>(),
            fields[6].parse::<usize>(),
            fields[7].parse::<usize>(),
            fields[8].parse::<usize>(),
            fields[9].parse::<usize>(),
            fields[10].parse::<usize>(),
            fields[11].parse::<u8>(),
        );
        let (Ok(query_len), Ok(query_start), Ok(query_end), Ok(target_len), Ok(target_start), Ok(target_end), Ok(matches), Ok(aligned_length), Ok(mapq)) = parsed else {
            stats.malformed_or_unknown += 1;
            continue;
        };
        if query_len == 0
            || target_len == 0
            || query_start > query_end
            || query_end > query_len
            || target_start > target_end
            || target_end > target_len
            || aligned_length == 0
            || guide_sequence.len() != target_len
        {
            stats.malformed_or_unknown += 1;
            continue;
        }
        let query_fraction = (query_end - query_start) as f32 / query_len as f32;
        let identity = matches as f32 / aligned_length as f32;
        if aligned_length < cli.min_alignment
            || query_fraction < cli.min_query_fraction
            || identity < cli.min_identity
        {
            stats.threshold_filtered += 1;
            continue;
        }
        let reverse = fields[4] == "-";
        let primary = !fields[12..].contains(&"tp:A:S");
        let (full_start, full_end) = project_full_span(
            query_len,
            query_start,
            query_end,
            reverse,
            target_len,
            target_start,
            target_end,
        );
        if full_end <= full_start {
            stats.malformed_or_unknown += 1;
            continue;
        }
        by_segment[segment_id].push(Anchor {
            segment_id,
            guide_name: fields[5].to_string(),
            reverse,
            full_start,
            full_end,
            identity,
            query_fraction,
            mapq,
            aligned_length,
            rank_score: identity as f64 * aligned_length as f64,
            primary,
        });
    }

    let mut by_guide: FxHashMap<String, Vec<Anchor>> = FxHashMap::default();
    for mut hits in by_segment {
        if hits.is_empty() {
            continue;
        }
        hits.sort_unstable_by(|left, right| {
            right
                .primary
                .cmp(&left.primary)
                .then_with(|| right.rank_score.total_cmp(&left.rank_score))
                .then_with(|| right.mapq.cmp(&left.mapq))
        });
        let best = hits[0].clone();
        if !best.primary || best.mapq < cli.min_mapq {
            stats.low_mapq_or_non_primary += 1;
            continue;
        }
        let ambiguous = hits.iter().skip(1).any(|hit| {
            !same_locus(&best, hit)
                && hit.rank_score >= best.rank_score * f64::from(cli.secondary_ratio)
        });
        if ambiguous {
            stats.ambiguous += 1;
            continue;
        }
        stats.accepted += 1;
        by_guide
            .entry(best.guide_name.clone())
            .or_default()
            .push(best);
    }
    Ok((by_guide, stats))
}

fn project_full_span(
    query_len: usize,
    query_start: usize,
    query_end: usize,
    reverse: bool,
    target_len: usize,
    target_start: usize,
    target_end: usize,
) -> (i64, i64) {
    let query_left = query_start as i64;
    let query_right = (query_len - query_end) as i64;
    let (start, end) = if reverse {
        (
            target_start as i64 - query_right,
            target_end as i64 + query_left,
        )
    } else {
        (
            target_start as i64 - query_left,
            target_end as i64 + query_right,
        )
    };
    (start.clamp(0, target_len as i64), end.clamp(0, target_len as i64))
}

fn same_locus(left: &Anchor, right: &Anchor) -> bool {
    if left.guide_name != right.guide_name || left.reverse != right.reverse {
        return false;
    }
    let intersection = left.full_end.min(right.full_end) - left.full_start.max(right.full_start);
    if intersection <= 0 {
        return false;
    }
    let smaller = (left.full_end - left.full_start).min(right.full_end - right.full_start);
    intersection * 10 >= smaller * 8
}

fn read_model_scores(
    path: &Path,
    name_to_id: &FxHashMap<String, usize>,
) -> Result<FxHashMap<EdgeKey, ModelEvidence>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );
    let mut evidence: FxHashMap<EdgeKey, ModelEvidence> = FxHashMap::default();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if line_index == 0 && fields.first().is_some_and(|field| *field == "source") {
            continue;
        }
        if fields.len() < 3 {
            bail!("model score line {} has fewer than three columns", line_index + 1);
        }
        let source = parse_oriented_label(fields[0], name_to_id).with_context(|| {
            format!("unknown model-score source {} on line {}", fields[0], line_index + 1)
        })?;
        let target = parse_oriented_label(fields[1], name_to_id).with_context(|| {
            format!("unknown model-score target {} on line {}", fields[1], line_index + 1)
        })?;
        let score: f64 = fields[2]
            .parse()
            .with_context(|| format!("invalid model score on line {}", line_index + 1))?;
        if !score.is_finite() {
            bail!("model score on line {} is not finite", line_index + 1);
        }
        let decision = fields.get(3).copied().unwrap_or("neutral").to_ascii_lowercase();
        let scorer = fields.get(4).copied().unwrap_or("external").to_string();
        let entry = evidence
            .entry(EdgeKey { source, target })
            .or_default();
        entry.score_sum += score;
        entry.votes += 1;
        entry.veto |= matches!(decision.as_str(), "veto" | "reject");
        entry.scorers.insert(scorer);
    }
    Ok(evidence)
}

fn parse_oriented_label(
    label: &str,
    name_to_id: &FxHashMap<String, usize>,
) -> Option<OrientedNode> {
    let (name, reverse) = if let Some(name) = label.strip_suffix('+') {
        (name, false)
    } else if let Some(name) = label.strip_suffix('-') {
        (name, true)
    } else {
        (label, false)
    };
    Some(OrientedNode {
        id: *name_to_id.get(name)?,
        reverse,
    })
}

fn build_candidates(
    graph: &GfaGraph,
    guides: &FxHashMap<String, Vec<u8>>,
    anchors_by_guide: &FxHashMap<String, Vec<Anchor>>,
    model_evidence: &FxHashMap<EdgeKey, ModelEvidence>,
    config: SelectionConfig,
) -> (Vec<RankedCandidate>, FxHashMap<&'static str, usize>) {
    let mut aggregates: FxHashMap<EdgeKey, CandidateAggregate> = FxHashMap::default();
    let mut rejections = FxHashMap::default();

    for (guide_name, anchors) in anchors_by_guide {
        let Some(guide_sequence) = guides.get(guide_name) else {
            continue;
        };
        let chain = monotonic_chain(anchors.clone());
        for pair in chain.windows(2) {
            match make_candidate(
                &pair[0],
                &pair[1],
                guide_sequence,
                graph,
                config,
            ) {
                Ok(candidate) => {
                    let key = candidate.key;
                    if let Some(aggregate) = aggregates.get_mut(&key) {
                        aggregate.guides.insert(candidate.guide_name.clone());
                        aggregate.min_gap = aggregate.min_gap.min(candidate.projected_gap);
                        aggregate.max_gap = aggregate.max_gap.max(candidate.projected_gap);
                        let gap_delta = aggregate.max_gap.saturating_sub(aggregate.min_gap);
                        if gap_delta > config.max_gap_disagreement as i64 {
                            aggregate.conflicting_gap = true;
                        }
                        if candidate.base_score > aggregate.best.base_score {
                            aggregate.best = candidate;
                        }
                    } else {
                        let mut support = FxHashSet::default();
                        support.insert(candidate.guide_name.clone());
                        aggregates.insert(
                            key,
                            CandidateAggregate {
                                min_gap: candidate.projected_gap,
                                max_gap: candidate.projected_gap,
                                best: candidate,
                                guides: support,
                                conflicting_gap: false,
                            },
                        );
                    }
                }
                Err(reason) => increment(&mut rejections, reason),
            }
        }
    }

    let mut ranked = Vec::with_capacity(aggregates.len());
    for aggregate in aggregates.into_values() {
        let guide_support = aggregate.guides.len();
        let model = model_evidence.get(&aggregate.best.key);
        let model_score = model.map(|value| value.score_sum / value.votes.max(1) as f64);
        let model_votes = model.map_or(0, |value| value.votes);
        let model_scorers = model.map_or(0, |value| value.scorers.len());
        let (eligible, reason) = if aggregate.conflicting_gap {
            (false, "guide_gap_conflict")
        } else if guide_support < config.min_guide_support {
            (false, "insufficient_guide_support")
        } else if model.is_some_and(|value| value.veto) {
            (false, "model_veto")
        } else if config
            .min_model_score
            .zip(model_score)
            .is_some_and(|(minimum, score)| score < minimum)
        {
            (false, "low_model_score")
        } else {
            (true, "ok")
        };
        if !eligible {
            increment(&mut rejections, reason);
        }
        let score = aggregate.best.base_score
            + guide_support as f64 * 120.0
            + model_score.unwrap_or(0.0) * config.model_weight;
        ranked.push(RankedCandidate {
            best: aggregate.best,
            guide_support,
            model_score,
            model_votes,
            model_scorers,
            score,
            eligible,
            reason: reason.to_string(),
            selected: false,
        });
    }
    ranked.sort_unstable_by(|left, right| {
        right
            .score
            .total_cmp(&left.score)
            .then_with(|| right.guide_support.cmp(&left.guide_support))
            .then_with(|| right.best.direct_graph.cmp(&left.best.direct_graph))
            .then_with(|| left.best.key.source.cmp(&right.best.key.source))
            .then_with(|| left.best.key.target.cmp(&right.best.key.target))
    });
    (ranked, rejections)
}

fn increment(counts: &mut FxHashMap<&'static str, usize>, reason: &'static str) {
    *counts.entry(reason).or_insert(0) += 1;
}

fn monotonic_chain(mut anchors: Vec<Anchor>) -> Vec<Anchor> {
    anchors.sort_unstable_by(|left, right| {
        left.full_start
            .cmp(&right.full_start)
            .then_with(|| right.full_end.cmp(&left.full_end))
            .then_with(|| right.rank_score.total_cmp(&left.rank_score))
    });
    let mut chain: Vec<Anchor> = Vec::with_capacity(anchors.len());
    for anchor in anchors {
        if chain
            .last()
            .is_some_and(|previous| anchor.full_end <= previous.full_end)
        {
            continue;
        }
        chain.push(anchor);
    }
    chain
}

fn make_candidate(
    left: &Anchor,
    right: &Anchor,
    guide_sequence: &[u8],
    graph: &GfaGraph,
    config: SelectionConfig,
) -> std::result::Result<Candidate, &'static str> {
    if left.segment_id == right.segment_id || right.full_end <= left.full_end {
        return Err("non_progressing_anchor");
    }
    let source = OrientedNode {
        id: left.segment_id,
        reverse: left.reverse,
    };
    let target = OrientedNode {
        id: right.segment_id,
        reverse: right.reverse,
    };
    let projected_gap = right.full_start - left.full_end;
    if projected_gap > config.max_gap as i64 {
        return Err("guide_gap_too_large");
    }
    if projected_gap < -(config.max_overlap as i64) {
        return Err("projected_overlap_too_large");
    }

    let left_coverage = graph.segments[left.segment_id].coverage.max(0.001);
    let right_coverage = graph.segments[right.segment_id].coverage.max(0.001);
    let coverage_ratio = left_coverage.min(right_coverage) / left_coverage.max(right_coverage);
    if coverage_ratio < config.min_coverage_ratio {
        return Err("coverage_incompatible");
    }

    let left_sequence = oriented_sequence(&graph.segments[left.segment_id], left.reverse);
    let right_sequence = oriented_sequence(&graph.segments[right.segment_id], right.reverse);
    let direct_overlap = graph.links.get(&(source, target)).copied();
    let mut overlap = 0_usize;
    let mut bridge = Vec::new();

    if let Some(graph_overlap) = direct_overlap {
        if graph_overlap > left_sequence.len() || graph_overlap > right_sequence.len() {
            return Err("invalid_graph_overlap");
        }
        overlap = graph_overlap;
    } else if projected_gap < 0 {
        let expected = (-projected_gap) as usize;
        let Some((validated_overlap, _identity)) = find_suffix_prefix_overlap(
            &left_sequence,
            &right_sequence,
            expected,
            config,
        ) else {
            return Err("unvalidated_overlap");
        };
        overlap = validated_overlap;
    } else if projected_gap > 0 {
        if config.mode == Mode::OverlapOnly {
            return Err("positive_gap_in_overlap_mode");
        }
        if config.mode == Mode::Conservative {
            let strict_gap = config.conservative_max_gap.min(config.max_gap);
            if projected_gap > strict_gap as i64
                || left.identity.min(right.identity) < 0.985
                || left.query_fraction.min(right.query_fraction) < 0.90
                || left.mapq.min(right.mapq) < 40
            {
                return Err("not_conservative_enough");
            }
        }
        let start = left.full_end as usize;
        let end = right.full_start as usize;
        if start > end || end > guide_sequence.len() {
            return Err("guide_gap_out_of_bounds");
        }
        bridge.extend_from_slice(&guide_sequence[start..end]);
        let ambiguous = bridge
            .iter()
            .filter(|base| !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        let n_fraction = ambiguous as f32 / bridge.len().max(1) as f32;
        if n_fraction > config.max_bridge_n_fraction {
            return Err("guide_gap_too_ambiguous");
        }
    }

    let identity = left.identity.min(right.identity);
    let query_fraction = left.query_fraction.min(right.query_fraction);
    let mapq = left.mapq.min(right.mapq);
    let aligned_length = left.aligned_length.min(right.aligned_length);
    let direct_graph = direct_overlap.is_some();
    let base_score = f64::from(identity) * 1_000.0
        + f64::from(query_fraction) * 350.0
        + f64::from(mapq) * 3.0
        + f64::from(coverage_ratio) * 250.0
        + (aligned_length as f64).ln_1p() * 25.0
        + if direct_graph { 500.0 } else { 0.0 }
        + (overlap as f64).ln_1p() * 35.0
        - bridge.len() as f64 * 0.35;

    Ok(Candidate {
        key: EdgeKey { source, target },
        guide_name: left.guide_name.clone(),
        projected_gap,
        overlap,
        bridge,
        direct_graph,
        identity,
        query_fraction,
        mapq,
        aligned_length,
        coverage_ratio,
        base_score,
    })
}

fn oriented_sequence(segment: &Segment, reverse: bool) -> Vec<u8> {
    if reverse {
        reverse_complement(&segment.sequence)
    } else {
        segment.sequence.clone()
    }
}

fn find_suffix_prefix_overlap(
    left: &[u8],
    right: &[u8],
    expected: usize,
    config: SelectionConfig,
) -> Option<(usize, f32)> {
    let maximum = config.max_overlap.min(left.len()).min(right.len());
    if maximum < config.min_overlap {
        return None;
    }
    let center = expected.clamp(config.min_overlap, maximum);
    let start = center
        .saturating_sub(config.overlap_slack)
        .max(config.min_overlap);
    let end = center.saturating_add(config.overlap_slack).min(maximum);
    let mut best: Option<(usize, f32, usize)> = None;
    for length in start..=end {
        let matches = left[left.len() - length..]
            .iter()
            .zip(&right[..length])
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        let identity = matches as f32 / length as f32;
        if identity < config.min_overlap_identity {
            continue;
        }
        let distance = length.abs_diff(expected);
        let replace = best.is_none_or(|(best_length, best_identity, best_distance)| {
            identity
                .total_cmp(&best_identity)
                .then_with(|| best_distance.cmp(&distance))
                .then_with(|| length.cmp(&best_length))
                == Ordering::Greater
        });
        if replace {
            best = Some((length, identity, distance));
        }
    }
    best.map(|(length, identity, _distance)| (length, identity))
}

fn select_path_cover(
    ranked: &mut [RankedCandidate],
    segment_count: usize,
) -> (
    Vec<Option<Candidate>>,
    Vec<Option<usize>>,
    Vec<Option<bool>>,
) {
    let mut successor: Vec<Option<Candidate>> = vec![None; segment_count];
    let mut predecessor = vec![None; segment_count];
    let mut orientation = vec![None; segment_count];
    let mut parent: Vec<usize> = (0..segment_count).collect();

    for candidate in ranked {
        if !candidate.eligible {
            continue;
        }
        let source = candidate.best.key.source;
        let target = candidate.best.key.target;
        if successor[source.id].is_some() || predecessor[target.id].is_some() {
            continue;
        }
        if orientation[source.id].is_some_and(|value| value != source.reverse)
            || orientation[target.id].is_some_and(|value| value != target.reverse)
        {
            continue;
        }
        let source_root = find(&mut parent, source.id);
        let target_root = find(&mut parent, target.id);
        if source_root == target_root {
            continue;
        }
        successor[source.id] = Some(candidate.best.clone());
        predecessor[target.id] = Some(source.id);
        orientation[source.id] = Some(source.reverse);
        orientation[target.id] = Some(target.reverse);
        union(&mut parent, source_root, target_root);
        candidate.selected = true;
    }
    (successor, predecessor, orientation)
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

fn build_output_records(
    segments: &[Segment],
    successor: &[Option<Candidate>],
    predecessor: &[Option<usize>],
    orientation: &[Option<bool>],
    min_length: usize,
) -> Vec<OutputRecord> {
    let mut used = vec![false; segments.len()];
    let mut records = Vec::new();
    let mut seen = FxHashSet::default();

    for start in 0..segments.len() {
        if predecessor[start].is_some() || used[start] {
            continue;
        }
        let record = extend_record(start, segments, successor, orientation, &mut used);
        push_record(record, min_length, &mut seen, &mut records);
    }
    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        let record = extend_record(start, segments, successor, orientation, &mut used);
        push_record(record, min_length, &mut seen, &mut records);
    }
    records.sort_unstable_by(|left, right| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    records
}

fn extend_record(
    start: usize,
    segments: &[Segment],
    successor: &[Option<Candidate>],
    orientation: &[Option<bool>],
    used: &mut [bool],
) -> OutputRecord {
    let start_reverse = orientation[start].unwrap_or(false);
    let mut sequence = oriented_sequence(&segments[start], start_reverse);
    let mut current = start;
    let mut unitigs = 0_usize;
    let mut guide_edges = 0_usize;
    let mut guide_bases = 0_usize;

    for _ in 0..segments.len() {
        if used[current] {
            break;
        }
        used[current] = true;
        unitigs += 1;
        let Some(edge) = successor[current].as_ref() else {
            break;
        };
        let target = edge.key.target;
        let target_sequence = oriented_sequence(&segments[target.id], target.reverse);
        if edge.overlap > 0 {
            let overlap = edge.overlap.min(target_sequence.len());
            sequence.extend_from_slice(&target_sequence[overlap..]);
        } else {
            sequence.extend_from_slice(&edge.bridge);
            sequence.extend_from_slice(&target_sequence);
            guide_bases += edge.bridge.len();
        }
        guide_edges += 1;
        current = target.id;
    }

    OutputRecord {
        sequence,
        unitigs,
        guide_edges,
        guide_bases,
    }
}

fn push_record(
    mut record: OutputRecord,
    min_length: usize,
    seen: &mut FxHashSet<Vec<u8>>,
    records: &mut Vec<OutputRecord>,
) {
    if record.sequence.len() < min_length {
        return;
    }
    let reverse = reverse_complement(&record.sequence);
    if reverse < record.sequence {
        record.sequence = reverse;
    }
    if seen.insert(record.sequence.clone()) {
        records.push(record);
    }
}

fn write_fasta(path: &Path, records: &[OutputRecord]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    for (index, record) in records.iter().enumerate() {
        writeln!(
            writer,
            ">proteinguide_{:06} len={} unitigs={} guide_edges={} guide_bases={}",
            index + 1,
            record.sequence.len(),
            record.unitigs,
            record.guide_edges,
            record.guide_bases
        )?;
        for chunk in record.sequence.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_report(path: &Path, ranked: &[RankedCandidate], segments: &[Segment]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(
        writer,
        "source\ttarget\tselected\teligible\treason\tscore\tguide\tguide_support\tprojected_gap\toverlap\tguide_bases\tdirect_graph\tidentity\tquery_fraction\tmapq\taligned_length\tcoverage_ratio\tmodel_score\tmodel_votes\tmodel_scorers"
    )?;
    for candidate in ranked {
        let source = oriented_label(candidate.best.key.source, segments);
        let target = oriented_label(candidate.best.key.target, segments);
        let model_score = candidate
            .model_score
            .map_or_else(String::new, |score| format!("{score:.6}"));
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{}\t{}\t{}",
            source,
            target,
            candidate.selected,
            candidate.eligible,
            candidate.reason,
            candidate.score,
            candidate.best.guide_name,
            candidate.guide_support,
            candidate.best.projected_gap,
            candidate.best.overlap,
            candidate.best.bridge.len(),
            candidate.best.direct_graph,
            candidate.best.identity,
            candidate.best.query_fraction,
            candidate.best.mapq,
            candidate.best.aligned_length,
            candidate.best.coverage_ratio,
            model_score,
            candidate.model_votes,
            candidate.model_scorers
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn oriented_label(node: OrientedNode, segments: &[Segment]) -> String {
    format!(
        "{}{}",
        segments[node.id].name,
        if node.reverse { '-' } else { '+' }
    )
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

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> SelectionConfig {
        SelectionConfig {
            mode: Mode::Conservative,
            min_coverage_ratio: 0.1,
            max_gap: 300,
            conservative_max_gap: 90,
            max_overlap: 100,
            min_overlap: 3,
            min_overlap_identity: 1.0,
            overlap_slack: 2,
            max_gap_disagreement: 30,
            max_bridge_n_fraction: 0.05,
            min_guide_support: 1,
            model_weight: 400.0,
            min_model_score: None,
        }
    }

    #[test]
    fn projects_query_tails_on_both_strands() {
        assert_eq!(
            project_full_span(100, 10, 90, false, 1_000, 200, 280),
            (190, 290)
        );
        assert_eq!(
            project_full_span(100, 5, 80, true, 1_000, 200, 280),
            (180, 285)
        );
    }

    #[test]
    fn validates_expected_suffix_prefix_overlap() {
        let overlap = find_suffix_prefix_overlap(b"AAAACCCC", b"CCCCGGGG", 4, config());
        assert_eq!(overlap, Some((4, 1.0)));
    }

    #[test]
    fn parses_oriented_segment_labels() {
        let mut names = FxHashMap::default();
        names.insert("u7".to_string(), 7);
        assert_eq!(
            parse_oriented_label("u7-", &names),
            Some(OrientedNode {
                id: 7,
                reverse: true
            })
        );
    }
}
