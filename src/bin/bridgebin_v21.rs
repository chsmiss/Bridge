use bridgeasm::bridgebin::{read_coverage_table, read_fasta, write_outputs};
use bridgeasm::bridgebin_quant::{quantify_bins, write_abundance_table};
use bridgeasm::bridgebin_reconcile::{read_marker_table, MarkerTable};
use bridgeasm::bridgebin_v2::{
    bin_contigs_v2, read_bio_feature_table, read_link_table, BioFeatureTable, BridgeBinV2Config,
    LinkTable,
};
use bridgeasm::bridgebin_v21::{
    read_pair_score_table, refine_bins_v21, BridgeBinV21Config, PairScoreTable,
};
use bridgeasm::bridgebin_v21_merge::{
    merge_bins_v21, BridgeBinV21MergeConfig,
};
use std::env;
use std::io;
use std::path::PathBuf;
use std::process;

#[derive(Debug)]
struct Cli {
    contigs: PathBuf,
    coverage: Option<PathBuf>,
    markers: Option<PathBuf>,
    bio_features: Option<PathBuf>,
    links: Option<PathBuf>,
    pair_scores: PathBuf,
    out_dir: PathBuf,
    v2: BridgeBinV2Config,
    v21: BridgeBinV21Config,
    merge: BridgeBinV21MergeConfig,
    enable_merge: bool,
    emit_unbinned: bool,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("bridgebin_v21: {error}");
        process::exit(1);
    }
}

fn run() -> io::Result<()> {
    let cli = parse_args().map_err(|message| io::Error::new(io::ErrorKind::InvalidInput, message))?;
    let contigs = read_fasta(&cli.contigs)?;
    let coverage = match cli.coverage.as_ref() {
        Some(path) => Some(read_coverage_table(path)?),
        None => None,
    };
    let markers: Option<MarkerTable> = match cli.markers.as_ref() {
        Some(path) => Some(read_marker_table(path)?),
        None => None,
    };
    let bio: Option<BioFeatureTable> = match cli.bio_features.as_ref() {
        Some(path) => Some(read_bio_feature_table(path)?),
        None => None,
    };
    let links: Option<LinkTable> = match cli.links.as_ref() {
        Some(path) => Some(read_link_table(path)?),
        None => None,
    };
    let pair_scores: PairScoreTable = read_pair_score_table(&cli.pair_scores)?;

    let (base, base_stats) = bin_contigs_v2(
        &contigs,
        coverage.as_ref(),
        markers.as_ref(),
        bio.as_ref(),
        links.as_ref(),
        &cli.v2,
    );
    eprintln!(
        "bridgebin_v21: v2 eligible={} candidates={} core_bins={} rescued_components={} unbinned={}",
        base_stats.eligible_contigs,
        base_stats.candidate_edges,
        base_stats.core_bins,
        base_stats.rescued_components,
        base_stats.unbinned_contigs
    );

    let (refined, refine_stats) = refine_bins_v21(
        &contigs,
        markers.as_ref(),
        base,
        &pair_scores,
        &cli.v21,
    );
    eprintln!(
        "bridgebin_v21: refine {} -> {} bins; conflicted={} split={} hard_negative_pairs={} marker_negative_pairs={} positive_pairs={} rescued_contigs={} ambiguous_residuals={}",
        refine_stats.input_bins,
        refine_stats.output_bins,
        refine_stats.conflicted_input_bins,
        refine_stats.split_bins,
        refine_stats.hard_negative_pairs,
        refine_stats.marker_negative_pairs,
        refine_stats.positive_pairs,
        refine_stats.rescued_contigs,
        refine_stats.ambiguous_residuals
    );

    let result = if cli.enable_merge {
        let (merged, merge_stats) = merge_bins_v21(
            &contigs,
            markers.as_ref(),
            refined,
            &pair_scores,
            &cli.v21,
            &cli.merge,
        );
        eprintln!(
            "bridgebin_v21: merge {} -> {} bins; candidates={} accepted={} hard_blocks={} marker_blocks={}",
            merge_stats.input_bins,
            merge_stats.output_bins,
            merge_stats.candidate_bin_pairs,
            merge_stats.accepted_merges,
            merge_stats.hard_blocked_bin_pairs,
            merge_stats.marker_blocked_bin_pairs
        );
        merged
    } else {
        refined
    };

    write_outputs(&contigs, &result, &cli.out_dir, cli.emit_unbinned)?;
    if let Some(table) = coverage.as_ref() {
        let abundance = quantify_bins(&result, table);
        write_abundance_table(&abundance, &cli.out_dir)?;
    }

    let total_bp: usize = contigs.iter().map(|contig| contig.seq.len()).sum();
    let binned_bp: usize = result
        .assignments
        .iter()
        .filter(|assignment| assignment.bin_index.is_some())
        .map(|assignment| assignment.length)
        .sum();
    eprintln!(
        "bridgebin_v21: final_bins={} assigned_bp={:.2}% pair_scores={}",
        result.bins.len(),
        if total_bp == 0 {
            0.0
        } else {
            100.0 * binned_bp as f64 / total_bp as f64
        },
        pair_scores.values.len()
    );
    Ok(())
}

fn parse_args() -> Result<Cli, String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|value| value == "-h" || value == "--help") {
        print_help();
        process::exit(0);
    }

    let mut contigs = None;
    let mut coverage = None;
    let mut markers = None;
    let mut bio_features = None;
    let mut links = None;
    let mut pair_scores = None;
    let mut out_dir = None;
    let mut v2 = BridgeBinV2Config::default();
    // Start v2.1 from the empirically safer balanced signed-graph regime.
    v2.core_min_attraction = 0.75;
    v2.component_min_attraction = 0.67;
    v2.rescue_min_attraction = 0.70;
    v2.rescue_margin = 0.03;
    let mut v21 = BridgeBinV21Config::default();
    let mut merge = BridgeBinV21MergeConfig::default();
    let mut enable_merge = true;
    let mut emit_unbinned = true;

    let mut index = 0usize;
    while index < args.len() {
        let key = &args[index];
        match key.as_str() {
            "--contigs" | "-c" => contigs = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--coverage" => coverage = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--markers" => markers = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--bio-features" => {
                bio_features = Some(PathBuf::from(value(&args, &mut index, key)?))
            }
            "--links" => links = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--pair-scores" => pair_scores = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--out-dir" | "-o" => out_dir = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--min-contig" => v2.min_contig_len = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v2-max-neighbors" => v2.max_neighbors = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v2-min-component-bp" => v2.min_component_bp = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v2-core-attraction" => v2.core_min_attraction = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v2-component-attraction" => v2.component_min_attraction = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v2-rescue-attraction" => v2.rescue_min_attraction = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v2-rescue-margin" => v2.rescue_margin = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-min-confidence" => v21.min_pair_confidence = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-split-max-same" => v21.split_max_same = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-join-min-same" => v21.join_min_same = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-rescue-min-same" => v21.rescue_min_same = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-rescue-margin" => v21.rescue_margin = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-min-pair-support" => v21.min_pair_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v21-top-pair-support" => v21.top_pair_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v21-min-subbin-bp" => v21.min_subbin_bp = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v21-merge-min-same" => merge.min_same = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-merge-min-support" => merge.min_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v21-merge-top-support" => merge.top_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--no-v21-merge" => enable_merge = false,
            "--no-unbinned" => emit_unbinned = false,
            unknown => return Err(format!("unknown argument '{unknown}'\n\nRun bridgebin_v21 --help for usage.")),
        }
        index += 1;
    }

    if v2.max_neighbors == 0 || v21.min_pair_support == 0 || v21.top_pair_support == 0 {
        return Err("neighbor/support counts must be positive".to_string());
    }
    if v21.split_max_same >= v21.join_min_same {
        return Err("--v21-split-max-same must be lower than --v21-join-min-same".to_string());
    }

    Ok(Cli {
        contigs: contigs.ok_or_else(|| "missing required --contigs <FASTA>".to_string())?,
        coverage,
        markers,
        bio_features,
        links,
        pair_scores: pair_scores.ok_or_else(|| "missing required --pair-scores <TSV>".to_string())?,
        out_dir: out_dir.ok_or_else(|| "missing required --out-dir <DIR>".to_string())?,
        v2,
        v21,
        merge,
        enable_merge,
        emit_unbinned,
    })
}

fn value(args: &[String], index: &mut usize, key: &str) -> Result<String, String> {
    *index += 1;
    args.get(*index)
        .cloned()
        .ok_or_else(|| format!("missing value for {key}"))
}

fn parse_usize(raw: String, key: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid integer for {key}: {raw}"))
}

fn parse_unit(raw: String, key: &str) -> Result<f64, String> {
    let value = raw
        .parse::<f64>()
        .map_err(|_| format!("invalid number for {key}: {raw}"))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(format!("{key} must be between 0 and 1"));
    }
    Ok(value)
}

fn print_help() {
    println!(
        "bridgebin_v21 - biological pair refinement on top of the signed BridgeBin v2 graph\n\n\
USAGE:\n  bridgebin_v21 --contigs <FASTA> --pair-scores <TSV> --out-dir <DIR> [OPTIONS]\n\n\
INPUTS:\n      --coverage <TSV>           Optional multi-sample coverage matrix\n      --markers <TSV>            Optional single-copy-marker table\n      --bio-features <TSV>       Optional DNA/gene/ESM-C/taxonomy contig features for v2\n      --links <TSV>              Optional graph/read/Hi-C signed evidence\n      --pair-scores <TSV>        Required pair head output: left right p_same confidence model\n\n\
V2.1 PAIR REFINEMENT:\n      --v21-min-confidence <P>       Pair confidence floor [0.80]\n      --v21-split-max-same <P>       Hard negative p_same ceiling [0.12]\n      --v21-join-min-same <P>        Positive within-bin edge floor [0.88]\n      --v21-rescue-min-same <P>      Residual rescue floor [0.84]\n      --v21-rescue-margin <P>        Best-vs-second rescue margin [0.08]\n      --v21-min-pair-support <N>     Required links from residual to bin [2]\n      --v21-min-subbin-bp <BP>       Keep split component if at least this large [20000]\n      --no-v21-merge                 Disable high-confidence cross-bin biological merge\n\n\
The pair-score file is intended to come from scripts/bridgebin_pair_head.py, whose\nfeatures can include DNABERT-S/GENERanno DNA embeddings, GENERanno CDS architecture,\nESM-C protein embeddings/repertoire, coverage/composition, taxonomy and physical links."
    );
}
