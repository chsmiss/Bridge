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
    v21_merge: BridgeBinV21MergeConfig,
    emit_unbinned: bool,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("bridgebin-v21: {error}");
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

    // Stage 1: conservative signed-graph core from cheap and physical evidence.
    let (initial, v2_stats) = bin_contigs_v2(
        &contigs,
        coverage.as_ref(),
        markers.as_ref(),
        bio.as_ref(),
        links.as_ref(),
        &cli.v2,
    );

    // Stage 2: split contaminated bins and conservatively rescue residual contigs using
    // calibrated biological pair evidence. Hard negatives and duplicated SCGs survive
    // all transitive joins inside this stage.
    let (refined, v21_stats) = refine_bins_v21(
        &contigs,
        markers.as_ref(),
        initial,
        &pair_scores,
        &cli.v21,
    );

    // Stage 3: reconnect over-split pure components. A bin merge requires several
    // independent high-p_same links; any confident low-p_same pair or shared SCG blocks
    // the merge, including through later transitive unions.
    let (result, merge_stats) = merge_bins_v21(
        &contigs,
        markers.as_ref(),
        refined,
        &pair_scores,
        &cli.v21,
        &cli.v21_merge,
    );

    write_outputs(&contigs, &result, &cli.out_dir, cli.emit_unbinned)?;
    if let Some(table) = coverage.as_ref() {
        let abundance = quantify_bins(&result, table);
        write_abundance_table(&abundance, &cli.out_dir)?;
    }

    let binned_bp: usize = result
        .assignments
        .iter()
        .filter(|assignment| assignment.bin_index.is_some())
        .map(|assignment| assignment.length)
        .sum();
    let total_bp: usize = contigs.iter().map(|contig| contig.seq.len()).sum();
    eprintln!(
        "bridgebin-v21: v2 core_bins={} rescued_components={} unbinned={} | split input_bins={} split_output={} conflicted_bins={} split_bins={} hard_negatives={} marker_negatives={} positive_pairs={} rescued_contigs={} ambiguous={} | merge input_bins={} output_bins={} candidates={} accepted={} hard_blocked={} marker_blocked={} | assigned_bp={:.2}%",
        v2_stats.core_bins,
        v2_stats.rescued_components,
        v2_stats.unbinned_contigs,
        v21_stats.input_bins,
        v21_stats.output_bins,
        v21_stats.conflicted_input_bins,
        v21_stats.split_bins,
        v21_stats.hard_negative_pairs,
        v21_stats.marker_negative_pairs,
        v21_stats.positive_pairs,
        v21_stats.rescued_contigs,
        v21_stats.ambiguous_residuals,
        merge_stats.input_bins,
        merge_stats.output_bins,
        merge_stats.candidate_bin_pairs,
        merge_stats.accepted_merges,
        merge_stats.hard_blocked_bin_pairs,
        merge_stats.marker_blocked_bin_pairs,
        if total_bp == 0 {
            0.0
        } else {
            100.0 * binned_bp as f64 / total_bp as f64
        }
    );
    Ok(())
}

fn parse_args() -> Result<Cli, String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        process::exit(0);
    }
    if args.iter().any(|arg| arg == "--version") {
        println!("bridgebin-v21 0.4.1-dev");
        process::exit(0);
    }

    let mut contigs = None;
    let mut coverage = None;
    let mut markers = None;
    let mut bio_features = None;
    let mut links = None;
    let mut pair_scores = None;
    let mut out_dir = None;
    let mut emit_unbinned = true;

    // Start from the non-oracle balanced v2 preset established by the Zymo
    // sensitivity run. Pair-level biology is a refinement layer, not an excuse
    // to loosen the signed-graph core further.
    let mut v2 = BridgeBinV2Config::default();
    v2.core_min_attraction = 0.74;
    v2.component_min_attraction = 0.66;
    v2.rescue_min_attraction = 0.68;
    v2.rescue_margin = 0.02;
    let mut v21 = BridgeBinV21Config::default();
    let mut v21_merge = BridgeBinV21MergeConfig::default();

    let mut index = 0usize;
    while index < args.len() {
        let key = &args[index];
        match key.as_str() {
            "--contigs" | "-c" => contigs = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--coverage" => coverage = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--markers" => markers = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--bio-features" => bio_features = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--links" => links = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--pair-scores" => pair_scores = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--out-dir" | "-o" => out_dir = Some(PathBuf::from(value(&args, &mut index, key)?)),
            "--min-contig" => v2.min_contig_len = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v2-max-neighbors" => v2.max_neighbors = parse_usize(value(&args, &mut index, key)?, key)?,
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
            "--v21-bin-merge-min-same" => v21_merge.min_same = parse_unit(value(&args, &mut index, key)?, key)?,
            "--v21-bin-merge-min-support" => v21_merge.min_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--v21-bin-merge-top-support" => v21_merge.top_support = parse_usize(value(&args, &mut index, key)?, key)?,
            "--no-unbinned" => emit_unbinned = false,
            unknown => return Err(format!("unknown argument {unknown:?}")),
        }
        index += 1;
    }

    if v2.max_neighbors == 0
        || v21.min_pair_support == 0
        || v21.top_pair_support == 0
        || v21_merge.min_support == 0
        || v21_merge.top_support == 0
    {
        return Err("neighbor/support counts must be positive".to_string());
    }
    if v21.split_max_same >= v21.join_min_same {
        return Err("--v21-split-max-same must be lower than --v21-join-min-same".to_string());
    }
    if v21_merge.min_same <= v21.split_max_same {
        return Err("--v21-bin-merge-min-same must exceed --v21-split-max-same".to_string());
    }

    Ok(Cli {
        contigs: contigs.ok_or_else(|| "missing --contigs".to_string())?,
        coverage,
        markers,
        bio_features,
        links,
        pair_scores: pair_scores.ok_or_else(|| "missing --pair-scores".to_string())?,
        out_dir: out_dir.ok_or_else(|| "missing --out-dir".to_string())?,
        v2,
        v21,
        v21_merge,
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
        return Err(format!("{key} must be in [0,1]"));
    }
    Ok(value)
}

fn print_help() {
    println!(
        "bridgebin-v21 - BridgeBin v2 signed graph plus calibrated biological pair refinement\n\n\
USAGE:\n  bridgebin-v21 --contigs <FASTA> --pair-scores <TSV> --out-dir <DIR> [OPTIONS]\n\n\
INPUTS:\n  -c, --contigs <FASTA>          assembled contigs\n      --coverage <TSV>           multi-sample coverage matrix\n      --markers <TSV>            single-copy marker hits\n      --bio-features <TSV>       contig-level taxonomy/gene/protein features for v2\n      --links <TSV>              physical must/cannot-link evidence for v2\n      --pair-scores <TSV>        left right p_same [confidence] [model]\n  -o, --out-dir <DIR>            output directory\n\n\
V2.1 PAIR REFINEMENT:\n      --v21-min-confidence <P>       minimum calibrated pair confidence [0.80]\n      --v21-split-max-same <P>       p_same at/below this is a hard negative [0.12]\n      --v21-join-min-same <P>        p_same at/above this joins split seeds [0.88]\n      --v21-rescue-min-same <P>      residual-to-bin posterior floor [0.84]\n      --v21-rescue-margin <P>        best-vs-second rescue margin [0.08]\n      --v21-min-pair-support <N>     independent supports required for rescue [2]\n      --v21-top-pair-support <N>     strongest supports aggregated per bin [8]\n      --v21-min-subbin-bp <BP>       keep split component above this size [20000]\n\n\
V2.1 CROSS-BIN MERGE:\n      --v21-bin-merge-min-same <P>   mean high-confidence p_same required [0.92]\n      --v21-bin-merge-min-support <N> independent cross-bin supports required [3]\n      --v21-bin-merge-top-support <N> strongest cross-bin supports aggregated [12]\n\n\
The pair scorer is intentionally external. DNA foundation models, GENERanno-derived\nfeatures, ESM-C, or a learned multimodal head can all emit the same calibrated p_same TSV.\nHard negatives and duplicated single-copy markers remain vetoes through transitive merges."
    );
}
