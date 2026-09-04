use bridgeasm::bridgebin::{
    bin_contigs, read_coverage_table, read_fasta, write_outputs, BridgeBinConfig,
};
use bridgeasm::bridgebin_quant::{quantify_bins, write_abundance_table};
use bridgeasm::bridgebin_reconcile::{
    read_marker_table, reconcile_bins, MarkerTable, ReconcileConfig,
};
use bridgeasm::bridgebin_v2::{
    bin_contigs_v2, read_bio_feature_table, read_link_table, BioFeatureTable, BridgeBinV2Config,
    LinkTable,
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
    out_dir: PathBuf,
    algorithm: String,
    config: BridgeBinConfig,
    reconcile: ReconcileConfig,
    v2: BridgeBinV2Config,
    emit_unbinned: bool,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("bridgebin: {err}");
        process::exit(1);
    }
}

fn run() -> io::Result<()> {
    let cli = parse_args().map_err(|msg| io::Error::new(io::ErrorKind::InvalidInput, msg))?;
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

    let result = match cli.algorithm.as_str() {
        "v0" => bin_contigs(&contigs, coverage.as_ref(), &cli.config),
        "v1" => {
            let initial = bin_contigs(&contigs, coverage.as_ref(), &cli.config);
            let (result, stats) = reconcile_bins(
                &contigs,
                coverage.as_ref(),
                markers.as_ref(),
                initial,
                &cli.config,
                &cli.reconcile,
            );
            eprintln!(
                "bridgebin: v1 reconciliation {} -> {} bins ({} merges, {} rescued contigs, {} marker-blocked comparisons)",
                stats.initial_bins,
                stats.final_bins,
                stats.merges,
                stats.rescued_contigs,
                stats.marker_blocked_pairs
            );
            result
        }
        "v2" => {
            let (result, stats) = bin_contigs_v2(
                &contigs,
                coverage.as_ref(),
                markers.as_ref(),
                bio.as_ref(),
                links.as_ref(),
                &cli.v2,
            );
            eprintln!(
                "bridgebin: v2 signed graph eligible={} candidates={} accepted_edges={} core_bins={} rescued_components={} marker_blocks={} taxonomy_blocks={} external_blocks={} unbinned={}",
                stats.eligible_contigs,
                stats.candidate_edges,
                stats.accepted_core_edges,
                stats.core_bins,
                stats.rescued_components,
                stats.marker_blocked_edges,
                stats.taxonomy_blocked_edges,
                stats.external_blocked_edges,
                stats.unbinned_contigs,
            );
            result
        }
        _ => unreachable!("algorithm validated by parse_args"),
    };

    write_outputs(&contigs, &result, &cli.out_dir, cli.emit_unbinned)?;

    if let Some(table) = coverage.as_ref() {
        let abundance = quantify_bins(&result, table);
        write_abundance_table(&abundance, &cli.out_dir)?;
        eprintln!(
            "bridgebin: quantified {} bins across {} samples",
            result.bins.len(),
            table.sample_names.len()
        );
    }

    let binned = result
        .assignments
        .iter()
        .filter(|a| a.bin_index.is_some())
        .count();
    let binned_bp: usize = result
        .assignments
        .iter()
        .filter(|a| a.bin_index.is_some())
        .map(|a| a.length)
        .sum();
    let total_bp: usize = contigs.iter().map(|c| c.seq.len()).sum();
    eprintln!(
        "bridgebin: algorithm={} {} contigs, {} bins, {} binned contigs, {:.2}% bp assigned",
        cli.algorithm,
        contigs.len(),
        result.bins.len(),
        binned,
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
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print_help();
        process::exit(0);
    }
    if args.iter().any(|a| a == "--version") {
        println!("bridgebin 0.4.0-dev");
        process::exit(0);
    }

    let mut contigs: Option<PathBuf> = None;
    let mut coverage: Option<PathBuf> = None;
    let mut markers: Option<PathBuf> = None;
    let mut bio_features: Option<PathBuf> = None;
    let mut links: Option<PathBuf> = None;
    let mut out_dir: Option<PathBuf> = None;
    let mut algorithm = "v2".to_string();
    let mut config = BridgeBinConfig::default();
    let mut reconcile = ReconcileConfig::default();
    let mut v2 = BridgeBinV2Config::default();
    let mut emit_unbinned = true;

    let mut i = 0usize;
    while i < args.len() {
        let key = &args[i];
        match key.as_str() {
            "--contigs" | "-c" => contigs = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--coverage" => coverage = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--markers" => markers = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--bio-features" => bio_features = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--links" => links = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--out-dir" | "-o" => out_dir = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--algorithm" => algorithm = value(&args, &mut i, key)?.to_ascii_lowercase(),
            "--min-contig" => {
                let parsed = parse_usize(value(&args, &mut i, key)?, key)?;
                config.min_contig_len = parsed;
                v2.min_contig_len = parsed;
            }
            "--seed-min-contig" => {
                config.seed_min_len = parse_usize(value(&args, &mut i, key)?, key)?
            }
            "--join-threshold" => {
                config.join_threshold = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--rescue-threshold" => {
                config.rescue_threshold = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--rescue-margin" => {
                config.rescue_margin = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--composition-weight" => {
                config.composition_weight = parse_nonnegative(value(&args, &mut i, key)?, key)?
            }
            "--coverage-weight" => {
                config.coverage_weight = parse_nonnegative(value(&args, &mut i, key)?, key)?
            }
            "--gc-weight" => config.gc_weight = parse_nonnegative(value(&args, &mut i, key)?, key)?,
            "--no-reconcile" => reconcile.enabled = false,
            "--reconcile-threshold" => {
                reconcile.merge_threshold = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--reconcile-margin" => {
                reconcile.merge_margin = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--reconcile-min-composition" => {
                reconcile.min_composition_score = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--reconcile-min-coverage" => {
                reconcile.min_coverage_score = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--same-coverage-min-composition" => {
                reconcile.same_coverage_min_composition = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--reconcile-max-gc-delta" => {
                reconcile.max_gc_delta = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--post-rescue-threshold" => {
                reconcile.post_rescue_threshold = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--post-rescue-margin" => {
                reconcile.post_rescue_margin = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--reconcile-max-merges" => {
                reconcile.max_merges = parse_usize(value(&args, &mut i, key)?, key)?
            }
            "--v2-max-neighbors" => {
                v2.max_neighbors = parse_usize(value(&args, &mut i, key)?, key)?
            }
            "--v2-min-component-bp" => {
                v2.min_component_bp = parse_usize(value(&args, &mut i, key)?, key)?
            }
            "--v2-core-attraction" => {
                v2.core_min_attraction = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-max-repulsion" => {
                v2.core_max_repulsion = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-component-attraction" => {
                v2.component_min_attraction = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-rescue-attraction" => {
                v2.rescue_min_attraction = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-rescue-margin" => {
                v2.rescue_margin = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-max-gc-delta" => {
                v2.max_gc_delta = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-min-composition" => {
                v2.min_component_composition = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-min-coverage" => {
                v2.min_component_coverage = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-taxonomy-confidence" => {
                v2.taxonomy_confidence = parse_unit(value(&args, &mut i, key)?, key)?
            }
            "--v2-soft-markers" => v2.hard_marker_veto = false,
            "--no-unbinned" => emit_unbinned = false,
            unknown => {
                return Err(format!(
                    "unknown argument '{unknown}'\n\nRun bridgebin --help for usage."
                ));
            }
        }
        i += 1;
    }

    let contigs = contigs.ok_or_else(|| "missing required --contigs <FASTA>".to_string())?;
    let out_dir = out_dir.ok_or_else(|| "missing required --out-dir <DIR>".to_string())?;
    if !matches!(algorithm.as_str(), "v0" | "v1" | "v2") {
        return Err("--algorithm must be one of v0, v1, v2".to_string());
    }
    if config.seed_min_len < config.min_contig_len {
        return Err("--seed-min-contig must be >= --min-contig".to_string());
    }
    if config.composition_weight + config.coverage_weight + config.gc_weight <= f64::EPSILON {
        return Err("at least one v0/v1 feature weight must be > 0".to_string());
    }
    if v2.max_neighbors == 0 {
        return Err("--v2-max-neighbors must be > 0".to_string());
    }

    Ok(Cli {
        contigs,
        coverage,
        markers,
        bio_features,
        links,
        out_dir,
        algorithm,
        config,
        reconcile,
        v2,
        emit_unbinned,
    })
}

fn value(args: &[String], i: &mut usize, key: &str) -> Result<String, String> {
    *i += 1;
    args.get(*i)
        .cloned()
        .ok_or_else(|| format!("missing value for {key}"))
}

fn parse_usize(raw: String, key: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid integer for {key}: {raw}"))
}

fn parse_unit(raw: String, key: &str) -> Result<f64, String> {
    let v = raw
        .parse::<f64>()
        .map_err(|_| format!("invalid number for {key}: {raw}"))?;
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        return Err(format!("{key} must be between 0 and 1"));
    }
    Ok(v)
}

fn parse_nonnegative(raw: String, key: &str) -> Result<f64, String> {
    let v = raw
        .parse::<f64>()
        .map_err(|_| format!("invalid number for {key}: {raw}"))?;
    if !v.is_finite() || v < 0.0 {
        return Err(format!("{key} must be finite and >= 0"));
    }
    Ok(v)
}

fn print_help() {
    println!(
        "bridgebin 0.4.0-dev - signed multimodal metagenomic binning\n\n\
USAGE:\n  bridgebin --contigs <FASTA> --out-dir <DIR> [OPTIONS]\n\n\
CORE:\n      --algorithm <v0|v1|v2>     Binning engine [default: v2]\n  -c, --contigs <FASTA>          Assembled contigs\n      --coverage <TSV>           Coverage matrix: contig sample1 [sample2 ...]\n      --markers <TSV>            Single-copy markers: contig marker (one hit/row)\n      --bio-features <TSV>       Optional taxonomy/gene/ESM-C feature table\n      --links <TSV>              Optional graph/read/Hi-C must/cannot-link evidence\n  -o, --out-dir <DIR>            Output directory\n      --min-contig <BP>          Minimum contig length [default: 1500]\n\n\
V2 SIGNED EVIDENCE GRAPH:\n      --v2-max-neighbors <N>         Sparse candidate neighbors [default: 64]\n      --v2-min-component-bp <BP>     Minimum core-bin size [default: 20000]\n      --v2-core-attraction <0..1>    Edge attraction floor [default: 0.80]\n      --v2-max-repulsion <0..1>      Allowed repulsion ceiling [default: 0.20]\n      --v2-component-attraction <N>  Component centroid attraction floor [default: 0.72]\n      --v2-rescue-attraction <N>     Residual-component rescue floor [default: 0.76]\n      --v2-rescue-margin <N>         Best-vs-second rescue margin [default: 0.05]\n      --v2-max-gc-delta <N>          Component GC delta ceiling [default: 0.075]\n      --v2-min-composition <N>       Component 5-mer similarity floor [default: 0.50]\n      --v2-min-coverage <N>          Component coverage similarity floor [default: 0.48]\n      --v2-taxonomy-confidence <N>   Confidence for taxonomy cannot-link [default: 0.90]\n      --v2-soft-markers              Do not make duplicate SCGs a hard veto\n\n\
BIO FEATURE TSV:\n  Header-based TSV. Supported columns: contig, taxonomy, taxonomy_confidence,\n  gene_profile, gene_confidence, esm_embedding/protein_embedding, protein_confidence.\n  Vector columns are comma-separated floats. Missing modalities are allowed.\n\n\
LINK TSV:\n  Header-based TSV with left/source, right/target, optional must_link, cannot_link,\n  and evidence_source. A cannot_link >= 0.95 is component-transitive hard negative\n  evidence. Positive graph/read links increase confidence but never force a merge.\n\n\
LEGACY V0/V1 OPTIONS:\n      --seed-min-contig, --join-threshold, --rescue-threshold, --rescue-margin,\n      --composition-weight, --coverage-weight, --gc-weight, and v1 reconcile flags\n      remain available for controlled benchmark comparisons.\n\n\
OUTPUT:\n  assignments.tsv, bins.tsv, bins/bin_XXXX.fa, abundance.tsv (with coverage),\n  and optionally unbinned.fa."
    );
}
