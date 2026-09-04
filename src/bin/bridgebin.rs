use bridgeasm::bridgebin::{
    bin_contigs, read_coverage_table, read_fasta, write_outputs, BridgeBinConfig,
};
use bridgeasm::bridgebin_quant::{quantify_bins, write_abundance_table};
use bridgeasm::bridgebin_reconcile::{
    read_marker_table, reconcile_bins, MarkerTable, ReconcileConfig,
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
    out_dir: PathBuf,
    config: BridgeBinConfig,
    reconcile: ReconcileConfig,
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

    let initial = bin_contigs(&contigs, coverage.as_ref(), &cli.config);
    let (result, stats) = reconcile_bins(
        &contigs,
        coverage.as_ref(),
        markers.as_ref(),
        initial,
        &cli.config,
        &cli.reconcile,
    );
    write_outputs(&contigs, &result, &cli.out_dir, cli.emit_unbinned)?;

    if cli.reconcile.enabled {
        eprintln!(
            "bridgebin: reconciliation {} -> {} bins ({} merges, {} rescued contigs, {} marker-blocked comparisons)",
            stats.initial_bins,
            stats.final_bins,
            stats.merges,
            stats.rescued_contigs,
            stats.marker_blocked_pairs
        );
    }

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
        "bridgebin: {} contigs, {} bins, {} binned contigs, {:.2}% bp assigned",
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
        println!("bridgebin 0.3.0");
        process::exit(0);
    }

    let mut contigs: Option<PathBuf> = None;
    let mut coverage: Option<PathBuf> = None;
    let mut markers: Option<PathBuf> = None;
    let mut out_dir: Option<PathBuf> = None;
    let mut config = BridgeBinConfig::default();
    let mut reconcile = ReconcileConfig::default();
    let mut emit_unbinned = true;

    let mut i = 0usize;
    while i < args.len() {
        let key = &args[i];
        match key.as_str() {
            "--contigs" | "-c" => contigs = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--coverage" => coverage = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--markers" => markers = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--out-dir" | "-o" => out_dir = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--min-contig" => config.min_contig_len = parse_usize(value(&args, &mut i, key)?, key)?,
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
    if config.seed_min_len < config.min_contig_len {
        return Err("--seed-min-contig must be >= --min-contig".to_string());
    }
    if config.composition_weight + config.coverage_weight + config.gc_weight <= f64::EPSILON {
        return Err("at least one feature weight must be > 0".to_string());
    }

    Ok(Cli {
        contigs,
        coverage,
        markers,
        out_dir,
        config,
        reconcile,
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
        "bridgebin 0.3.0 - evidence-aware metagenomic binning and quantification\n\n\
USAGE:\n  bridgebin --contigs <FASTA> --out-dir <DIR> [OPTIONS]\n\n\
INPUT:\n  -c, --contigs <FASTA>          Assembled contigs\n      --coverage <TSV>           Coverage matrix: contig sample1 [sample2 ...]\n      --markers <TSV>            Optional single-copy markers: contig marker (one hit/row)\n  -o, --out-dir <DIR>            Output directory\n\n\
INITIAL BINNING:\n      --min-contig <BP>          Minimum contig length [default: 1500]\n      --seed-min-contig <BP>     Minimum length for seed contigs [default: 2500]\n      --join-threshold <0..1>    Seed-to-bin similarity threshold [default: 0.76]\n      --rescue-threshold <0..1>  First-pass rescue threshold [default: 0.70]\n      --rescue-margin <0..1>     Required best-vs-second margin [default: 0.025]\n      --composition-weight <N>   Composition weight [default: 0.45]\n      --coverage-weight <N>      Multi-sample coverage weight [default: 0.50]\n      --gc-weight <N>            GC-content weight [default: 0.05]\n\n\
BIN RECONCILIATION:\n      --no-reconcile                     Disable v1 bin-bin reconciliation\n      --reconcile-threshold <0..1>       Reciprocal bin merge threshold [default: 0.72]\n      --reconcile-margin <0..1>          Best-vs-second merge margin [default: 0.015]\n      --reconcile-min-composition <0..1> Hard 5-mer composition floor [default: 0.58]\n      --reconcile-min-coverage <0..1>    Hard coverage-consistency floor [default: 0.62]\n      --same-coverage-min-composition <0..1> Composition floor when coverage is nearly identical [default: 0.82]\n      --reconcile-max-gc-delta <0..1>    Hard GC delta [default: 0.055]\n      --post-rescue-threshold <0..1>     Rescue after bin merges [default: 0.72]\n      --post-rescue-margin <0..1>        Post-merge rescue margin [default: 0.03]\n      --reconcile-max-merges <N>         Safety cap [default: 256]\n\n\
QUANTIFICATION:\n  With --coverage, BridgeBin writes abundance.tsv containing length-weighted\n  median depth, length-weighted mean depth, and relative abundance per bin/sample.\n\n\
OUTPUT:\n  assignments.tsv, bins.tsv, bins/bin_XXXX.fa, abundance.tsv (with coverage),\n  and optionally unbinned.fa.\n\n\
BridgeBin v1 uses conservative seed bins followed by reciprocal bin-bin\nreconciliation. Optional single-copy marker conflicts are hard negative evidence."
    );
}
