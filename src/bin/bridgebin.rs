use bridgeasm::bridgebin::{
    bin_contigs, read_coverage_table, read_fasta, write_outputs, BridgeBinConfig,
};
use std::env;
use std::io;
use std::path::PathBuf;
use std::process;

#[derive(Debug)]
struct Cli {
    contigs: PathBuf,
    coverage: Option<PathBuf>,
    out_dir: PathBuf,
    config: BridgeBinConfig,
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

    let result = bin_contigs(&contigs, coverage.as_ref(), &cli.config);
    write_outputs(&contigs, &result, &cli.out_dir, cli.emit_unbinned)?;

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
        println!("bridgebin 0.1.0");
        process::exit(0);
    }

    let mut contigs: Option<PathBuf> = None;
    let mut coverage: Option<PathBuf> = None;
    let mut out_dir: Option<PathBuf> = None;
    let mut config = BridgeBinConfig::default();
    let mut emit_unbinned = true;

    let mut i = 0usize;
    while i < args.len() {
        let key = &args[i];
        match key.as_str() {
            "--contigs" | "-c" => contigs = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--coverage" => coverage = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--out-dir" | "-o" => out_dir = Some(PathBuf::from(value(&args, &mut i, key)?)),
            "--min-contig" => {
                config.min_contig_len = parse_usize(value(&args, &mut i, key)?, key)?
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
            "--gc-weight" => {
                config.gc_weight = parse_nonnegative(value(&args, &mut i, key)?, key)?
            }
            "--no-unbinned" => emit_unbinned = false,
            unknown => {
                return Err(format!(
                    "unknown argument '{unknown}'\n\nRun bridgebin --help for usage."
                ))
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
        out_dir,
        config,
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
        "bridgebin 0.1.0 - evidence-aware metagenomic contig binning\n\n\
USAGE:\n  bridgebin --contigs <FASTA> --out-dir <DIR> [OPTIONS]\n\n\
INPUT:\n  -c, --contigs <FASTA>          Assembled contigs\n      --coverage <TSV>           Optional coverage matrix: contig sample1 [sample2 ...]\n  -o, --out-dir <DIR>            Output directory\n\n\
BINNING:\n      --min-contig <BP>          Minimum contig length [default: 1500]\n      --seed-min-contig <BP>     Minimum length for seed contigs [default: 2500]\n      --join-threshold <0..1>    Seed-to-bin similarity threshold [default: 0.76]\n      --rescue-threshold <0..1>  Short-contig rescue threshold [default: 0.70]\n      --rescue-margin <0..1>     Required best-vs-second margin [default: 0.025]\n      --composition-weight <N>   Canonical TNF weight [default: 0.45]\n      --coverage-weight <N>      Multi-sample coverage weight [default: 0.50]\n      --gc-weight <N>            GC-content weight [default: 0.05]\n      --no-unbinned              Do not write unbinned.fa\n\n\
OUTPUT:\n  assignments.tsv, bins.tsv, bins/bin_XXXX.fa, and optionally unbinned.fa\n\n\
The v0 algorithm uses exact centroid search with composition + differential coverage.\n\
It is intentionally simple and deterministic; graph/read/protein evidence are planned next."
    );
}
