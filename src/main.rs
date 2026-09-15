use anyhow::Result;
use bridgeasm::assembler::{assemble, AssembleConfig};
use bridgeasm::dna::MAX_K;
use bridgeasm::multik::{write_multik_outputs, MultiKConfig};
use bridgeasm::multik_fast::build_multik_graph_fast;
use bridgeasm::output::write_outputs;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "bridgeasm",
    version,
    about = "Evidence-aware short-read metagenome assembler"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Assemble paired or single-end FASTQ reads with one fixed k.
    Assemble {
        #[arg(short = '1', long)]
        read1: PathBuf,
        #[arg(short = '2', long)]
        read2: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short = 'k', long, default_value_t = 31)]
        k: usize,
        #[arg(long, default_value_t = 2)]
        min_count: u32,
        #[arg(long, default_value_t = 16)]
        mercy_max_kmers: usize,
        #[arg(long, default_value_t = 1)]
        mercy_min_support: u16,
        #[arg(long, default_value_t = 25.0)]
        mercy_min_quality: f32,
        #[arg(long, default_value_t = 2)]
        min_read_support: u32,
        #[arg(long, default_value_t = 2)]
        min_pair_support: u32,
        #[arg(long, default_value_t = 5)]
        min_primary_support: u32,
        #[arg(long, default_value_t = 0.75)]
        primary_dominance: f32,
        /// Use dominant same-read unitig triplets to extend the conservative path cover.
        #[arg(long, default_value_t = false)]
        threaded_path_cover: bool,
        /// Build a major-strain path cover by greedily matching high-confidence graph edges.
        #[arg(long, default_value_t = false)]
        major_path_cover: bool,
        /// Minimum support fraction at the weaker endpoint of a major-path edge.
        #[arg(long, default_value_t = 0.20)]
        path_cover_secondary_dominance: f32,
        #[arg(long, default_value_t = 200)]
        min_contig_length: usize,
        #[arg(long, default_value_t = 100)]
        scaffold_gap_bases: usize,
        #[arg(long)]
        max_pairs: Option<usize>,
        #[arg(short = 't', long, default_value_t = 1)]
        threads: usize,
    },
    /// Stage34: one FASTQ pass, parallel k layers, exact cross-k projections.
    Multik {
        #[arg(short = '1', long)]
        read1: PathBuf,
        #[arg(short = '2', long)]
        read2: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
        /// Comma-separated k values, for example 21,31,41,55.
        #[arg(long, value_delimiter = ',', default_value = "21,31,41,55")]
        ks: Vec<usize>,
        #[arg(long, default_value_t = 2)]
        min_count: u32,
        /// Count support once per original physical fragment for each k-mer.
        #[arg(long, default_value_t = 2)]
        min_fragment_support: u32,
        #[arg(long, default_value_t = 20.0)]
        min_mean_quality: f32,
        #[arg(long)]
        max_pairs: Option<usize>,
        /// Maximum lower-k unique walk when probing a high-k dead-end rescue.
        #[arg(long, default_value_t = 500)]
        max_rescue_bases: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Assemble {
            read1,
            read2,
            output,
            k,
            min_count,
            mercy_max_kmers,
            mercy_min_support,
            mercy_min_quality,
            min_read_support,
            min_pair_support,
            min_primary_support,
            primary_dominance,
            threaded_path_cover,
            major_path_cover,
            path_cover_secondary_dominance,
            min_contig_length,
            scaffold_gap_bases,
            max_pairs,
            threads,
        } => {
            if k == 0 || k > MAX_K {
                anyhow::bail!("k must be in 1..={MAX_K}");
            }
            if !(0.0..=60.0).contains(&mercy_min_quality) {
                anyhow::bail!("mercy minimum quality must be in 0..=60");
            }
            if !(0.0..=1.0).contains(&path_cover_secondary_dominance)
                || path_cover_secondary_dominance > primary_dominance
            {
                anyhow::bail!("path-cover secondary dominance must be in 0..=primary-dominance");
            }
            let config = AssembleConfig {
                read1,
                read2,
                output_dir: output.clone(),
                k,
                min_count,
                mercy_max_kmers,
                mercy_min_support,
                mercy_min_quality,
                min_read_support,
                min_pair_support,
                min_primary_support,
                primary_dominance,
                threaded_path_cover,
                major_path_cover,
                path_cover_secondary_dominance,
                min_contig_length,
                scaffold_gap_bases,
                max_pairs,
                threads,
            };
            let product = assemble(&config)?;
            write_outputs(&product, &output)?;
            eprintln!(
                "assembled {} primary contigs (N50 {}, total {} bp), {} bubbles / {} haplotigs",
                product.stats.primary_contigs,
                product.stats.primary_n50,
                product.stats.primary_bases,
                product.stats.simple_bubbles,
                product.stats.haplotigs
            );
        }
        Command::Multik {
            read1,
            read2,
            output,
            ks,
            min_count,
            min_fragment_support,
            min_mean_quality,
            max_pairs,
            max_rescue_bases,
        } => {
            let config = MultiKConfig {
                read1,
                read2,
                output_dir: output.clone(),
                ks,
                min_count,
                min_fragment_support,
                min_mean_quality,
                max_pairs,
                max_rescue_bases,
            };
            let graph = build_multik_graph_fast(&config)?;
            write_multik_outputs(&graph, &output)?;
            eprintln!(
                "built {} multi-k layers from {} physical read pairs in {:.3}s; {} cross-k rescue candidates",
                graph.layers.len(),
                graph.summary.read_pairs,
                graph.summary.timings_seconds.total_seconds,
                graph.rescue_candidates.len()
            );
        }
    }
    Ok(())
}
