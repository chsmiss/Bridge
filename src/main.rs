use anyhow::Result;
use bridgeasm::assembler::{assemble, AssembleConfig};
use bridgeasm::dna::MAX_K;
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
    /// Assemble paired or single-end FASTQ reads.
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
        #[arg(long, default_value_t = 200)]
        min_contig_length: usize,
        #[arg(long)]
        max_pairs: Option<usize>,
        #[arg(short = 't', long, default_value_t = 1)]
        threads: usize,
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
            min_contig_length,
            max_pairs,
            threads,
        } => {
            if k == 0 || k > MAX_K {
                anyhow::bail!("k must be in 1..={MAX_K}");
            }
            if !(0.0..=60.0).contains(&mercy_min_quality) {
                anyhow::bail!("mercy minimum quality must be in 0..=60");
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
                min_contig_length,
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
    }
    Ok(())
}
