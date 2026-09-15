use anyhow::Result;
use bridgeasm::assembler::{assemble, AssembleConfig};
use bridgeasm::dna::MAX_K;
use bridgeasm::multik::MultiKConfig;
use bridgeasm::multik_hybrid::{run_multik_hybrid, write_hybrid_outputs, HybridConfig};
use bridgeasm::multik_hybrid_safe::enforce_hybrid_safe_gate;
use bridgeasm::multik_stream::write_streaming_outputs;
use bridgeasm::multik_v4::run_multik_v4;
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
        #[arg(long, default_value_value_t = 5)]
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
    /// Stage34 v4: four-pass compact global multi-k graph.
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
        /// Worker threads for Stage34 fixed-memory scans and cross-k analysis.
        #[arg(short = 't', long, default_value_t = 4)]
        threads: usize,
    },
    /// Stage34 v5: compact global backbone plus conservative local adaptive high-k refinement.
    MultikHybrid {
        #[arg(short = '1', long)]
        read1: PathBuf,
        #[arg(short = '2', long)]
        read2: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
        /// Global compact backbone k. Must fit the u128 backbone representation.
        #[arg(long, default_value_t = 31)]
        backbone_k: usize,
        /// First local high-k candidate. Larger k values are chosen adaptively from routed read lengths.
        #[arg(long, default_value_t = 55)]
        local_start_k: usize,
        /// Hard cap for local wide-key refinement. The v5 wide key supports k <= 255.
        #[arg(long, default_value_t = 255)]
        max_local_k: usize,
        #[arg(long, default_value_t = 2)]
        min_count: u32,
        #[arg(long, default_value_t = 2)]
        min_fragment_support: u32,
        #[arg(long, default_value_t = 20.0)]
        min_mean_quality: f32,
        /// Raw-graph state radius around branch/dead-end seeds.
        #[arg(long, default_value_t = 64)]
        neighborhood_radius: usize,
        /// Maximum raw oriented states in one local neighborhood.
        #[arg(long, default_value_t = 4096)]
        max_neighborhood_states: usize,
        /// Maximum adaptive neighborhoods. Limited to 64 so routing uses one u64 mask.
        #[arg(long, default_value_t = 64)]
        max_neighborhoods: usize,
        #[arg(long)]
        max_pairs: Option<usize>,
        #[arg(short = 't', long, default_value_t = 4)]
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
            threads,
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
            let run = run_multik_v4(&config, threads)?;
            write_streaming_outputs(&run, &output)?;
            eprintln!(
                "built {} four-pass packed-edge multi-k layers from {} physical read pairs in {:.3}s using {} threads; {} cross-k rescue candidates",
                run.summary.layers.len(),
                run.summary.read_pairs,
                run.summary.timings_seconds.total_seconds,
                threads.max(1),
                run.rescue_candidates.len()
            );
        }
        Command::MultikHybrid {
            read1,
            read2,
            output,
            backbone_k,
            local_start_k,
            max_local_k,
            min_count,
            min_fragment_support,
            min_mean_quality,
            neighborhood_radius,
            max_neighborhood_states,
            max_neighborhoods,
            max_pairs,
            threads,
        } => {
            let config = HybridConfig {
                read1,
                read2,
                output_dir: output.clone(),
                backbone_k,
                local_start_k,
                max_local_k,
                min_count,
                min_fragment_support,
                min_mean_quality,
                neighborhood_radius,
                max_neighborhood_states,
                max_neighborhoods,
                max_pairs,
            };
            let mut run = run_multik_hybrid(&config, threads)?;
            enforce_hybrid_safe_gate(&mut run);
            write_hybrid_outputs(&run, &output)?;
            let resolved = run
                .summary
                .local
                .iter()
                .filter(|item| item.resolved)
                .count();
            let above_55 = run
                .summary
                .local
                .iter()
                .filter_map(|item| item.selected_k)
                .filter(|&k| k > 55)
                .count();
            eprintln!(
                "hybrid Stage34 built k={} backbone from {} pairs and refined {} neighborhoods in {:.3}s using {} threads; {} connectivity-safe resolved, {} selected k>55",
                run.summary.backbone_k,
                run.summary.read_pairs,
                run.summary.neighborhoods,
                run.summary.total_seconds,
                threads.max(1),
                resolved,
                above_55
            );
        }
    }
    Ok(())
}
