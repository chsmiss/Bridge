use anyhow::Result;
use bridgeasm::dna::{base_bits, canonical_kmers, KmerKey};
use bridgeasm::fastq::for_each_pair;
use clap::{Parser, ValueEnum};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::mem::size_of;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Generic,
    Compact,
}

#[derive(Debug, Parser)]
struct Args {
    #[arg(short = '1', long)]
    read1: PathBuf,
    #[arg(short = '2', long)]
    read2: Option<PathBuf>,
    #[arg(short = 'k', long, default_value_t = 21)]
    k: usize,
    #[arg(long, value_enum)]
    mode: Mode,
    #[arg(long)]
    max_pairs: Option<usize>,
}

#[derive(Debug, Serialize)]
struct Summary {
    mode: String,
    k: usize,
    read_pairs: usize,
    observations: u64,
    distinct: usize,
    seconds: f64,
    key_bytes: usize,
}

#[inline]
fn add_compact(sequence: &[u8], k: usize, counts: &mut FxHashMap<u128, u32>) -> u64 {
    if sequence.len() < k || k == 0 || k > 63 {
        return 0;
    }
    let mask = (1_u128 << (2 * k)) - 1;
    let mut forward = 0_u128;
    let mut reverse = 0_u128;
    let mut valid = 0_usize;
    let mut observations = 0_u64;
    for &base in sequence {
        let Some(bits) = base_bits(base) else {
            forward = 0;
            reverse = 0;
            valid = 0;
            continue;
        };
        let complement = 3_u8 - bits;
        forward = ((forward << 2) | u128::from(bits)) & mask;
        reverse = (reverse >> 2) | (u128::from(complement) << (2 * (k - 1)));
        valid = valid.saturating_add(1).min(k);
        if valid >= k {
            let key = forward.min(reverse);
            let entry = counts.entry(key).or_insert(0);
            *entry = entry.saturating_add(1);
            observations += 1;
        }
    }
    observations
}

fn main() -> Result<()> {
    let args = Args::parse();
    let started = Instant::now();
    let mut observations = 0_u64;
    match args.mode {
        Mode::Generic => {
            let mut counts: FxHashMap<KmerKey, u32> = FxHashMap::default();
            let read_pairs = for_each_pair(
                &args.read1,
                args.read2.as_deref(),
                args.max_pairs,
                |_index, left, right| {
                    for item in canonical_kmers(&left.sequence, args.k)? {
                        let entry = counts.entry(item.key).or_insert(0);
                        *entry = entry.saturating_add(1);
                        observations += 1;
                    }
                    if let Some(right) = right {
                        for item in canonical_kmers(&right.sequence, args.k)? {
                            let entry = counts.entry(item.key).or_insert(0);
                            *entry = entry.saturating_add(1);
                            observations += 1;
                        }
                    }
                    Ok(())
                },
            )?;
            println!(
                "{}",
                serde_json::to_string_pretty(&Summary {
                    mode: "generic".into(),
                    k: args.k,
                    read_pairs,
                    observations,
                    distinct: counts.len(),
                    seconds: started.elapsed().as_secs_f64(),
                    key_bytes: size_of::<KmerKey>(),
                })?
            );
        }
        Mode::Compact => {
            let mut counts: FxHashMap<u128, u32> = FxHashMap::default();
            let read_pairs = for_each_pair(
                &args.read1,
                args.read2.as_deref(),
                args.max_pairs,
                |_index, left, right| {
                    observations += add_compact(&left.sequence, args.k, &mut counts);
                    if let Some(right) = right {
                        observations += add_compact(&right.sequence, args.k, &mut counts);
                    }
                    Ok(())
                },
            )?;
            println!(
                "{}",
                serde_json::to_string_pretty(&Summary {
                    mode: "compact".into(),
                    k: args.k,
                    read_pairs,
                    observations,
                    distinct: counts.len(),
                    seconds: started.elapsed().as_secs_f64(),
                    key_bytes: size_of::<u128>(),
                })?
            );
        }
    }
    Ok(())
}
