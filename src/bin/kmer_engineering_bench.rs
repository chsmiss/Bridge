use anyhow::{bail, Result};
use bridgeasm::dna::{base_bits, complement_bits, KmerKey};
use bridgeasm::fastq::for_each_pair;
use clap::{Parser, ValueEnum};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    RepeatedWide,
    SinglePassCompact,
}

#[derive(Parser, Debug)]
struct Args {
    #[arg(short = '1', long)]
    read1: PathBuf,
    #[arg(short = '2', long)]
    read2: Option<PathBuf>,
    #[arg(long, default_value = "21,25,31,41,55")]
    k: String,
    #[arg(long, value_enum)]
    mode: Mode,
    #[arg(long)]
    max_pairs: Option<usize>,
}

#[derive(Debug, Serialize)]
struct KSummary {
    k: usize,
    observations: u64,
    distinct: usize,
    key_bytes: usize,
}

#[derive(Debug, Serialize)]
struct Summary {
    mode: String,
    read_pairs: usize,
    elapsed_seconds: f64,
    k: Vec<KSummary>,
}

#[derive(Debug)]
struct WideCounter {
    k: usize,
    forward: KmerKey,
    reverse: KmerKey,
    valid: usize,
    observations: u64,
    counts: FxHashMap<KmerKey, u32>,
}

impl WideCounter {
    fn new(k: usize) -> Self {
        Self {
            k,
            forward: KmerKey::ZERO,
            reverse: KmerKey::ZERO,
            valid: 0,
            observations: 0,
            counts: FxHashMap::default(),
        }
    }

    fn reset(&mut self) {
        self.forward = KmerKey::ZERO;
        self.reverse = KmerKey::ZERO;
        self.valid = 0;
    }

    fn push(&mut self, bits: u8) {
        self.forward.shift_left_append(bits, self.k);
        self.reverse.shift_right_prepend_complement(bits, self.k);
        self.valid = self.valid.saturating_add(1).min(self.k);
        if self.valid < self.k {
            return;
        }
        let key = self.forward.min(self.reverse);
        let entry = self.counts.entry(key).or_insert(0);
        *entry = entry.saturating_add(1);
        self.observations += 1;
    }

    fn record(&mut self, sequence: &[u8]) {
        self.reset();
        for &base in sequence {
            match base_bits(base) {
                Some(bits) => self.push(bits),
                None => self.reset(),
            }
        }
    }
}

#[derive(Debug)]
enum CompactCounter {
    U64 {
        k: usize,
        forward: u64,
        reverse: u64,
        valid: usize,
        mask: u64,
        observations: u64,
        counts: FxHashMap<u64, u32>,
    },
    U128 {
        k: usize,
        forward: u128,
        reverse: u128,
        valid: usize,
        mask: u128,
        observations: u64,
        counts: FxHashMap<u128, u32>,
    },
}

impl CompactCounter {
    fn new(k: usize) -> Result<Self> {
        if k == 0 || k > 63 {
            bail!("compact benchmark supports k in 1..=63, got {k}");
        }
        if k <= 31 {
            let bits = 2 * k;
            let mask = if bits == 64 {
                u64::MAX
            } else {
                (1_u64 << bits) - 1
            };
            Ok(Self::U64 {
                k,
                forward: 0,
                reverse: 0,
                valid: 0,
                mask,
                observations: 0,
                counts: FxHashMap::default(),
            })
        } else {
            let bits = 2 * k;
            let mask = if bits == 128 {
                u128::MAX
            } else {
                (1_u128 << bits) - 1
            };
            Ok(Self::U128 {
                k,
                forward: 0,
                reverse: 0,
                valid: 0,
                mask,
                observations: 0,
                counts: FxHashMap::default(),
            })
        }
    }

    fn reset(&mut self) {
        match self {
            Self::U64 {
                forward,
                reverse,
                valid,
                ..
            } => {
                *forward = 0;
                *reverse = 0;
                *valid = 0;
            }
            Self::U128 {
                forward,
                reverse,
                valid,
                ..
            } => {
                *forward = 0;
                *reverse = 0;
                *valid = 0;
            }
        }
    }

    fn push(&mut self, bits: u8) {
        match self {
            Self::U64 {
                k,
                forward,
                reverse,
                valid,
                mask,
                observations,
                counts,
            } => {
                *forward = ((*forward << 2) | u64::from(bits)) & *mask;
                *reverse = (*reverse >> 2) | (u64::from(complement_bits(bits)) << (2 * (*k - 1)));
                *valid = valid.saturating_add(1).min(*k);
                if *valid >= *k {
                    let key = (*forward).min(*reverse);
                    let entry = counts.entry(key).or_insert(0);
                    *entry = entry.saturating_add(1);
                    *observations += 1;
                }
            }
            Self::U128 {
                k,
                forward,
                reverse,
                valid,
                mask,
                observations,
                counts,
            } => {
                *forward = ((*forward << 2) | u128::from(bits)) & *mask;
                *reverse = (*reverse >> 2) | (u128::from(complement_bits(bits)) << (2 * (*k - 1)));
                *valid = valid.saturating_add(1).min(*k);
                if *valid >= *k {
                    let key = (*forward).min(*reverse);
                    let entry = counts.entry(key).or_insert(0);
                    *entry = entry.saturating_add(1);
                    *observations += 1;
                }
            }
        }
    }

    fn k(&self) -> usize {
        match self {
            Self::U64 { k, .. } | Self::U128 { k, .. } => *k,
        }
    }

    fn observations(&self) -> u64 {
        match self {
            Self::U64 { observations, .. } | Self::U128 { observations, .. } => *observations,
        }
    }

    fn distinct(&self) -> usize {
        match self {
            Self::U64 { counts, .. } => counts.len(),
            Self::U128 { counts, .. } => counts.len(),
        }
    }

    fn key_bytes(&self) -> usize {
        match self {
            Self::U64 { .. } => std::mem::size_of::<u64>(),
            Self::U128 { .. } => std::mem::size_of::<u128>(),
        }
    }
}

fn parse_ks(text: &str) -> Result<Vec<usize>> {
    let mut ks = Vec::new();
    for raw in text.split(',') {
        let k: usize = raw.trim().parse()?;
        if k == 0 || k > 63 {
            bail!("benchmark k must be in 1..=63, got {k}");
        }
        if !ks.contains(&k) {
            ks.push(k);
        }
    }
    ks.sort_unstable();
    if ks.is_empty() {
        bail!("at least one k is required");
    }
    Ok(ks)
}

fn repeated_wide(args: &Args, ks: &[usize]) -> Result<Summary> {
    let started = Instant::now();
    let mut read_pairs = 0;
    let mut summaries = Vec::new();
    for &k in ks {
        let mut counter = WideCounter::new(k);
        read_pairs = for_each_pair(
            &args.read1,
            args.read2.as_deref(),
            args.max_pairs,
            |_index, left, right| {
                counter.record(&left.sequence);
                if let Some(right) = right {
                    counter.record(&right.sequence);
                }
                Ok(())
            },
        )?;
        summaries.push(KSummary {
            k,
            observations: counter.observations,
            distinct: counter.counts.len(),
            key_bytes: std::mem::size_of::<KmerKey>(),
        });
    }
    Ok(Summary {
        mode: "repeated-wide".to_string(),
        read_pairs,
        elapsed_seconds: started.elapsed().as_secs_f64(),
        k: summaries,
    })
}

fn single_pass_compact(args: &Args, ks: &[usize]) -> Result<Summary> {
    let started = Instant::now();
    let mut counters: Vec<CompactCounter> = ks
        .iter()
        .copied()
        .map(CompactCounter::new)
        .collect::<Result<_>>()?;
    let read_pairs = for_each_pair(
        &args.read1,
        args.read2.as_deref(),
        args.max_pairs,
        |_index, left, right| {
            for sequence in std::iter::once(left.sequence.as_slice())
                .chain(right.as_ref().map(|record| record.sequence.as_slice()))
            {
                for counter in &mut counters {
                    counter.reset();
                }
                for &base in sequence {
                    let Some(bits) = base_bits(base) else {
                        for counter in &mut counters {
                            counter.reset();
                        }
                        continue;
                    };
                    for counter in &mut counters {
                        counter.push(bits);
                    }
                }
            }
            Ok(())
        },
    )?;
    let summaries = counters
        .iter()
        .map(|counter| KSummary {
            k: counter.k(),
            observations: counter.observations(),
            distinct: counter.distinct(),
            key_bytes: counter.key_bytes(),
        })
        .collect();
    Ok(Summary {
        mode: "single-pass-compact".to_string(),
        read_pairs,
        elapsed_seconds: started.elapsed().as_secs_f64(),
        k: summaries,
    })
}

fn main() -> Result<()> {
    let args = Args::parse();
    let ks = parse_ks(&args.k)?;
    let summary = match args.mode {
        Mode::RepeatedWide => repeated_wide(&args, &ks)?,
        Mode::SinglePassCompact => single_pass_compact(&args, &ks)?,
    };
    println!("{}", serde_json::to_string_pretty(&summary)?);
    Ok(())
}
