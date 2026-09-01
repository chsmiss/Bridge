use anyhow::{bail, Result};
use bridgeasm::dna::base_bits;
use bridgeasm::fastq::for_each_pair;
use clap::{Parser, ValueEnum};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Repeated,
    OnePass,
}

#[derive(Debug, Parser)]
struct Args {
    #[arg(short = '1', long)]
    read1: PathBuf,
    #[arg(short = '2', long)]
    read2: Option<PathBuf>,
    #[arg(long, default_value = "21,31,41,55")]
    k: String,
    #[arg(long, value_enum)]
    mode: Mode,
    #[arg(long)]
    max_pairs: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
struct Roller {
    k: usize,
    mask: u128,
    forward: u128,
    reverse: u128,
    valid: usize,
}

impl Roller {
    fn new(k: usize) -> Result<Self> {
        if k == 0 || k > 63 {
            bail!("benchmark roller supports k in 1..=63; got {k}");
        }
        let bits = 2 * k;
        let mask = if bits == 128 {
            u128::MAX
        } else {
            (1_u128 << bits) - 1
        };
        Ok(Self {
            k,
            mask,
            forward: 0,
            reverse: 0,
            valid: 0,
        })
    }

    #[inline]
    fn reset(&mut self) {
        self.forward = 0;
        self.reverse = 0;
        self.valid = 0;
    }

    #[inline]
    fn push(&mut self, base: u8) -> Option<u128> {
        let Some(bits) = base_bits(base) else {
            self.reset();
            return None;
        };
        let complement = 3_u8 - bits;
        self.forward = ((self.forward << 2) | u128::from(bits)) & self.mask;
        self.reverse = (self.reverse >> 2) | (u128::from(complement) << (2 * (self.k - 1)));
        self.valid = self.valid.saturating_add(1).min(self.k);
        (self.valid >= self.k).then_some(self.forward.min(self.reverse))
    }
}

#[derive(Debug, Default)]
struct Counter {
    roller: Option<Roller>,
    observations: u64,
    counts: FxHashMap<u128, u32>,
}

impl Counter {
    fn new(k: usize) -> Result<Self> {
        Ok(Self {
            roller: Some(Roller::new(k)?),
            observations: 0,
            counts: FxHashMap::default(),
        })
    }

    fn reset_roller(&mut self) {
        if let Some(roller) = &mut self.roller {
            roller.reset();
        }
    }

    fn add_sequence(&mut self, sequence: &[u8]) {
        let roller = self.roller.as_mut().expect("roller exists");
        roller.reset();
        for &base in sequence {
            if let Some(key) = roller.push(base) {
                self.observations += 1;
                let entry = self.counts.entry(key).or_insert(0);
                *entry = entry.saturating_add(1);
            }
        }
    }
}

#[derive(Debug, Serialize)]
struct KSummary {
    k: usize,
    observations: u64,
    distinct: usize,
}

#[derive(Debug, Serialize)]
struct Summary {
    mode: String,
    read_pairs: usize,
    seconds: f64,
    k: Vec<KSummary>,
}

fn parse_ks(value: &str) -> Result<Vec<usize>> {
    let mut ks = value
        .split(',')
        .filter(|part| !part.trim().is_empty())
        .map(|part| part.trim().parse::<usize>().map_err(Into::into))
        .collect::<Result<Vec<_>>>()?;
    ks.sort_unstable();
    ks.dedup();
    if ks.is_empty() {
        bail!("at least one k is required");
    }
    for &k in &ks {
        if k == 0 || k > 63 {
            bail!("benchmark supports k in 1..=63; got {k}");
        }
    }
    Ok(ks)
}

fn repeated(args: &Args, ks: &[usize]) -> Result<Summary> {
    let started = Instant::now();
    let mut summaries = Vec::new();
    let mut read_pairs = 0;
    for &k in ks {
        let mut counter = Counter::new(k)?;
        read_pairs = for_each_pair(
            &args.read1,
            args.read2.as_deref(),
            args.max_pairs,
            |_index, left, right| {
                counter.add_sequence(&left.sequence);
                if let Some(right) = right {
                    counter.add_sequence(&right.sequence);
                }
                counter.reset_roller();
                Ok(())
            },
        )?;
        summaries.push(KSummary {
            k,
            observations: counter.observations,
            distinct: counter.counts.len(),
        });
    }
    Ok(Summary {
        mode: "repeated".to_string(),
        read_pairs,
        seconds: started.elapsed().as_secs_f64(),
        k: summaries,
    })
}

fn one_pass(args: &Args, ks: &[usize]) -> Result<Summary> {
    let started = Instant::now();
    let mut counters = ks
        .iter()
        .map(|&k| Counter::new(k))
        .collect::<Result<Vec<_>>>()?;
    let read_pairs = for_each_pair(
        &args.read1,
        args.read2.as_deref(),
        args.max_pairs,
        |_index, left, right| {
            for counter in &mut counters {
                counter.add_sequence(&left.sequence);
                if let Some(right) = right.as_ref() {
                    counter.add_sequence(&right.sequence);
                }
                counter.reset_roller();
            }
            Ok(())
        },
    )?;
    Ok(Summary {
        mode: "one-pass".to_string(),
        read_pairs,
        seconds: started.elapsed().as_secs_f64(),
        k: ks
            .iter()
            .copied()
            .zip(counters)
            .map(|(k, counter)| KSummary {
                k,
                observations: counter.observations,
                distinct: counter.counts.len(),
            })
            .collect(),
    })
}

fn main() -> Result<()> {
    let args = Args::parse();
    let ks = parse_ks(&args.k)?;
    let summary = match args.mode {
        Mode::Repeated => repeated(&args, &ks)?,
        Mode::OnePass => one_pass(&args, &ks)?,
    };
    println!("{}", serde_json::to_string_pretty(&summary)?);
    Ok(())
}
