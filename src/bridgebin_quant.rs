use crate::bridgebin::{BinningResult, CoverageTable};
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct BinAbundance {
    pub bin_index: usize,
    pub sample_index: usize,
    pub sample_name: String,
    pub robust_depth: f64,
    pub mean_depth: f64,
    pub relative_abundance: f64,
    pub covered_bp: usize,
    pub covered_contigs: usize,
}

pub fn quantify_bins(result: &BinningResult, coverage: &CoverageTable) -> Vec<BinAbundance> {
    let n_bins = result.bins.len();
    let n_samples = coverage.sample_names.len();
    let mut values = vec![vec![Vec::<(f64, usize)>::new(); n_samples]; n_bins];

    for assignment in &result.assignments {
        let Some(bin_index) = assignment.bin_index else {
            continue;
        };
        let Some(row) = coverage.values.get(&assignment.contig_id) else {
            continue;
        };
        if row.len() != n_samples {
            continue;
        }
        for (sample_index, depth) in row.iter().copied().enumerate() {
            values[bin_index][sample_index].push((depth, assignment.length));
        }
    }

    let mut robust = vec![vec![0.0f64; n_samples]; n_bins];
    let mut means = vec![vec![0.0f64; n_samples]; n_bins];
    let mut covered_bp = vec![vec![0usize; n_samples]; n_bins];
    let mut covered_contigs = vec![vec![0usize; n_samples]; n_bins];

    for bin_index in 0..n_bins {
        for sample_index in 0..n_samples {
            let observations = &mut values[bin_index][sample_index];
            robust[bin_index][sample_index] = weighted_median(observations);
            let bp: usize = observations.iter().map(|(_, len)| *len).sum();
            let weighted_sum: f64 = observations
                .iter()
                .map(|(depth, len)| *depth * *len as f64)
                .sum();
            means[bin_index][sample_index] = if bp == 0 {
                0.0
            } else {
                weighted_sum / bp as f64
            };
            covered_bp[bin_index][sample_index] = bp;
            covered_contigs[bin_index][sample_index] = observations.len();
        }
    }

    let totals: Vec<f64> = (0..n_samples)
        .map(|sample_index| {
            (0..n_bins)
                .map(|bin_index| robust[bin_index][sample_index])
                .sum()
        })
        .collect();

    let mut out = Vec::with_capacity(n_bins * n_samples);
    for bin_index in 0..n_bins {
        for sample_index in 0..n_samples {
            let total = totals[sample_index];
            out.push(BinAbundance {
                bin_index,
                sample_index,
                sample_name: coverage.sample_names[sample_index].clone(),
                robust_depth: robust[bin_index][sample_index],
                mean_depth: means[bin_index][sample_index],
                relative_abundance: if total > 0.0 {
                    robust[bin_index][sample_index] / total
                } else {
                    0.0
                },
                covered_bp: covered_bp[bin_index][sample_index],
                covered_contigs: covered_contigs[bin_index][sample_index],
            });
        }
    }
    out
}

pub fn write_abundance_table<P: AsRef<Path>>(
    rows: &[BinAbundance],
    out_dir: P,
) -> io::Result<()> {
    let out_dir = out_dir.as_ref();
    fs::create_dir_all(out_dir)?;
    let mut writer = BufWriter::new(File::create(out_dir.join("abundance.tsv"))?);
    writeln!(
        writer,
        "bin\tsample\trobust_depth\tmean_depth\trelative_abundance\tcovered_bp\tcovered_contigs"
    )?;
    for row in rows {
        writeln!(
            writer,
            "bin_{:04}\t{}\t{:.6}\t{:.6}\t{:.8}\t{}\t{}",
            row.bin_index + 1,
            row.sample_name,
            row.robust_depth,
            row.mean_depth,
            row.relative_abundance,
            row.covered_bp,
            row.covered_contigs
        )?;
    }
    Ok(())
}

fn weighted_median(observations: &mut [(f64, usize)]) -> f64 {
    if observations.is_empty() {
        return 0.0;
    }
    observations.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    let total_weight: usize = observations.iter().map(|(_, weight)| *weight).sum();
    if total_weight == 0 {
        return 0.0;
    }
    let target = total_weight.div_ceil(2);
    let mut cumulative = 0usize;
    for (value, weight) in observations.iter().copied() {
        cumulative = cumulative.saturating_add(weight);
        if cumulative >= target {
            return value;
        }
    }
    observations.last().map(|x| x.0).unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridgebin::{Assignment, BinSummary};
    use std::collections::HashMap;

    #[test]
    fn quantifies_depth_and_relative_abundance() {
        let result = BinningResult {
            assignments: vec![
                Assignment {
                    contig_id: "a".into(),
                    bin_index: Some(0),
                    score: 1.0,
                    length: 2_000,
                },
                Assignment {
                    contig_id: "b".into(),
                    bin_index: Some(0),
                    score: 0.9,
                    length: 1_000,
                },
                Assignment {
                    contig_id: "c".into(),
                    bin_index: Some(1),
                    score: 1.0,
                    length: 3_000,
                },
            ],
            bins: vec![
                BinSummary {
                    bin_index: 0,
                    contig_count: 2,
                    total_bp: 3_000,
                    mean_gc: 0.5,
                },
                BinSummary {
                    bin_index: 1,
                    contig_count: 1,
                    total_bp: 3_000,
                    mean_gc: 0.5,
                },
            ],
        };
        let coverage = CoverageTable {
            sample_names: vec!["s1".into(), "s2".into()],
            values: HashMap::from([
                ("a".into(), vec![10.0, 20.0]),
                ("b".into(), vec![12.0, 18.0]),
                ("c".into(), vec![30.0, 10.0]),
            ]),
        };

        let rows = quantify_bins(&result, &coverage);
        let b0s1 = rows
            .iter()
            .find(|r| r.bin_index == 0 && r.sample_index == 0)
            .unwrap();
        let b1s1 = rows
            .iter()
            .find(|r| r.bin_index == 1 && r.sample_index == 0)
            .unwrap();
        assert!((b0s1.robust_depth - 10.0).abs() < 1e-12);
        assert!((b0s1.mean_depth - (32_000.0 / 3_000.0)).abs() < 1e-12);
        assert!((b0s1.relative_abundance - 0.25).abs() < 1e-12);
        assert!((b1s1.relative_abundance - 0.75).abs() < 1e-12);
    }
}
