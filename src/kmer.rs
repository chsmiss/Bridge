use crate::dna::{canonical_kmers, KmerKey};
use crate::fastq::for_each_pair;
use anyhow::Result;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::path::Path;

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct KmerEvidence {
    pub count: u32,
    pub fragment_count: u32,
    pub forward_count: u32,
    pub reverse_count: u32,
    pub quality_sum: u64,
}

impl KmerEvidence {
    pub fn mean_quality(self, k: usize) -> f32 {
        if self.count == 0 || k == 0 {
            return 0.0;
        }
        self.quality_sum as f32 / (self.count as f32 * k as f32)
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct KmerCountSummary {
    pub k: usize,
    pub read_pairs: usize,
    pub observations: u64,
    pub distinct: usize,
    pub abundance_candidates: usize,
    pub solid: usize,
    pub rescued: usize,
}

#[derive(Debug)]
pub struct KmerSet {
    pub evidence: FxHashMap<KmerKey, KmerEvidence>,
    pub retained: FxHashSet<KmerKey>,
    pub rescued: FxHashSet<KmerKey>,
    pub summary: KmerCountSummary,
}

pub fn count_and_filter(
    read1: &Path,
    read2: Option<&Path>,
    k: usize,
    min_count: u32,
    mercy_max_kmers: usize,
    mercy_min_support: u16,
    max_pairs: Option<usize>,
) -> Result<KmerSet> {
    // Keep the public API stable while making the production min-count=2 path
    // quality- and independent-fragment-aware. Tests and explicitly permissive
    // runs using min-count=1 retain the historical behavior.
    let min_fragment_support = if min_count <= 1 { 1 } else { 2 };
    let min_mean_quality = if min_count <= 1 { 0.0 } else { 20.0 };
    let mercy_min_quality = if min_count <= 1 { 0.0 } else { 25.0 };

    let mut evidence: FxHashMap<KmerKey, KmerEvidence> = FxHashMap::default();
    let mut observations = 0_u64;

    let read_pairs = for_each_pair(read1, read2, max_pairs, |_index, left, right| {
        let mut fragment_keys = Vec::with_capacity(
            left.sequence.len().saturating_sub(k).saturating_add(1)
                * if right.is_some() { 2 } else { 1 },
        );
        count_record(
            &left.sequence,
            &left.quality,
            k,
            &mut evidence,
            &mut observations,
            &mut fragment_keys,
        )?;
        if let Some(right) = right {
            count_record(
                &right.sequence,
                &right.quality,
                k,
                &mut evidence,
                &mut observations,
                &mut fragment_keys,
            )?;
        }
        fragment_keys.sort_unstable();
        fragment_keys.dedup();
        for key in fragment_keys {
            if let Some(entry) = evidence.get_mut(&key) {
                entry.fragment_count = entry.fragment_count.saturating_add(1);
            }
        }
        Ok(())
    })?;

    let abundance_candidates = evidence
        .values()
        .filter(|value| value.count >= min_count)
        .count();
    let solid: FxHashSet<KmerKey> = evidence
        .iter()
        .filter_map(|(key, value)| {
            (value.count >= min_count
                && value.fragment_count >= min_fragment_support
                && value.mean_quality(k) >= min_mean_quality)
                .then_some(*key)
        })
        .collect();

    let mut mercy_support: FxHashMap<KmerKey, u16> = FxHashMap::default();
    if mercy_max_kmers > 0 {
        for_each_pair(read1, read2, max_pairs, |_index, left, right| {
            collect_mercy_candidates(
                &left.sequence,
                k,
                &solid,
                mercy_max_kmers,
                &mut mercy_support,
            )?;
            if let Some(right) = right {
                collect_mercy_candidates(
                    &right.sequence,
                    k,
                    &solid,
                    mercy_max_kmers,
                    &mut mercy_support,
                )?;
            }
            Ok(())
        })?;
    }

    let rescued: FxHashSet<KmerKey> = mercy_support
        .into_iter()
        .filter_map(|(key, support)| {
            let value = evidence.get(&key)?;
            (support >= mercy_min_support && value.mean_quality(k) >= mercy_min_quality)
                .then_some(key)
        })
        .collect();
    let mut retained = solid;
    retained.extend(rescued.iter().copied());

    let summary = KmerCountSummary {
        k,
        read_pairs,
        observations,
        distinct: evidence.len(),
        abundance_candidates,
        solid: retained.len().saturating_sub(rescued.len()),
        rescued: rescued.len(),
    };

    Ok(KmerSet {
        evidence,
        retained,
        rescued,
        summary,
    })
}

fn count_record(
    sequence: &[u8],
    quality: &[u8],
    k: usize,
    evidence: &mut FxHashMap<KmerKey, KmerEvidence>,
    observations: &mut u64,
    fragment_keys: &mut Vec<KmerKey>,
) -> Result<()> {
    let kmers = canonical_kmers(sequence, k)?;
    let mut quality_prefix = vec![0_u64; quality.len() + 1];
    for (index, value) in quality.iter().enumerate() {
        let phred = value.saturating_sub(33) as u64;
        quality_prefix[index + 1] = quality_prefix[index] + phred;
    }

    for item in kmers {
        let window_quality = quality_prefix[item.position + k] - quality_prefix[item.position];
        let entry = evidence.entry(item.key).or_default();
        entry.count = entry.count.saturating_add(1);
        if item.reverse {
            entry.reverse_count = entry.reverse_count.saturating_add(1);
        } else {
            entry.forward_count = entry.forward_count.saturating_add(1);
        }
        entry.quality_sum = entry.quality_sum.saturating_add(window_quality);
        fragment_keys.push(item.key);
        *observations += 1;
    }
    Ok(())
}

fn collect_mercy_candidates(
    sequence: &[u8],
    k: usize,
    solid: &FxHashSet<KmerKey>,
    mercy_max_kmers: usize,
    support: &mut FxHashMap<KmerKey, u16>,
) -> Result<()> {
    let kmers = canonical_kmers(sequence, k)?;
    if kmers.len() < 3 {
        return Ok(());
    }

    let mut left_solid: Option<usize> = None;
    for (index, item) in kmers.iter().enumerate() {
        if solid.contains(&item.key) {
            if let Some(left) = left_solid {
                let weak_count = index.saturating_sub(left + 1);
                if weak_count > 0 && weak_count <= mercy_max_kmers {
                    for weak in &kmers[left + 1..index] {
                        let entry = support.entry(weak.key).or_insert(0);
                        *entry = entry.saturating_add(1);
                    }
                }
            }
            left_solid = Some(index);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evidence_quality_mean() {
        let evidence = KmerEvidence {
            count: 2,
            fragment_count: 2,
            quality_sum: 31 * 30 * 2,
            ..KmerEvidence::default()
        };
        assert!((evidence.mean_quality(31) - 30.0).abs() < 1e-6);
    }
}
