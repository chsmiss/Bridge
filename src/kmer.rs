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
    pub min_count: u32,
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

#[derive(Clone, Copy, Debug)]
pub struct KmerFilterConfig {
    pub k: usize,
    pub min_count: u32,
    pub mercy_max_kmers: usize,
    pub mercy_min_support: u16,
    pub mercy_min_quality: f32,
    pub max_pairs: Option<usize>,
}

pub fn count_and_filter(
    read1: &Path,
    read2: Option<&Path>,
    config: KmerFilterConfig,
) -> Result<KmerSet> {
    let KmerFilterConfig {
        k,
        min_count,
        mercy_max_kmers,
        mercy_min_support,
        mercy_min_quality,
        max_pairs,
    } = config;
    // Keep the public API stable while making the production min-count=2 path
    // quality- and independent-fragment-aware. Tests and explicitly permissive
    // runs using min-count=1 retain the historical behavior.
    let min_fragment_support = if min_count <= 1 { 1 } else { 2 };
    let min_mean_quality = if min_count <= 1 { 0.0 } else { 20.0 };
    let mercy_min_quality = if min_count <= 1 {
        0.0
    } else {
        mercy_min_quality
    };

    // Experimental, reference-free recovery gate. This is intentionally an
    // environment knob while the idea is benchmarked; zero preserves the
    // production behavior exactly. A weak terminal chain is considered only
    // when the opposite mate contains at least one solid k-mer. The ordinary
    // mercy quality/support filters are still applied afterwards.
    let mate_terminal_mercy_kmers = std::env::var("BRIDGEASM_MATE_TERMINAL_MERCY_KMERS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);

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
    if mercy_max_kmers > 0 || mate_terminal_mercy_kmers > 0 {
        for_each_pair(read1, read2, max_pairs, |_index, left, right| {
            // Count support once per physical fragment, even if the same weak
            // k-mer occurs more than once or is present in both mates.
            let mut fragment_candidates: FxHashSet<KmerKey> = FxHashSet::default();
            if mercy_max_kmers > 0 {
                collect_mercy_candidates(
                    &left.sequence,
                    k,
                    &solid,
                    mercy_max_kmers,
                    &mut fragment_candidates,
                )?;
                if let Some(right) = right.as_ref() {
                    collect_mercy_candidates(
                        &right.sequence,
                        k,
                        &solid,
                        mercy_max_kmers,
                        &mut fragment_candidates,
                    )?;
                }
            }

            if mate_terminal_mercy_kmers > 0 {
                if let Some(right) = right.as_ref() {
                    let left_has_solid = record_has_solid(&left.sequence, k, &solid)?;
                    let right_has_solid = record_has_solid(&right.sequence, k, &solid)?;
                    if right_has_solid {
                        collect_mate_anchored_candidates(
                            &left.sequence,
                            k,
                            &solid,
                            mate_terminal_mercy_kmers,
                            &mut fragment_candidates,
                        )?;
                    }
                    if left_has_solid {
                        collect_mate_anchored_candidates(
                            &right.sequence,
                            k,
                            &solid,
                            mate_terminal_mercy_kmers,
                            &mut fragment_candidates,
                        )?;
                    }
                }
            }

            for key in fragment_candidates {
                let entry = mercy_support.entry(key).or_insert(0);
                *entry = entry.saturating_add(1);
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
        min_count,
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
    candidates: &mut FxHashSet<KmerKey>,
) -> Result<()> {
    let kmers = canonical_kmers(sequence, k)?;
    if kmers.len() < 3 {
        return Ok(());
    }

    let mut left_solid: Option<usize> = None;
    let mut previous_position: Option<usize> = None;
    for (index, item) in kmers.iter().enumerate() {
        // canonical_kmers skips windows containing N or another invalid base.
        // A mercy bridge must stay within one continuous read interval.
        if previous_position.is_some_and(|position| item.position != position + 1) {
            left_solid = None;
        }
        if solid.contains(&item.key) {
            if let Some(left) = left_solid {
                let weak_count = index.saturating_sub(left + 1);
                let continuous_span = item.position - kmers[left].position == index - left;
                if continuous_span && weak_count > 0 && weak_count <= mercy_max_kmers {
                    candidates.extend(kmers[left + 1..index].iter().map(|weak| weak.key));
                }
            }
            left_solid = Some(index);
        }
        previous_position = Some(item.position);
    }
    Ok(())
}

fn record_has_solid(sequence: &[u8], k: usize, solid: &FxHashSet<KmerKey>) -> Result<bool> {
    Ok(canonical_kmers(sequence, k)?
        .iter()
        .any(|item| solid.contains(&item.key)))
}

fn collect_mate_anchored_candidates(
    sequence: &[u8],
    k: usize,
    solid: &FxHashSet<KmerKey>,
    max_weak_kmers: usize,
    candidates: &mut FxHashSet<KmerKey>,
) -> Result<()> {
    if max_weak_kmers == 0 {
        return Ok(());
    }
    let kmers = canonical_kmers(sequence, k)?;
    if kmers.is_empty() {
        return Ok(());
    }

    // Treat stretches separated by N/invalid bases independently. Within one
    // continuous run, rescue only terminal weak chains adjacent to the first
    // or last solid k-mer. If a whole run has no local solid anchor, it may be
    // rescued only when the entire run fits the configured bound; the caller
    // has already required a solid anchor in the opposite mate.
    let mut run_start = 0_usize;
    while run_start < kmers.len() {
        let mut run_end = run_start + 1;
        while run_end < kmers.len() && kmers[run_end].position == kmers[run_end - 1].position + 1 {
            run_end += 1;
        }

        let first_solid = (run_start..run_end).find(|&index| solid.contains(&kmers[index].key));
        let last_solid = (run_start..run_end)
            .rev()
            .find(|&index| solid.contains(&kmers[index].key));

        match (first_solid, last_solid) {
            (Some(first), Some(last)) => {
                let left_start = first.saturating_sub(max_weak_kmers).max(run_start);
                candidates.extend(kmers[left_start..first].iter().map(|item| item.key));
                let right_end = (last + 1 + max_weak_kmers).min(run_end);
                candidates.extend(kmers[last + 1..right_end].iter().map(|item| item.key));
            }
            (None, None) if run_end - run_start <= max_weak_kmers => {
                candidates.extend(kmers[run_start..run_end].iter().map(|item| item.key));
            }
            _ => {}
        }
        run_start = run_end;
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

    #[test]
    fn mercy_does_not_bridge_across_invalid_bases() {
        let sequence = b"AAAAACCCCCNNNNNGGGGGTTTTT";
        let k = 5;
        let kmers = canonical_kmers(sequence, k).unwrap();
        let mut solid = FxHashSet::default();
        solid.insert(kmers.first().unwrap().key);
        solid.insert(kmers.last().unwrap().key);
        let mut candidates = FxHashSet::default();
        collect_mercy_candidates(sequence, k, &solid, 100, &mut candidates).unwrap();
        assert!(candidates.is_empty());
    }

    #[test]
    fn mercy_candidates_are_unique_within_one_fragment() {
        let sequence = b"AAAAACAAAACAAAACAAAAA";
        let k = 5;
        let kmers = canonical_kmers(sequence, k).unwrap();
        let mut solid = FxHashSet::default();
        solid.insert(kmers.first().unwrap().key);
        solid.insert(kmers.last().unwrap().key);
        let mut candidates = FxHashSet::default();
        collect_mercy_candidates(sequence, k, &solid, 100, &mut candidates).unwrap();
        let unique_from_slice: FxHashSet<_> = kmers[1..kmers.len() - 1]
            .iter()
            .map(|item| item.key)
            .collect();
        assert_eq!(candidates, unique_from_slice);
    }

    #[test]
    fn mate_anchored_candidates_recover_terminal_weak_chain_only() {
        let sequence = b"AAAACCCCGGGGTTTT";
        let k = 5;
        let kmers = canonical_kmers(sequence, k).unwrap();
        let mut solid = FxHashSet::default();
        solid.insert(kmers[4].key);
        solid.insert(kmers[7].key);
        let mut candidates = FxHashSet::default();
        collect_mate_anchored_candidates(sequence, k, &solid, 2, &mut candidates).unwrap();
        let expected: FxHashSet<_> = [kmers[2].key, kmers[3].key, kmers[8].key, kmers[9].key]
            .into_iter()
            .collect();
        assert_eq!(candidates, expected);
    }

    #[test]
    fn mate_anchored_candidates_can_rescue_short_fully_weak_run() {
        let sequence = b"AAAACCCCGG";
        let k = 5;
        let solid = FxHashSet::default();
        let kmers = canonical_kmers(sequence, k).unwrap();
        let mut candidates = FxHashSet::default();
        collect_mate_anchored_candidates(sequence, k, &solid, kmers.len(), &mut candidates)
            .unwrap();
        let expected: FxHashSet<_> = kmers.iter().map(|item| item.key).collect();
        assert_eq!(candidates, expected);
    }
}
