from pathlib import Path

path = Path("src/kmer.rs")
text = path.read_text()

old = '''    let mut mercy_support: FxHashMap<KmerKey, u16> = FxHashMap::default();
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
'''
new = '''    let mut mercy_support: FxHashMap<KmerKey, u16> = FxHashMap::default();
    if mercy_max_kmers > 0 {
        for_each_pair(read1, read2, max_pairs, |_index, left, right| {
            // Count support once per physical fragment, even if the same weak
            // k-mer occurs more than once or is present in both mates.
            let mut fragment_candidates: FxHashSet<KmerKey> = FxHashSet::default();
            collect_mercy_candidates(
                &left.sequence,
                k,
                &solid,
                mercy_max_kmers,
                &mut fragment_candidates,
            )?;
            if let Some(right) = right {
                collect_mercy_candidates(
                    &right.sequence,
                    k,
                    &solid,
                    mercy_max_kmers,
                    &mut fragment_candidates,
                )?;
            }
            for key in fragment_candidates {
                let entry = mercy_support.entry(key).or_insert(0);
                *entry = entry.saturating_add(1);
            }
            Ok(())
        })?;
    }
'''
if old not in text:
    raise SystemExit("mercy support block not found")
text = text.replace(old, new, 1)

old = '''fn collect_mercy_candidates(
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
'''
new = '''fn collect_mercy_candidates(
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
'''
if old not in text:
    raise SystemExit("collect_mercy_candidates block not found")
text = text.replace(old, new, 1)

marker = '''    fn evidence_quality_mean() {
        let evidence = KmerEvidence {
            count: 2,
            fragment_count: 2,
            quality_sum: 31 * 30 * 2,
            ..KmerEvidence::default()
        };
        assert!((evidence.mean_quality(31) - 30.0).abs() < 1e-6);
    }
'''
addition = marker + '''

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
'''
if marker not in text:
    raise SystemExit("test insertion marker not found")
text = text.replace(marker, addition, 1)

path.write_text(text)
