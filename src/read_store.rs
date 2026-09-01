use crate::dna::{base_bits, bits_base};
use crate::fastq::{for_each_pair, FastqRecord, PairSource};
use anyhow::Result;
use std::path::Path;

#[derive(Debug, Default)]
pub struct PackedReadStore {
    bases: Vec<u64>,
    valid: Vec<u64>,
    qualities: Vec<u8>,
    offsets: Vec<u64>,
    lengths: Vec<u32>,
    pair_count: usize,
    paired: bool,
    total_bases: usize,
}

impl PackedReadStore {
    pub fn load(
        read1: &Path,
        read2: Option<&Path>,
        max_pairs: Option<usize>,
    ) -> Result<Self> {
        let mut store = Self {
            paired: read2.is_some(),
            ..Self::default()
        };
        let pair_count = for_each_pair(read1, read2, max_pairs, |_index, left, right| {
            store.push_record(&left.sequence, &left.quality);
            if let Some(right) = right {
                store.push_record(&right.sequence, &right.quality);
            }
            Ok(())
        })?;
        store.pair_count = pair_count;
        Ok(store)
    }

    fn push_record(&mut self, sequence: &[u8], quality: &[u8]) {
        debug_assert_eq!(sequence.len(), quality.len());
        self.offsets.push(self.total_bases as u64);
        self.lengths.push(sequence.len() as u32);
        self.qualities.extend_from_slice(quality);
        for &base in sequence {
            let index = self.total_bases;
            let base_word = (index * 2) >> 6;
            let base_shift = (index * 2) & 63;
            if self.bases.len() <= base_word {
                self.bases.push(0);
            }
            let valid_word = index >> 6;
            let valid_shift = index & 63;
            if self.valid.len() <= valid_word {
                self.valid.push(0);
            }
            if let Some(bits) = base_bits(base) {
                self.bases[base_word] |= u64::from(bits) << base_shift;
                self.valid[valid_word] |= 1_u64 << valid_shift;
            }
            self.total_bases += 1;
        }
    }

    #[inline]
    fn base_at(&self, index: usize) -> u8 {
        let valid_word = index >> 6;
        let valid_shift = index & 63;
        if ((self.valid[valid_word] >> valid_shift) & 1) == 0 {
            return b'N';
        }
        let base_word = (index * 2) >> 6;
        let base_shift = (index * 2) & 63;
        bits_base(((self.bases[base_word] >> base_shift) & 3) as u8)
    }

    fn decode_record(&self, record_index: usize) -> FastqRecord {
        let start = self.offsets[record_index] as usize;
        let length = self.lengths[record_index] as usize;
        let end = start + length;
        let mut sequence = Vec::with_capacity(length);
        for index in start..end {
            sequence.push(self.base_at(index));
        }
        FastqRecord {
            id: Vec::new(),
            sequence,
            quality: self.qualities[start..end].to_vec(),
        }
    }

    pub fn pair_count(&self) -> usize {
        self.pair_count
    }

    pub fn record_count(&self) -> usize {
        self.lengths.len()
    }

    pub fn total_bases(&self) -> usize {
        self.total_bases
    }

    pub fn packed_bytes(&self) -> usize {
        self.bases.capacity() * std::mem::size_of::<u64>()
            + self.valid.capacity() * std::mem::size_of::<u64>()
            + self.qualities.capacity() * std::mem::size_of::<u8>()
            + self.offsets.capacity() * std::mem::size_of::<u64>()
            + self.lengths.capacity() * std::mem::size_of::<u32>()
    }
}

impl PairSource for PackedReadStore {
    fn for_each_pair_dyn(
        &self,
        max_pairs: Option<usize>,
        callback: &mut dyn FnMut(usize, FastqRecord, Option<FastqRecord>) -> Result<()>,
    ) -> Result<usize> {
        let limit = max_pairs.unwrap_or(self.pair_count).min(self.pair_count);
        for pair_index in 0..limit {
            if self.paired {
                let left = self.decode_record(pair_index * 2);
                let right = self.decode_record(pair_index * 2 + 1);
                callback(pair_index, left, Some(right))?;
            } else {
                callback(pair_index, self.decode_record(pair_index), None)?;
            }
        }
        Ok(limit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn round_trips_paired_reads_and_invalid_bases() {
        let mut left = tempfile::NamedTempFile::new().unwrap();
        let mut right = tempfile::NamedTempFile::new().unwrap();
        writeln!(left, "@l1\nACNT\n+\nIJKL").unwrap();
        writeln!(right, "@r1\nTGCA\n+\nMNOP").unwrap();
        let store = PackedReadStore::load(left.path(), Some(right.path()), None).unwrap();
        assert_eq!(store.pair_count(), 1);
        assert_eq!(store.record_count(), 2);
        let mut seen = Vec::new();
        let mut callback = |_index, left: FastqRecord, right: Option<FastqRecord>| {
            seen.push((left.sequence, left.quality));
            let right = right.unwrap();
            seen.push((right.sequence, right.quality));
            Ok(())
        };
        store.for_each_pair_dyn(None, &mut callback).unwrap();
        assert_eq!(seen[0], (b"ACNT".to_vec(), b"IJKL".to_vec()));
        assert_eq!(seen[1], (b"TGCA".to_vec(), b"MNOP".to_vec()));
    }
}
