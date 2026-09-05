use anyhow::{bail, Result};
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

/// Illumina 150 bp reads leave at least three start positions at k=147.
pub const MAX_K: usize = 147;
const WORDS: usize = 5;

#[inline]
pub fn base_bits(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
pub fn bits_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

#[inline]
pub fn complement_bits(bits: u8) -> u8 {
    3 - (bits & 0b11)
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Exact packed DNA key supporting k <= 147. Bits above 2*k are always zero.
/// Words are little-endian: words[0] stores the least-significant 64 bits.
#[derive(Clone, Copy, Default, Eq, PartialEq)]
pub struct KmerKey {
    pub words: [u64; WORDS],
}

impl Hash for KmerKey {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Most production short-read k-mers use only one or two of the five
        // words. Hashing all five words wastes memory bandwidth on guaranteed
        // trailing zeroes. Equal keys always have the same highest non-zero
        // word, so omitting those zero suffixes preserves the Hash contract.
        let mut used = WORDS;
        while used > 1 && self.words[used - 1] == 0 {
            used -= 1;
        }
        for word in &self.words[..used] {
            state.write_u64(*word);
        }
    }
}

impl Ord for KmerKey {
    fn cmp(&self, other: &Self) -> Ordering {
        for index in (0..WORDS).rev() {
            match self.words[index].cmp(&other.words[index]) {
                Ordering::Equal => {}
                ordering => return ordering,
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for KmerKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Debug for KmerKey {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("KmerKey(")?;
        for word in self.words.iter().rev() {
            write!(formatter, "{word:016x}")?;
        }
        formatter.write_str(")")
    }
}

impl KmerKey {
    pub const ZERO: Self = Self { words: [0; WORDS] };

    pub fn from_sequence(sequence: &[u8]) -> Result<Self> {
        if sequence.is_empty() || sequence.len() > MAX_K {
            bail!("k-mer length must be in 1..={MAX_K}");
        }
        let mut value = Self::ZERO;
        for &base in sequence {
            let bits = base_bits(base)
                .ok_or_else(|| anyhow::anyhow!("invalid DNA base in k-mer: {}", base as char))?;
            value.shift_left_two();
            value.words[0] |= u64::from(bits);
        }
        value.mask_to_k(sequence.len());
        Ok(value)
    }

    pub fn reverse_complement(self, k: usize) -> Self {
        let mut output = Self::ZERO;
        for position in 0..k {
            let source_offset = 2 * position;
            let bits = self.get_two_bits(source_offset);
            output.shift_left_two();
            output.words[0] |= u64::from(complement_bits(bits));
        }
        output.mask_to_k(k);
        output
    }

    pub fn canonical(self, k: usize) -> (Self, bool) {
        let reverse = self.reverse_complement(k);
        if reverse < self {
            (reverse, true)
        } else {
            (self, false)
        }
    }

    pub fn to_sequence(self, k: usize) -> Vec<u8> {
        let mut sequence = vec![b'A'; k];
        for (position, base) in sequence.iter_mut().enumerate() {
            let offset = 2 * (k - 1 - position);
            *base = bits_base(self.get_two_bits(offset));
        }
        sequence
    }

    #[inline]
    pub fn last_base(self) -> u8 {
        bits_base((self.words[0] & 0b11) as u8)
    }

    #[inline]
    pub fn first_base(self, k: usize) -> u8 {
        bits_base(self.get_two_bits(2 * (k - 1)))
    }

    #[inline]
    pub fn shift_left_append(&mut self, bits: u8, k: usize) {
        self.shift_left_two();
        self.words[0] |= u64::from(bits & 0b11);
        self.mask_to_k(k);
    }

    #[inline]
    pub fn shift_right_prepend_complement(&mut self, bits: u8, k: usize) {
        self.shift_right_two();
        self.set_two_bits(2 * (k - 1), complement_bits(bits));
        self.mask_to_k(k);
    }

    #[inline]
    fn shift_left_two(&mut self) {
        let mut carry = 0_u64;
        for word in &mut self.words {
            let next_carry = *word >> 62;
            *word = (*word << 2) | carry;
            carry = next_carry;
        }
    }

    #[inline]
    fn shift_right_two(&mut self) {
        let mut carry = 0_u64;
        for word in self.words.iter_mut().rev() {
            let next_carry = *word & 0b11;
            *word = (*word >> 2) | (carry << 62);
            carry = next_carry;
        }
    }

    #[inline]
    fn get_two_bits(self, offset: usize) -> u8 {
        let word = offset / 64;
        let shift = offset % 64;
        if shift <= 62 {
            ((self.words[word] >> shift) & 0b11) as u8
        } else {
            let low = (self.words[word] >> 63) & 0b1;
            let high = if word + 1 < WORDS {
                (self.words[word + 1] & 0b1) << 1
            } else {
                0
            };
            (low | high) as u8
        }
    }

    #[inline]
    fn set_two_bits(&mut self, offset: usize, bits: u8) {
        let word = offset / 64;
        let shift = offset % 64;
        let bits = u64::from(bits & 0b11);
        if shift <= 62 {
            self.words[word] &= !(0b11_u64 << shift);
            self.words[word] |= bits << shift;
        } else {
            self.words[word] &= !(1_u64 << 63);
            self.words[word] |= (bits & 1) << 63;
            if word + 1 < WORDS {
                self.words[word + 1] &= !1_u64;
                self.words[word + 1] |= (bits >> 1) & 1;
            }
        }
    }

    #[inline]
    fn mask_to_k(&mut self, k: usize) {
        let used_bits = 2 * k;
        let full_words = used_bits / 64;
        let remainder = used_bits % 64;
        if remainder == 0 {
            for word in self.words.iter_mut().skip(full_words) {
                *word = 0;
            }
        } else {
            let mask = (1_u64 << remainder) - 1;
            self.words[full_words] &= mask;
            for word in self.words.iter_mut().skip(full_words + 1) {
                *word = 0;
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrientedKmer {
    pub key: KmerKey,
    /// false means the read-forward encoding equals the canonical key;
    /// true means the read traverses the reverse-complement orientation.
    pub reverse: bool,
    pub position: usize,
}

/// Visit exact canonical k-mers in O(read_length) time without allocating a
/// per-read output vector. Invalid bases reset the roller.
pub fn for_each_canonical_kmer<F>(sequence: &[u8], k: usize, mut callback: F) -> Result<usize>
where
    F: FnMut(OrientedKmer),
{
    if k == 0 || k > MAX_K {
        bail!("k must be in 1..={MAX_K}");
    }
    if sequence.len() < k {
        return Ok(0);
    }

    let mut forward = KmerKey::ZERO;
    let mut reverse = KmerKey::ZERO;
    let mut valid = 0_usize;
    let mut emitted = 0_usize;

    for (index, &base) in sequence.iter().enumerate() {
        let Some(bits) = base_bits(base) else {
            forward = KmerKey::ZERO;
            reverse = KmerKey::ZERO;
            valid = 0;
            continue;
        };

        if valid < k {
            forward.shift_left_append(bits, k);
            reverse.set_two_bits(2 * valid, complement_bits(bits));
            valid += 1;
        } else {
            forward.shift_left_append(bits, k);
            reverse.shift_right_prepend_complement(bits, k);
        }

        if valid >= k {
            let position = index + 1 - k;
            if reverse < forward {
                callback(OrientedKmer {
                    key: reverse,
                    reverse: true,
                    position,
                });
            } else {
                callback(OrientedKmer {
                    key: forward,
                    reverse: false,
                    position,
                });
            }
            emitted += 1;
        }
    }

    Ok(emitted)
}

/// Collect exact canonical k-mers. Kept for algorithms that need random access
/// to neighboring k-mers; hot counting paths should prefer for_each_canonical_kmer.
pub fn canonical_kmers(sequence: &[u8], k: usize) -> Result<Vec<OrientedKmer>> {
    let mut output = Vec::with_capacity(sequence.len().saturating_sub(k).saturating_add(1));
    for_each_canonical_kmer(sequence, k, |item| output.push(item))?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_round_trip() {
        let sequence = b"ACGTTGCA";
        assert_eq!(reverse_complement(&reverse_complement(sequence)), sequence);
    }

    #[test]
    fn packed_key_supports_long_kmers() {
        for k in [1, 21, 31, 41, 63, 71, 91, 111, 127, 147] {
            let sequence: Vec<u8> = (0..k).map(|index| b"ACGT"[index % 4]).collect();
            let key = KmerKey::from_sequence(&sequence).unwrap();
            assert_eq!(key.to_sequence(k), sequence);
            let rc = key.reverse_complement(k);
            assert_eq!(rc.to_sequence(k), reverse_complement(&sequence));
        }
    }

    #[test]
    fn rolling_matches_direct_encoding() {
        let sequence: Vec<u8> = (0..300)
            .map(|index| b"ACGTGCACT"[(index * 7 + index / 5) % 9])
            .collect();
        for k in [21, 31, 41, 91, 127, 147] {
            let rolling = canonical_kmers(&sequence, k).unwrap();
            for item in rolling {
                let direct = KmerKey::from_sequence(&sequence[item.position..item.position + k])
                    .unwrap()
                    .canonical(k);
                assert_eq!(item.key, direct.0);
                assert_eq!(item.reverse, direct.1);
            }
        }
    }

    #[test]
    fn visitor_matches_collector() {
        let sequence = b"ACGTACGTNNACGTACGTACGT";
        for k in [3, 5, 7] {
            let collected = canonical_kmers(sequence, k).unwrap();
            let mut visited = Vec::new();
            let emitted = for_each_canonical_kmer(sequence, k, |item| visited.push(item)).unwrap();
            assert_eq!(emitted, collected.len());
            assert_eq!(visited, collected);
        }
    }
}
