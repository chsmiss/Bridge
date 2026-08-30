from pathlib import Path

main = Path("src/main.rs")
text = main.read_text()
text = text.replace(
'''        #[arg(long, default_value_t = 1)]
        mercy_min_support: u16,
''',
'''        #[arg(long, default_value_t = 1)]
        mercy_min_support: u16,
        #[arg(long, default_value_t = 25.0)]
        mercy_min_quality: f32,
''',
1,
)
text = text.replace(
'''            mercy_min_support,
            min_read_support,
''',
'''            mercy_min_support,
            mercy_min_quality,
            min_read_support,
''',
1,
)
text = text.replace(
'''            if k == 0 || k > MAX_K {
                anyhow::bail!("k must be in 1..={MAX_K}");
            }
''',
'''            if k == 0 || k > MAX_K {
                anyhow::bail!("k must be in 1..={MAX_K}");
            }
            if !(0.0..=60.0).contains(&mercy_min_quality) {
                anyhow::bail!("mercy minimum quality must be in 0..=60");
            }
''',
1,
)
text = text.replace(
'''                mercy_min_support,
                min_read_support,
''',
'''                mercy_min_support,
                mercy_min_quality,
                min_read_support,
''',
1,
)
main.write_text(text)

assembler = Path("src/assembler.rs")
text = assembler.read_text()
text = text.replace(
'''use crate::kmer::{count_and_filter, KmerCountSummary};
''',
'''use crate::kmer::{count_and_filter, KmerCountSummary, KmerFilterConfig};
''',
1,
)
text = text.replace(
'''    pub mercy_min_support: u16,
    pub min_read_support: u32,
''',
'''    pub mercy_min_support: u16,
    pub mercy_min_quality: f32,
    pub min_read_support: u32,
''',
1,
)
old_call = '''    let kmer_set = count_and_filter(
        &config.read1,
        config.read2.as_deref(),
        config.k,
        config.min_count,
        config.mercy_max_kmers,
        config.mercy_min_support,
        config.max_pairs,
    )?;
'''
new_call = '''    let kmer_set = count_and_filter(
        &config.read1,
        config.read2.as_deref(),
        KmerFilterConfig {
            k: config.k,
            min_count: config.min_count,
            mercy_max_kmers: config.mercy_max_kmers,
            mercy_min_support: config.mercy_min_support,
            mercy_min_quality: config.mercy_min_quality,
            max_pairs: config.max_pairs,
        },
    )?;
'''
if old_call not in text:
    raise SystemExit("count_and_filter call not found")
text = text.replace(old_call, new_call, 1)
assembler.write_text(text)

kmer = Path("src/kmer.rs")
text = kmer.read_text()
marker = '''#[derive(Debug)]
pub struct KmerSet {
    pub evidence: FxHashMap<KmerKey, KmerEvidence>,
    pub retained: FxHashSet<KmerKey>,
    pub rescued: FxHashSet<KmerKey>,
    pub summary: KmerCountSummary,
}

'''
insert = marker + '''#[derive(Clone, Copy, Debug)]
pub struct KmerFilterConfig {
    pub k: usize,
    pub min_count: u32,
    pub mercy_max_kmers: usize,
    pub mercy_min_support: u16,
    pub mercy_min_quality: f32,
    pub max_pairs: Option<usize>,
}

'''
if marker not in text:
    raise SystemExit("KmerSet marker not found")
text = text.replace(marker, insert, 1)
old_signature = '''pub fn count_and_filter(
    read1: &Path,
    read2: Option<&Path>,
    k: usize,
    min_count: u32,
    mercy_max_kmers: usize,
    mercy_min_support: u16,
    max_pairs: Option<usize>,
) -> Result<KmerSet> {
'''
new_signature = '''pub fn count_and_filter(
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
'''
if old_signature not in text:
    raise SystemExit("count_and_filter signature not found")
text = text.replace(old_signature, new_signature, 1)
text = text.replace(
'''    let mercy_min_quality = if min_count <= 1 { 0.0 } else { 25.0 };
''',
'''    let mercy_min_quality = if min_count <= 1 {
        0.0
    } else {
        mercy_min_quality
    };
''',
1,
)
kmer.write_text(text)

tests = Path("tests/integration.rs")
text = tests.read_text()
text = text.replace(
'''        mercy_min_support: 1,
        min_read_support: 1,
''',
'''        mercy_min_support: 1,
        mercy_min_quality: 25.0,
        min_read_support: 1,
''',
1,
)
tests.write_text(text)
