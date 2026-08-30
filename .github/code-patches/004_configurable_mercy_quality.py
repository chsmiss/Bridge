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
'''    pub mercy_min_support: u16,
    pub min_read_support: u32,
''',
'''    pub mercy_min_support: u16,
    pub mercy_min_quality: f32,
    pub min_read_support: u32,
''',
1,
)
text = text.replace(
'''        config.mercy_min_support,
        config.max_pairs,
''',
'''        config.mercy_min_support,
        config.mercy_min_quality,
        config.max_pairs,
''',
1,
)
assembler.write_text(text)

kmer = Path("src/kmer.rs")
text = kmer.read_text()
text = text.replace(
'''    mercy_max_kmers: usize,
    mercy_min_support: u16,
    max_pairs: Option<usize>,
''',
'''    mercy_max_kmers: usize,
    mercy_min_support: u16,
    mercy_min_quality: f32,
    max_pairs: Option<usize>,
''',
1,
)
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
