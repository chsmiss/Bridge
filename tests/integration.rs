use bridgeasm::assembler::{assemble, AssembleConfig};
use std::fs::File;
use std::io::Write;

fn write_reads(path: &std::path::Path, genome: &[u8], read_length: usize, step: usize) {
    let mut file = File::create(path).unwrap();
    append_reads(&mut file, "r", genome, read_length, step, 1);
}

fn append_reads(
    file: &mut File,
    label: &str,
    genome: &[u8],
    read_length: usize,
    step: usize,
    copies: usize,
) {
    for copy in 0..copies {
        for (index, start) in (0..=genome.len() - read_length).step_by(step).enumerate() {
            let sequence = &genome[start..start + read_length];
            writeln!(file, "@{label}_{copy}_{index}").unwrap();
            writeln!(file, "{}", String::from_utf8_lossy(sequence)).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(sequence.len())).unwrap();
        }
    }
}

fn config(reads: std::path::PathBuf, output: std::path::PathBuf, k: usize) -> AssembleConfig {
    AssembleConfig {
        read1: reads,
        read2: None,
        output_dir: output,
        k,
        min_count: 1,
        mercy_max_kmers: 0,
        mercy_min_support: 1,
        min_read_support: 1,
        min_pair_support: 1,
        min_primary_support: 2,
        primary_dominance: 0.70,
        min_contig_length: 20,
        max_pairs: None,
        threads: 1,
    }
}

#[test]
fn assembles_a_linear_genome() {
    let directory = tempfile::tempdir().unwrap();
    let reads = directory.path().join("reads.fastq");
    let output = directory.path().join("out");
    let genome = b"ACGTTGCAAGTCGATCGTACCTGACTGATCGTAGCTAGCTACGATCGATGCTAGCATCGATCGTACGATGCTAGCTAGCATGCTAGCATCGATCGTAGCTA";
    write_reads(&reads, genome, 60, 5);

    let product = assemble(&config(reads, output, 21)).unwrap();

    assert!(product.stats.primary_contigs >= 1);
    assert!(product.stats.largest_primary >= genome.len() - 5);
}

#[test]
fn supports_k_greater_than_31() {
    let directory = tempfile::tempdir().unwrap();
    let reads = directory.path().join("reads.fastq");
    let output = directory.path().join("out");
    let genome = b"ACGTTGCAAGTCGATCGTACCTGACTGATCGTAGCTAGCTACGATCGATGCTAGCATCGATCGTACGATGCTAGCTAGCATGCTAGCATCGATCGTAGCTAACGATCGTACGATCG";
    write_reads(&reads, genome, 90, 3);

    let product = assemble(&config(reads, output, 41)).unwrap();

    assert!(product.stats.primary_contigs >= 1);
    assert!(product.stats.graph.canonical_nodes > 0);
}

#[test]
fn preserves_a_supported_strain_bubble() {
    let directory = tempfile::tempdir().unwrap();
    let reads = directory.path().join("mixture.fastq");
    let output = directory.path().join("out");
    let major = b"ACGTTGCAAGTCGATCGTACCTGACTGATCGTAGCTAGCTACGATCGATGCTAGCATCGATCGTACGATGCTAGCTAGCATGCTAGCATCGATCGTAGCTAACGATCGTACGATCGATGCACTGATCGTAGCATCGATGCTAGCTAGCATCGATCGTACGATCGATGCATCGATCGTAGCATCGATCGTACTGACTGATCGTAGCTAGCATCGATCGTACGATGCTAGCATCGATCGTA";
    let mut minor = major.to_vec();
    minor[125] = if minor[125] == b'A' { b'C' } else { b'A' };
    let mut file = File::create(&reads).unwrap();
    append_reads(&mut file, "major", major, 90, 4, 4);
    append_reads(&mut file, "minor", &minor, 90, 4, 1);
    drop(file);

    let product = assemble(&config(reads, output, 31)).unwrap();

    assert!(product.stats.simple_bubbles >= 1);
    assert!(product.stats.variant_alleles >= 2);
    let min_coverage = product
        .bubble_alleles
        .iter()
        .map(|allele| allele.mean_coverage)
        .fold(f32::INFINITY, f32::min);
    let max_coverage = product
        .bubble_alleles
        .iter()
        .map(|allele| allele.mean_coverage)
        .fold(0.0_f32, f32::max);
    assert!(min_coverage < max_coverage * 0.5);
}
