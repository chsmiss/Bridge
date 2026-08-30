use bridgeasm::assembler::{assemble, AssembleConfig};
use std::fs::File;
use std::io::Write;

fn write_reads(path: &std::path::Path, genome: &[u8], read_length: usize, step: usize) {
    let mut file = File::create(path).unwrap();
    for (index, start) in (0..=genome.len() - read_length).step_by(step).enumerate() {
        let sequence = &genome[start..start + read_length];
        writeln!(file, "@r{index}").unwrap();
        writeln!(file, "{}", String::from_utf8_lossy(sequence)).unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "{}", "I".repeat(sequence.len())).unwrap();
    }
}

#[test]
fn assembles_a_linear_genome() {
    let directory = tempfile::tempdir().unwrap();
    let reads = directory.path().join("reads.fastq");
    let output = directory.path().join("out");
    let genome = b"ACGTTGCAAGTCGATCGTACCTGACTGATCGTAGCTAGCTACGATCGATGCTAGCATCGATCGTACGATGCTAGCTAGCATGCTAGCATCGATCGTAGCTA";
    write_reads(&reads, genome, 60, 5);

    let product = assemble(&AssembleConfig {
        read1: reads,
        read2: None,
        output_dir: output,
        k: 21,
        min_count: 1,
        mercy_max_kmers: 0,
        mercy_min_support: 1,
        min_read_support: 1,
        min_pair_support: 1,
        min_contig_length: 20,
        max_pairs: None,
        threads: 1,
    })
    .unwrap();

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

    let product = assemble(&AssembleConfig {
        read1: reads,
        read2: None,
        output_dir: output,
        k: 41,
        min_count: 1,
        mercy_max_kmers: 0,
        mercy_min_support: 1,
        min_read_support: 1,
        min_pair_support: 1,
        min_contig_length: 20,
        max_pairs: None,
        threads: 1,
    })
    .unwrap();

    assert!(product.stats.primary_contigs >= 1);
    assert!(product.stats.graph.canonical_nodes > 0);
}
