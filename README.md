# BridgeAsm

BridgeAsm is an experimental Rust assembler for Illumina paired-end metagenomes. Its design goal is to increase recoverable sequence and strain preservation while requiring positive read or fragment evidence for ambiguous joins.

> **Status:** research prototype. It is not yet a production replacement for MEGAHIT or metaSPAdes. Current Zymo subset tests show substantially lower assembly error rates, but lower genome fraction and contiguity.

## Current algorithm

1. Stream FASTQ/FASTQ.GZ without retaining all reads.
2. Count exact canonical packed k-mers (`k <= 127`).
3. Apply quality- and independent-fragment-aware solid-k-mer filtering.
4. Rescue short weak paths between solid anchors with anchored mercy.
5. Build an edge-centric bidirected de Bruijn graph using canonical nodes plus an orientation bit.
6. Compact maximal non-branching paths into oriented unitigs.
7. Thread complete reads and paired fragments through unitigs.
8. Preserve alternate bubble alleles; collapse a primary bubble path only when alternatives have bilateral physical support.
9. Emit evidence-supported primary walks and optional paired-end N-gap scaffolds as separate products.

Long k-mer keys are construction-time objects. Persistent graph traversal uses dense integer IDs.

## Build

```bash
cargo build --release
```

## Run

```bash
target/release/bridgeasm assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o result \
  -k 31 \
  --min-count 2 \
  --mercy-max-kmers 16 \
  --min-read-support 2 \
  --min-pair-support 2 \
  --threads 8
```

Main outputs:

- `primary_contigs.fasta`: deduplicated exact primary walks.
- `primary_scaffolds.fasta`: pair-supported scaffolds with explicit `N` gaps.
- `primary_scaffolds.agp`: scaffold component and gap provenance.
- `unitigs.fasta`: compacted graph unitigs.
- `variants.fasta` and `haplotigs.fasta`: preserved alternate paths.
- `assembly.gfa`: unitigs and evidence tags.
- `run_profile.json`: graph, evidence, assembly, and timing statistics.

## Tests and benchmarks

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
scripts/benchmark_synthetic.sh
```

The repository includes deterministic GitHub Actions workflows for:

- synthetic major/minor strain regression;
- a real ERR2935805 Zymo Log subset against MEGAHIT and metaSPAdes;
- MetaQUAST evaluation against the Zymo reference collection.

See [`docs/benchmark_zymo_subset.md`](docs/benchmark_zymo_subset.md) for measured results and negative ablations.

### Full Zymo Log benchmark

```bash
scripts/download_ena_run.sh ERR2935805 data/ERR2935805

THREADS=16 scripts/benchmark_zymo_log.sh \
  data/ERR2935805/*_1.fastq.gz \
  data/ERR2935805/*_2.fastq.gz \
  benchmark/zymo_log \
  references/zymo.fasta
```

A fifth argument limits the run to the first N read pairs for capacity testing.

## Development principles

- Positive physical evidence may support a join; absence of evidence does not.
- Minor and low-depth paths are not deleted merely because coverage is lower or a path disappears at larger k.
- Weak-edge rescue must improve reference recovery without increasing switches or misassemblies; otherwise it is rejected.
- Exact contigs and N-gap scaffolds remain separate products.
- Every optimization must preserve deterministic synthetic regression results.

## Near-term roadmap

1. Memory-bounded partitioned counting and direct compacted-DBG construction.
2. Unitig-native sparse indexing for parallel read/pair threading.
3. Fragment-class phasing instead of global pairwise-link materialization.
4. Bounded local path resolution for high-coverage graph breaks.
5. Repeat/copy-number flow only after measured hard-region attribution.
6. Protein/DNA model priors only for candidate reranking after physical candidate generation.
