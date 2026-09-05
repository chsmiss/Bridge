# BridgeAsm

BridgeAsm is an experimental Rust assembler for Illumina paired-end metagenomes. Its design goal is to increase recoverable sequence and strain preservation while requiring positive read or fragment evidence for ambiguous joins.

> **Status:** research prototype. It is not yet a production replacement for MEGAHIT or metaSPAdes. Current Zymo tests show substantially lower assembly error rates, but lower contiguity and slightly lower genome fraction than MEGAHIT at 500,000 read pairs.

## Current algorithm

1. Stream FASTQ/FASTQ.GZ without retaining all reads.
2. Count exact canonical packed k-mers (`k <= 147`), covering standard 150 bp Illumina reads.
3. Apply quality- and independent-fragment-aware solid-k-mer filtering.
4. Rescue short weak paths between solid anchors with anchored mercy.
5. Build an edge-centric bidirected de Bruijn graph using canonical nodes plus an orientation bit.
6. Compact maximal non-branching paths into oriented unitigs.
7. Thread complete reads and paired fragments through unitigs.
8. Preserve alternate bubble alleles; collapse a primary bubble path only when alternatives have bilateral physical support.
9. Emit evidence-supported primary walks and optional paired-end N-gap scaffolds as separate products.

The `k=147` limit is intended for 150 bp reads. A dataset can only use k values no longer than its actual reads; ERR2935805 contains 101 bp mates, so its largest meaningful odd k is 99. Long k-mer keys are construction-time objects. Persistent graph traversal uses dense integer IDs.

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
- real ERR2935805 Zymo Log subsets against MEGAHIT and metaSPAdes;
- MetaQUAST evaluation against the Zymo reference collection;
- a 500,000-pair cross-k recovery benchmark;
- a focused medium-k, mercy, and primary-emission sweep intended to close the remaining GF gap.

See [`docs/benchmark_zymo_subset.md`](docs/benchmark_zymo_subset.md) and [`docs/benchmarks/zymo_log_500k_recovery.md`](docs/benchmarks/zymo_log_500k_recovery.md) for measured results and negative ablations.

### Current 500k result

On the deterministic first 500,000 ERR2935805 read pairs:

- best BridgeAsm single-k primary GF: **4.740%** at k=31;
- exact cross-k recovery-union GF: **5.095%**;
- MEGAHIT GF: **5.461%**;
- BridgeAsm retained substantially fewer mismatches, indels, and misassemblies;
- the remaining unique-coverage gap is about **0.366 percentage points**, concentrated mainly in the Pseudomonas reference.

The cross-k union is a diagnostic recovery catalog, not yet a production primary assembly, because overlapping paths from several k values produce a high duplication ratio.

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
- Weak-edge rescue must improve reference recovery without an unacceptable increase in switches or misassemblies; otherwise it is rejected.
- Exact contigs and N-gap scaffolds remain separate products.
- Cross-k recovery outputs must be clustered and deduplicated without reference truth before promotion to a production assembly.
- Every optimization must preserve deterministic synthetic regression results.

## Near-term roadmap

1. Close the remaining 500k GF gap using medium-k complementarity, calibrated mercy, and evidence-supported emission.
2. Replace exact-only cross-k deduplication with reference-free containment and overlap clustering.
3. Add memory-bounded partitioned counting and direct compacted-DBG construction.
4. Add unitig-native sparse indexing for parallel read/pair threading.
5. Replace global pairwise links with fragment-class phasing.
6. Add bounded local path resolution for high-coverage graph breaks.
7. Use repeat/copy-number flow and biological model priors only after measured failure attribution.
