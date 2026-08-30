# BridgeAsm

BridgeAsm is an experimental Rust assembler for Illumina paired-end metagenomes. Its design goal is to improve completeness and strain preservation without accepting unsupported joins merely to increase N50.

> Status: research prototype. It is not yet a replacement for MEGAHIT or metaSPAdes on production-scale data.

## Current algorithm

1. Stream FASTQ/FASTQ.GZ without retaining all reads.
2. Count exact canonical packed k-mers (`k <= 127`).
3. Retain solid k-mers and optionally rescue short weak paths between solid anchors (anchored mercy).
4. Build a bidirected de Bruijn graph using canonical nodes plus an orientation bit.
5. Compact maximal non-branching paths into unitigs.
6. Thread full reads and read pairs through oriented unitigs.
7. Emit maximal primary walks only through uniquely supported direct-read transitions.
8. Keep graph and evidence outputs for later strain-aware resolution.

The long k-mer representation is used only during construction. Persistent graph traversal uses dense integer IDs.

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
  --threads 8
```

Outputs:

- `primary_contigs.fasta`: deduplicated evidence-supported primary walks.
- `unitigs.fasta`: oriented compacted graph unitigs.
- `assembly.gfa`: unitigs and observed read/pair transitions.
- `run_profile.json`: graph, assembly, and stage timing statistics.

## Synthetic benchmark

```bash
scripts/benchmark_synthetic.sh
```

## Zymo Log benchmark

Download ERR2935805 from ENA:

```bash
scripts/download_ena_run.sh ERR2935805 data/ERR2935805
```

Run BridgeAsm and, when installed, MEGAHIT and MetaQUAST:

```bash
THREADS=16 scripts/benchmark_zymo_log.sh \
  data/ERR2935805/*_1.fastq.gz \
  data/ERR2935805/*_2.fastq.gz \
  benchmark/zymo_log \
  references/zymo.fasta
```

A fifth argument limits BridgeAsm to the first N read pairs for capacity smoke testing.

## Development principles

- Positive physical evidence may support a join; absence of evidence does not.
- Minor/low-depth paths are not deleted merely because they disappear at larger k.
- Exact contigs and N-gap scaffolds must remain separate products.
- Every optimization must preserve deterministic synthetic regression results.

## Near-term roadmap

- Memory-bounded partitioned k-mer counting and direct compacted-DBG construction.
- Unitig-native sparse seed index for parallel read/pair threading.
- Fragment-class phasing rather than global pairwise-link materialization.
- Evidence-aware bubble/repeat resolution and local assembly of hard regions.
- Protein/DNA foundation-model priors only after physical candidate generation.
