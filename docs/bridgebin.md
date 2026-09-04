# BridgeBin v0

BridgeBin is the binning arm of Bridge. The first version is deliberately small, deterministic, and dependency-light so it can serve as a measurable baseline before adding learned biological priors.

## Current model

For each contig, BridgeBin computes:

- canonical tetranucleotide frequencies (reverse-complement invariant);
- GC fraction;
- optional single- or multi-sample coverage vectors.

Long seed contigs are processed from longest to shortest. Each seed is compared against every current bin centroid and joins the best bin only when the weighted similarity exceeds a conservative threshold. Shorter eligible contigs are then rescued only if the best bin clears both a score threshold and a best-vs-second margin. This avoids the worst single-linkage chaining failure mode.

The current implementation is intentionally an exact centroid search. If the number of inferred bins approaches the number of contigs, worst-case runtime is quadratic. That makes v0 useful as an algorithmic baseline and for oracle-style experiments, but not yet the final large-coassembly implementation.

## Build

The binary lives at `src/bin/bridgebin.rs`, so Cargo discovers it automatically:

```bash
cargo build --release --bin bridgebin
```

## Input

```bash
bridgebin \
  --contigs assembly.fa \
  --coverage depth.tsv \
  --out-dir bridgebin_out
```

Coverage is a whitespace-delimited matrix. A header is recommended:

```text
contig sample_1 sample_2 sample_3
contig_1 31.2 5.4 18.1
contig_2 30.8 5.1 17.7
contig_3 7.3 25.0 3.2
```

Without a coverage table, the coverage term is automatically omitted and the remaining weights are renormalized.

## Output

- `assignments.tsv`: contig-to-bin assignments and assignment scores;
- `bins.tsv`: bin-level size/GC summary;
- `bins/bin_XXXX.fa`: one FASTA per inferred bin;
- `unbinned.fa`: filtered or ambiguous contigs unless `--no-unbinned` is used.

## Near-term roadmap

1. Add BAM/CRAM depth extraction and robust multi-sample abundance normalization.
2. Replace full centroid scan with sparse candidate generation while retaining an exact/oracle mode.
3. Add assembly-graph and paired-read linkage as positive/negative constraints.
4. Add marker-gene contamination constraints and conservative bin splitting.
5. Add ORF/protein embeddings as the Biological Brain prior.
6. Benchmark with CAMI/AMBER plus CheckM2 and GUNC against MetaBAT2, SemiBin2, VAMB-multi, COMEBin, and QuickBin.

The intended direction is not just post-assembly clustering. BridgeBin should eventually provide genome-identity probabilities back to Bridge's ambiguous graph resolution so assembly and binning can inform each other.
