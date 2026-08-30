# ERR2935805 Zymo Log subset benchmark

This report records deterministic short-read smoke tests for the experimental BridgeAsm Rust assembler. It is a development benchmark, not a claim of production superiority over MEGAHIT or metaSPAdes.

## Data and evaluation

- ENA run: `ERR2935805` (Zymo Log Illumina paired-end reads).
- Fast smoke: first 20,000 read pairs.
- Reference smoke: first 50,000 read pairs.
- Minimum reported sequence length: 200 bp.
- Reference evaluation: MetaQUAST against the ZymoBIOMICS reference collection downloaded by the workflow.
- CPU measurements are from GitHub-hosted runners and are intended for regression comparison, not hardware-independent ranking.

## Production graph policy

The current default uses:

- quality- and independent-fragment-aware solid k-mers;
- anchored mercy rescue for short weak paths spanned by reads between solid anchors;
- ordinary `(k+1)`-mer edges at the configured abundance threshold;
- once-observed edges when **both endpoints independently passed the solid-k-mer quality and fragment gates**;
- full-read evidence to rank ambiguous graph exits rather than treating a weak edge as an automatic join;
- preservation of alternate bubble alleles;
- primary bubble collapse only when multiple alternatives have bilateral physical support.

This singleton-solid-edge policy is deliberately narrower than retaining arbitrary singleton k-mers. Both endpoint nodes must already be supported by independent fragments and base-quality evidence. The edge only keeps the candidate path in the graph; primary emission still requires positive transition evidence.

## Edge-policy ablation on the 50k subset

| edge policy | k | GF | total bp | NA50 | misassemblies | mismatches/100 kb | indels/100 kb |
|---|---:|---:|---:|---:|---:|---:|---:|
| strict threshold + mercy | 21 | 1.460% | 1,065,725 | 265 | 7 | 3.95 | 1.06 |
| solid-endpoint singleton edges | 21 | **1.507%** | **1,101,806** | **268** | **6** | **3.09** | 1.10 |
| unique non-branching singleton edges | 21 | 1.508% | 1,102,282 | 268 | 7 | 3.82 | 1.10 |
| strict threshold + mercy | 31 | 1.370% | 991,461 | 242 | 0 | 2.26 | 1.05 |
| solid-endpoint singleton edges | 31 | **1.427%** | **1,013,008** | 242 | **0** | **0.91** | **0.85** |
| unique non-branching singleton edges | 31 | 1.427% | 1,014,133 | 242 | 0 | 0.91 | 0.85 |

The reference-based result favors the solid-endpoint rule. It increased genome fraction for both k values and preserved the zero-misassembly k=31 result. Restricting rescue to unique non-branching edges did not improve correctness and slightly worsened the k=21 error profile. The broad solid-endpoint rule is therefore retained, with ambiguity deferred to read-supported primary-walk selection.

## 20k read-pair continuity smoke

| assembler/configuration | sequences >=200 bp | total bp | N50 | largest |
|---|---:|---:|---:|---:|
| BridgeAsm k=21 | 577 | 159,662 | 264 | 760 |
| BridgeAsm k=31 | 375 | 93,458 | 247 | 834 |
| BridgeAsm k=51 | 35 | 8,815 | 244 | 798 |
| MEGAHIT | 515 | 211,780 | 398 | 2,076 |
| metaSPAdes | 2,550 | 791,838 | 313 | 1,800 |

Approximate peak RSS on the same hosted runner was 220-244 MB for BridgeAsm, 267 MB for MEGAHIT, and 942 MB for metaSPAdes.

## 50k read-pair MetaQUAST reference smoke

| assembler/configuration | genome fraction | total assembled bp | N50 | NA50 | misassemblies | mismatches/100 kb | indels/100 kb | duplication ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BridgeAsm k=21 | 1.507% | 1,101,806 | 326 | 268 | 6 | 3.09 | 1.10 | 1.022 |
| BridgeAsm k=31 | 1.427% | 1,013,008 | 287 | 242 | 0 | 0.91 | 0.85 | 1.022 |
| MEGAHIT | 2.284% | 686,693 | 349 | 270 | 79 | 38.06 | 6.38 | 1.098 |
| metaSPAdes | 3.169% | 1,029,504 | 337 | 302 | 47 | 49.62 | 6.91 | 1.178 |

## Interpretation

BridgeAsm currently occupies the high-precision, lower-recall end of the assembly tradeoff:

- k=31 produced zero reported misassemblies and substantially fewer base errors in this subset;
- k=21 recovered more reference sequence but introduced several misassemblies;
- MEGAHIT and metaSPAdes recovered more reference sequence and longer paths, but with many more errors on this small, low-depth subset;
- solid-endpoint weak-edge rescue produced a modest but measurable recall gain without degrading the k=31 correctness result;
- the remaining completeness gap cannot be closed by admitting weak edges alone.

The next useful algorithmic work is:

1. graph-to-primary emission that uses physically supported walks without duplicating constituent unitigs;
2. bounded local path resolution around high-coverage graph breaks;
3. unitig-native sparse indexing and memory-bounded parallel threading;
4. fragment-class phasing rather than global pairwise-link materialization;
5. reference-free candidate generation followed by calibrated read/pair evidence.

Large-k graphs and pretrained biological priors should only be revisited after a hard-region audit shows that they address a measured, recoverable failure class.

## Reproduction

```bash
# Real-read continuity smoke
# GitHub Actions: .github/workflows/real-data-smoke.yml

# Reference-based smoke
# GitHub Actions: .github/workflows/zymo-reference-smoke.yml

# Full local benchmark when sufficient RAM and disk are available
THREADS=16 scripts/benchmark_zymo_log.sh \
  data/ERR2935805/*_1.fastq.gz \
  data/ERR2935805/*_2.fastq.gz \
  benchmark/zymo_log \
  references/zymo.fasta
```
