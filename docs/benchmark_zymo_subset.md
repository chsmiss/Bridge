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
- strict observed `(k+1)`-mer edge abundance thresholds;
- anchored mercy rescue for short weak paths spanned by reads between solid anchors;
- preservation of alternate bubble alleles;
- primary-walk collapse only when multiple alternatives are physically supported on both flanks;
- pair-supported N-gap scaffolds kept separate from exact contigs.

Two broader weak-edge rules were tested and rejected:

1. Retain every once-observed edge whose endpoints are both solid.
2. Retain only once-observed solid-endpoint edges that are unique missing continuations.

The second rule changed k=21 genome fraction from 1.507% to 1.508% on the 50k subset, while misassemblies increased from 6 to 7 and mismatches from 3.09 to 3.82 per 100 kb. It produced no substantive k=31 benefit. The production graph therefore remains strict edge threshold plus anchored mercy.

## 20k read-pair continuity smoke

| assembler/configuration | sequences >=200 bp | total bp | N50 | largest |
|---|---:|---:|---:|---:|
| BridgeAsm k=21 | 689 | 182,206 | 259 | 681 |
| BridgeAsm k=31 | 333 | 85,794 | 250 | 834 |
| BridgeAsm k=51 | 32 | 8,041 | 234 | 798 |
| MEGAHIT | 515 | 211,780 | 398 | 2,076 |
| metaSPAdes | 2,550 | 791,838 | 313 | 1,800 |

Approximate peak RSS on the same hosted runner was 220-244 MB for BridgeAsm, 267 MB for MEGAHIT, and 942 MB for metaSPAdes.

## 50k read-pair MetaQUAST reference smoke

The table below records the completed reference smoke used to evaluate the rejected weak-edge experiment. The strict production run is regenerated automatically by `zymo-reference-smoke.yml`; the weak-edge comparison is included here because it establishes the rollback decision.

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
- simply admitting weak singleton edges did not close the completeness gap safely.

The next useful algorithmic work is not a broader weak-edge rule. Higher-value directions are:

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
