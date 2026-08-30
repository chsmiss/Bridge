# Zymo Log 50k-pair mercy baseline

## Provenance

- Dataset: ERR2935805, deterministic first 50,000 paired reads.
- Reference: ZymoBIOMICS official reference genomes, Zenodo record 3935737, `ZymoBIOMICS.STD.refseq.v2.zip`.
- BridgeAsm commit: `0f081d623afc2320f258d2eb659685e99e0924dd`.
- Threads: 2.
- Minimum reported contig: 200 bp.
- Evaluator: MetaQUAST 5.3.0 against the combined official reference set.

This run predates commit `8c3a54821aff5a94a966010c1bc10c783e817618`, which changed mercy support to count independent physical fragments and prohibited rescue across discontinuous read intervals. It is retained as a frozen baseline and must not be mixed with later fragment-aware mercy results.

## Results

| Assembly | Total bp | N50 | GF (%) | Misassemblies | Mismatches / 100 kbp | Indels / 100 kbp | Wall | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BridgeAsm default, mercy 16/support 1 | 973,686 | 318 | 1.370 | 0 | 2.26 | 0.62 | 5.99 s | 623 MiB |
| BridgeAsm long mercy 64/support 1 | 1,263,970 | 339 | 1.735 | 6 | 7.77 | 1.11 | 6.56 s | 849 MiB |
| BridgeAsm long mercy 64/support 2 | 558,287 | 277 | 0.796 | 0 | 2.69 | 1.08 | 6.18 s | 615 MiB |
| MEGAHIT | 1,650,352 | 623 | 2.284 | 79 | 36.26 | 3.17 | 29.04 s | 296 MiB |
| metaSPAdes | 2,305,711 | 579 | 3.169 | 47 | 50.34 | 5.40 | 24.44 s | 245 MiB |

## Graph diagnostics

| BridgeAsm mode | Retained solid | Rescued k-mers | Unitigs | Unitig edges | Branching unitigs | Unitig N50 | Direct transitions | Primary transitions |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Default | 1,695,488 | 58,240 | 35,312 | 4,044 | 1,725 | 168 | 2,649 | 64 |
| Long mercy support 1 | 1,695,488 | 284,766 | 39,924 | 23,660 | 10,518 | 188 | 14,786 | 285 |
| Long mercy support 2 | 1,695,488 | 3,542 | 48,503 | 2,734 | 1,149 | 114 | 1,869 | 48 |

## Interpretation

1. The default mode occupies the high-precision end of the Pareto frontier: zero reported misassemblies and low base-error rates, but substantially lower genome fraction and continuity.
2. Long mercy with one supporting occurrence recovers real reference sequence: genome fraction increases by about 27% relative to default. It also introduces six misassemblies and triples the mismatch rate, so it is not suitable as the production default.
3. The old support-2 implementation retained very little sequence and is superseded by fragment-aware mercy semantics. Its low recall must not be interpreted as evidence against independent-fragment rescue.
4. MEGAHIT and metaSPAdes recover more reference sequence on this sparse subset, but at much higher reported structural and base-error counts.
5. Long mercy creates both useful connectivity and many additional branches. Future rescue must separate high-quality, independently supported weak paths from single-fragment noise rather than simply increasing the maximum mercy span.

## Next experiment

Run the same reference benchmark after fragment-aware mercy support with long-mercy quality thresholds Q25, Q30, and Q35. Promote a setting only if it increases genome fraction while retaining a clearly better correctness profile than the comparison assemblers.
