# Zymo Log 500k-pair cross-k recovery benchmark

## Provenance

- Dataset: ERR2935805, deterministic first 500,000 read pairs.
- Observed read length: 101 bp per mate.
- References: ZymoBIOMICS reference collection, Zenodo record 3935737.
- BridgeAsm workflow commit: `26a5a8edb1f99489353d5c9a35cea1f73012e190`.
- GitHub Actions run: `33288385719`.
- Evaluator: MetaQUAST 5.3.0, minimum reported sequence length 200 bp.
- Threads: 2 for BridgeAsm and MEGAHIT.

This is a development subset, not the complete 47.8-million-pair ERR2935805 run.

## Primary and cross-k recovery results

| Assembly | GF (%) | N50 | NA50 | Largest | Misassemblies | Mismatches / 100 kbp | Indels / 100 kbp | Duplication |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BridgeAsm k21 primary | 4.687 | 1,188 | 1,188 | 7,274 | 10 | 2.04 | 0.85 | 1.860 |
| BridgeAsm k31 primary | 4.740 | 1,997 | 1,997 | 12,441 | 4 | 2.22 | 0.75 | 1.880 |
| BridgeAsm k51 primary | 4.183 | 1,543 | 1,543 | 8,759 | 0 | 1.77 | 0.43 | 2.001 |
| BridgeAsm k91 primary | 1.017 | 239 | 239 | 2,146 | 0 | 1.47 | 0.93 | 1.098 |
| Cross-k primary exact union | 5.086 | 1,483 | 1,483 | 12,441 | 14 | 1.97 | 0.66 | 5.172 |
| Cross-k all-exact union | 5.095 | 1,265 | 1,265 | 12,441 | 18 | 1.97 | 0.65 | 6.417 |
| MEGAHIT | **5.461** | **269,356** | **237,827** | **545,716** | 128 | 27.31 | 4.18 | 1.036 |

## Resource snapshot

| Assembly | Wall | Peak RSS |
|---|---:|---:|
| BridgeAsm k21 | ~75 s | ~2.78 GiB |
| BridgeAsm k31 | ~64 s | ~2.46 GiB |
| BridgeAsm k51 | ~51 s | ~2.37 GiB |
| BridgeAsm k91 | ~21 s | ~1.27 GiB |
| MEGAHIT | ~115 s | ~0.33 GiB |

## Interpretation

1. Increasing from 50k to 500k pairs makes the completeness comparison much more meaningful. BridgeAsm k31 reaches 4.740% GF versus MEGAHIT at 5.461%.
2. Exact cross-k union reaches 5.095%, leaving a 0.366-percentage-point gap to MEGAHIT while retaining much lower base-error and structural-error counts.
3. The cross-k union is a diagnostic recovery catalog, not a production primary assembly. Its duplication ratio is too high because overlapping paths from several k values are retained independently.
4. k21 and k31 are genuinely complementary. k51 contributes little unique sequence and k91 is coverage-limited at this sampling depth.
5. The 101-bp Zymo Log reads cannot produce k147 windows. k147 support is retained for 150-bp datasets, while this dataset's largest meaningful odd k is 99.
6. The remaining GF gap is concentrated mainly in Pseudomonas aeruginosa. BridgeAsm's all-exact union covers about 0.69 Mb of that reference versus about 0.99 Mb for MEGAHIT; Listeria monocytogenes recovery is already comparable.
7. BridgeAsm emits more aligned sequence than MEGAHIT but less unique reference coverage. The immediate problem is therefore recovery and nonredundant path selection, not raw FASTA volume.

## Follow-up experiments

- Focused long-mercy and emission-threshold sweep at k21, k25, k29, k31 and k41.
- Separate pair-supported N-gap scaffold evaluation; scaffold GF must not be reported as gap-free contig GF.
- If a recovery union exceeds MEGAHIT GF, add reference-free containment/overlap clustering before considering it a production output.
- Validate any selected setting on independent mock communities before claiming general superiority.
