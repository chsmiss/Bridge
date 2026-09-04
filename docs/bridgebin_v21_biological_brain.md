# BridgeBin v2.1 — Biological Brain

BridgeBin v2.1 separates two problems that were conflated in the original binner:

1. **candidate generation / partition mechanics** — can the algorithm expose and act on the right ambiguous relationships without catastrophic transitive merges?
2. **genome-identity inference** — can biological representations tell whether two ambiguous contigs are from the same genome?

The production architecture is therefore:

```text
assembly contigs
  -> conservative signed v2 core
  -> sparse truth-free candidate mining
  -> Biological Brain feature extraction
  -> calibrated P(same genome | pair evidence)
  -> split contaminated bins
  -> constrained cross-bin merge
  -> residual rescue
  -> final bins
```

Truth is never used in this production path.

## Why v1 was replaced

The v1 family effectively reduced heterogeneous evidence to one scalar score followed by greedy/reconciliation merges. That is structurally dangerous: enough positive coverage/composition evidence can numerically cancel a decisive contradiction such as a duplicated single-copy marker.

v2/v2.1 keeps positive and negative evidence distinct. A trusted cannot-link remains a veto after transitive unions.

Examples of negative/veto evidence:

- duplicated trusted single-copy marker,
- calibrated pair model with very low `p_same`,
- explicit assembly/path mutually-exclusive constraint,
- high-confidence incompatible taxonomy when available.

Examples of positive evidence:

- multi-sample abundance agreement,
- nucleotide composition,
- physical assembly/read/long-read links,
- DNA foundation-model identity,
- coding architecture,
- gene/protein repertoire,
- calibrated pair-model support.

A positive edge never deletes a hard negative.

## v2.1 pair-score interface

The Rust refinement layer consumes a small stable TSV schema:

```text
left    right    p_same    confidence    model
```

`p_same` is the calibrated probability/decision score that two contigs originate from the same genome. Model inference remains outside Rust so DNA language models, GENERanno, ESM-C, or future models can be swapped without rewriting the partitioner.

Current v2.1 defaults:

```text
split_max_same     0.12
join_min_same      0.88
rescue_min_same    0.84
rescue_margin      0.08
min_pair_support   2
cross-bin merge    p_same >= 0.92, >=3 supports
```

Validation-derived pair-model thresholds can override these defaults.

## Three refinement stages

### 1. Split

Inside a current bin, a calibrated low-`p_same` pair or duplicated SCG becomes a hard cannot-link. High-`p_same` edges can connect seeds only while respecting every hard conflict.

Small ambiguous fragments are not forced into a component merely to maximize assigned bp.

### 2. Constrained cross-bin merge

Over-split pure v2 components may be reconnected only when several independent high-`p_same` cross-bin links agree.

Any confident low-`p_same` pair or shared SCG blocks the bin merge. The veto is component-level and therefore survives transitive merges.

### 3. Residual rescue

An unbinned contig is compared with multiple anchors in candidate bins. Assignment requires:

- enough independent pair supports,
- high aggregated posterior,
- a best-vs-second margin,
- no marker/hard-negative contradiction.

Otherwise it remains unbinned.

## Sparse Biological Brain candidate mining

Running an expensive foundation model on every contig pair is unnecessary and scales as `O(N^2)`.

`bridgebin_candidate_pairs.py` mines only:

- `within_bin_anchor` / `within_bin_neighbor` / `within_bin_contrast` pairs to find hidden mixtures,
- `cross_bin_merge` pairs between nearby current bins,
- `residual_rescue` pairs from unbinned contigs to candidate-bin anchors.

Candidate selection uses cheap coverage, canonical 4-mer composition, and GC only for **recall**. Cheap similarity never forces a final merge.

### Diverse anchors

Using only the longest contigs is unsafe. In the Zymo v2 fixture, all ten longest anchors of one mixed E. coli/Salmonella bin came from E. coli, hiding the Salmonella subpopulation from the Biological Brain.

The current miner starts with the longest contig, then performs a cheap-feature farthest-point pass over a long-contig pool. On the deterministic Zymo reference-fragment fixture this increased benchmark-only candidate exposure from:

```text
mixed-bin cross-genome exposure       100%
fragmented-genome cross-bin exposure   85.7%
residual same-genome exposure         100%
```

to:

```text
mixed-bin cross-genome exposure       100%
fragmented-genome cross-bin exposure  100%
residual same-genome exposure         100%
```

with about 43.7k candidate pairs, well below the default 250k budget. These numbers are truth-only diagnostics; the miner itself does not see truth.

## Biological Brain modalities

### DNA identity

Preferred role: organism/genome identity, especially when coverage is missing or shared by close genomes.

Adapters currently support the DNABERT-S authors' public checkpoint (`zhihan1996/DNABERT-S`) and generic compatible DNA models. Long contigs are windowed and encoded in both orientations before aggregation.

A cheap canonical 5-mer embedding is always evaluated as a baseline: a DNA foundation model should not become a default dependency unless it adds hard-negative separation beyond k-mer composition.

### GENERanno coding architecture

The public prokaryotic GENERanno CDS annotator is treated as a strand-aware nucleotide CDS classifier, not as a gene-function annotator.

The adapter reconstructs supported ORFs and emits:

- CDS coordinates and strand,
- translated proteins,
- coding density,
- ORF density and length distribution,
- strand balance,
- contig-level coding-architecture vector.

Functional repertoire is intentionally delegated to protein-family/ESM-C processing rather than invented from CDS labels.

### ESM-C

A naive mean over all protein embeddings can be dominated by conserved housekeeping proteins and make related genomes look artificially similar.

BridgeBin therefore keeps two possible protein views:

1. pooled mean/std representation as a generic signal,
2. **protein repertoire**: an unsupervised spherical-k-means prototype codebook followed by a contig-level TF-IDF histogram.

The repertoire view upweights rare/accessory protein content relative to ubiquitous protein families.

For production, a fixed broad codebook is preferable to fitting one independently in every sample.

## Pair head

`bridgebin_pair_head.py` currently learns a calibrated same-genome classifier from similarities and explicit modality-presence features.

Separate modalities:

```text
dna
gene
architecture
protein
repertoire
taxonomy
coverage
composition
gc
physical support
```

Vector similarities use:

```text
max(0, cosine)
```

not `(cos + 1) / 2`; orthogonal representations are zero positive evidence rather than 0.5 similarity.

False merges receive a larger training weight by default.

## Leakage-safe validation

Randomly splitting pair rows is invalid because contigs from the same genome can appear in both train and validation sets.

If `left_genome/right_genome` metadata is present, the pair head uses **whole-genome hold-out**:

- choose held-out genomes,
- remove every pair touching those genomes from training,
- evaluate only on pairs involving held-out genomes.

The model stores:

- validation protocol,
- high-precision recommended join/split thresholds,
- false-merge rate,
- learned modality weights.

Training datasets should deliberately oversample hard negatives such as different genomes with nearly identical coverage profiles.

## Zymo structural upper-bound experiment

The deterministic Zymo reference-fragment benchmark deliberately gives E. coli and Salmonella identical six-sample coverage profiles.

Baseline signed v2-balanced:

```text
bins             15
weighted purity  0.8455
bp recall        0.8341
F1               0.8398
ARI              0.5769
HQ-like genomes  0
binned bp        0.8624
```

A **truth-derived oracle pair table** was used only to test the capability of the downstream refinement algorithm. It is not a biological model and must never be reported as such.

Oracle pair evidence, split + rescue, cross-bin merge disabled:

```text
bins             16
weighted purity  0.99385
bp recall        0.97166
F1               0.98263
ARI              0.96659
HQ-like genomes  8 / 8
binned bp        1.00000
```

Same oracle pair evidence with constrained cross-bin merge:

```text
bins             10
weighted purity  0.99385
bp recall        0.98376
F1               0.98878
ARI              0.98116
HQ-like genomes  8 / 8
binned bp        1.00000
```

The merge stage accepted six candidate bin merges while 86 bin pairs were blocked by hard-negative evidence.

**Interpretation:** the current partition mechanics are capable of recovering the intended genomes when pair-level genome identity is correct. The principal research bottleneck has moved to learning and calibrating real pair evidence.

## What is not yet proven

The following must remain separate from the oracle upper bound:

- real DNABERT-S gain over 5-mer composition,
- real ESM-C/protein-repertoire gain,
- real GENERanno coding-architecture gain,
- performance on actual assembler contigs rather than reference fragments,
- strain-rich natural communities,
- HiFi/ONT-specific weighting,
- generalization to unseen taxa and unseen abundance regimes.

## Required next benchmarks

A production candidate should be accepted only after the following ladder:

1. **representation test**
   - frozen DNA/protein feature AUC and nearest-neighbor genome accuracy,
   - compare directly with canonical 5-mer,
   - explicitly report close-genome pairs such as E. coli/Salmonella.
2. **pair-model test**
   - whole-genome hold-out,
   - false-merge rate,
   - recall at >=99% and >=99.5% precision,
   - calibration curve.
3. **candidate-recall test**
   - mixed-bin exposure,
   - fragmented-genome cross-bin exposure,
   - residual exposure,
   - candidate budget.
4. **real contig MAG benchmark**
   - Bridge, MEGAHIT, metaSPAdes contigs,
   - compare BridgeBin with SemiBin2, COMEBin, VAMB/TaxVAMB, MetaBAT2 and relevant long-read tools,
   - CheckM2/CheckM/AMBER-style completeness/contamination plus ARI/purity when truth exists.
5. **technology regimes**
   - Illumina single sample,
   - Illumina multi-sample,
   - PacBio HiFi,
   - ONT,
   - hybrid/co-assembly.

## Production entry point

`bridgebin_biobrain_pipeline.py` orchestrates the no-truth path and writes an auditable manifest. Expensive model inference may be run inline or supplied as precomputed TSVs.

The oracle fixture is intentionally not accepted by that production entry point.
