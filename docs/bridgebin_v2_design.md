# BridgeBin v2: signed multimodal genome binning

Status: experimental redesign on `codex/bridgebin-v2-signed-evidence`.

## Why v1 is being replaced

BridgeBin v1 used a useful but incomplete idea: start with conservative bins, then greedily reconcile reciprocal-best bin pairs. On the Zymo reference benchmark this increased recall/completeness but created catastrophic cross-genome merges (for example Listeria + Enterococcus), reducing purity, ARI, F1, and HQ-like MAG count.

The failure mode is structural, not just a threshold problem. A single scalar similarity lets strong positive evidence compensate for biologically decisive negative evidence. Transitive greedy merging can then make a locally plausible error irreversible.

BridgeBin v2 changes the problem definition:

> Binning is a sparse signed-graph partitioning problem with heterogeneous, confidence-weighted evidence and explicit conflicts.

Positive and negative evidence are retained separately until the partition decision.

## Lessons imported from current binners

### MetaBAT2 / QuickBin: cheap rejection before expensive scoring

Coverage and composition remain strong first-line filters. QuickBin's staged GC/depth -> low-order k-mer -> higher-order k-mer -> small neural network cascade is particularly attractive as a systems design: most impossible pairs should be rejected before an expensive model runs. Read/assembly links should relax or strengthen a candidate, not force a merge.

### SemiBin2: representation learning is not the final clustering algorithm

SemiBin2 learns an embedding from self-supervised must-link/cannot-link examples, builds a neighborhood graph, clusters the graph, then uses single-copy genes to recognize and split mixed bins. The important lesson is that learned embeddings and biological quality control are separate stages.

### COMEBin: multi-view contrastive learning plus partition ensembles

COMEBin creates multiple sequence views, trains with contrastive learning, constructs a nearest-neighbor graph, runs many Leiden parameter settings, and selects biologically plausible bins. This motivates multi-view augmentation and multiple candidate partitions instead of trusting one global radius.

### VAMB / TaxVAMB: latent spaces and optional modalities

VAMB uses a VAE over TNF and abundance with adaptive clustering. TaxVAMB demonstrates that taxonomic information can materially improve purity and recovery when taxonomy is reliable, but can degrade performance in environments with noisy/incomplete taxonomy. BridgeBin must therefore treat taxonomy as confidence-gated optional evidence, never as an unconditional truth source.

### MetaDecoder / MetaCAT: unknown-K probabilistic residual clustering

MetaDecoder and the newer MetaCAT avoid assuming a fixed number of genomes. MetaCAT combines sparse affinity graphs, marker-seeded label propagation, and a sparse weighted Dirichlet-process Gaussian mixture for large/complex residual clusters. BridgeBin should eventually use probabilistic responsibilities for residual/ambiguous contigs rather than force every contig through one hard threshold.

### Binny / marker-driven refiners: biological quality belongs inside the loop

Single-copy markers, lineage-specific expectations, rRNA/mobile-region handling, and iterative split/recluster are effective because they detect failure modes that coverage/composition alone cannot see. Marker evidence should be lineage-aware and copy-number-aware rather than a generic duplicate-count heuristic.

### Long-read binners (SemiBin2 long-read, LorBin, LRBinner)

Long-read contigs have stronger sequence/gene features but variable cluster density and can contain chimeric joins. Strong methods generate multiple DBSCAN/BIRCH-like candidate clusters and select by biological quality. BridgeBin should use the same evidence model for short and long reads but different candidate-generation and split policies.

## Cross-domain formulation

BridgeBin v2 borrows from constrained clustering, correlation clustering/minimum-cost multicut, entity resolution, and data association.

For candidate contig pair `(i,j)`, retain two quantities:

- `A_ij`: attraction / evidence that the pair is from the same genome.
- `R_ij`: repulsion / evidence that the pair cannot or should not share a genome.

This is intentionally different from a single similarity score.

A future global objective is a signed correlation-clustering / multicut problem:

```
min partition P
    sum positive-edge-costs cut by P
  + sum negative-edge-costs kept inside P
  + bin-quality penalties(P)
```

The number of bins does not need to be specified in advance.

The first Rust implementation uses a simpler conservative constrained-agglomeration approximation so we can validate the evidence model before introducing a heavy optimizer.

## Evidence model

### 1. Intrinsic sequence evidence

Per contig/window:

- canonical 3/4/5-mer frequencies (v2 prototype starts with 5-mer)
- GC and GC-compensated composition
- entropy / low-complexity masks
- later: codon usage, amino-acid composition
- later: mask or separately model rRNA, CRISPR, transposons and other mobile/repetitive sequence before computing genome-signature statistics

### 2. Abundance evidence

Per sample:

- mean depth and, when available, depth variance
- log-depth ratio
- cross-sample correlation/covariance
- effective number of informative samples

Weights should be adaptive: one sample mainly constrains absolute abundance; many independent samples provide a strong differential-coverage genome fingerprint.

### 3. Physical linkage evidence

The `--links` interface carries pair evidence from:

- assembly graph adjacency
- paired-end links
- read threading
- long reads spanning multiple contigs
- Hi-C links
- future optical/link-read technologies

A positive link increases attraction but does not force a merge. A high-confidence mutually exclusive phase/path can be emitted as a `cannot_link` and remains valid across transitive component merges.

### 4. Marker / taxonomy evidence

- duplicated trustworthy single-copy markers are strong negative evidence
- future: use lineage-specific marker copy-number expectations rather than one universal marker set
- high-confidence incompatible taxonomy is a cannot-link
- taxonomy that is low-confidence or missing contributes nothing
- compatible taxonomy may add only a modest positive prior

### 5. Gene and protein-model evidence

BridgeBin should not run large foundation models inside the Rust binary. Python/GPU adapters produce per-contig features, and Rust consumes those features through `--bio-features`.

Supported v2 prototype fields:

```
contig	taxonomy	taxonomy_confidence	gene_profile	gene_confidence	esm_embedding	protein_confidence
```

Vector fields are comma-separated. Missing modalities are allowed.

Recommended biological feature pipeline:

1. call CDS/ORFs with Prodigal/FragGeneScan or GENERanno;
2. retain gene coordinates, coding density, strand/order and annotation confidence;
3. encode translated ORFs with ESM-C;
4. aggregate proteins as a set, not only a naive global mean;
5. emit gene-family/profile features and compact protein embedding features;
6. optionally emit taxonomy posteriors and model confidence.

### Why ESM-C should be an adjudicator, not the primary binner

Housekeeping proteins are shared across many genomes and close taxa can have similar protein-language-model embeddings. A single mean ESM embedding is therefore not a safe genome identifier. It is more useful for hard candidate edges after coverage/composition have already narrowed the possibilities.

For a contig with ORFs `p_1 ... p_m`, future aggregation should retain distributional structure, for example:

```
protein_signature = [mean(z_p), variance(z_p), prototype_histogram(z_p), marker_specific(z_p)]
```

A small learned projection/head can then estimate pair compatibility. The ESM-C base can remain frozen initially. This also reuses the existing Biological Brain model infrastructure without coupling checkpoint code to the Rust assembler/binner.

### Why GENERanno is useful

GENERanno is attractive for two reasons:

- direct prokaryotic/metagenomic CDS annotation from nucleotide sequence;
- contextual DNA representations that can supplement short-contig composition and gene-boundary evidence.

The first experiment should use its CDS/gene outputs and confidence rather than assume its hidden representation is automatically a genome-identity embedding.

## v2 partition algorithm (prototype)

### Step A: build sparse candidates

For each eligible contig, use cheap GC + coverage screening and keep only the best `K` neighbors. Add every explicit physical-link pair even if it falls outside the neighborhood search.

The current implementation still scans candidate distances exactly before retaining a bounded neighborhood. This is deliberate for correctness on the Zymo development benchmark. Production scaling should replace it with GC/depth bins or HNSW/ANN.

### Step B: score signed edges

For each candidate pair compute:

```
PairEvidence {
    attraction,
    repulsion,
    composition,
    coverage,
    gc,
    gene,
    protein,
    marker_conflict,
    taxonomy_conflict,
    external_cannot_link,
}
```

Positive physical links boost attraction multiplicatively. Strong conflicts cannot be averaged away.

### Step C: build high-purity core components

Sort candidate edges by `attraction - lambda * repulsion` and merge only when:

- edge attraction is high;
- edge repulsion is low;
- no hard marker/taxonomy/external conflict exists;
- the two entire components remain compatible by centroid composition/coverage/GC;
- no hard cannot-link exists between any member of the two components.

The final condition prevents an A-B, B-C chain from silently violating a known A-C conflict.

### Step D: conservative residual rescue

Small components are not automatically emitted as bins. They are attached to an established core only when:

- no component-level conflict exists;
- best attraction exceeds a rescue threshold;
- best-vs-second margin is large enough.

Otherwise they remain unbinned. This is purity-first by design.

### Step E (planned): split mixed cores

Before any later merge/reconciliation, a core should be split if it shows:

- duplicated lineage-aware SCGs;
- multiple high-confidence taxonomic modes;
- multimodal coverage;
- incompatible graph phases;
- bimodal gene/protein signatures.

Candidate split engines, in priority order for experiments:

1. weighted seeded KMeans for small marker-defined mixtures;
2. variable-density DBSCAN/HDBSCAN for long-read clusters;
3. sparse weighted DPGMM for large residual clusters;
4. global signed multicut once pair calibration is reliable.

## Short-read, long-read and hybrid modes

### Illumina / short read

Prioritize:

1. differential coverage
2. paired-end/read-thread/assembly-graph linkage
3. 4/5-mer composition on sufficiently long contigs
4. marker/gene evidence
5. ESM-C/GeneAnno adjudication for hard edges

### HiFi / ONT long read

Prioritize:

1. long-read physical span and assembly topology
2. longer-window composition/codon/gene evidence
3. variable-density candidate clustering
4. window-level chimera detection before whole-contig assignment
5. error-aware ORF confidence for noisy ONT

### Hybrid

Use all modalities with confidence gating and modality dropout during learned representation training.

## Training strategy for Biological Brain features

Do not start by training a model to predict bins. Start with pairwise edge calibration.

Create truth-labelled candidate pairs from Zymo/CAMI references:

- positive: same reference genome
- hard negative: different genomes with similar coverage/composition
- very hard negative: same genus/species-complex/strain-like pairs

Train/evaluate:

```
P(same_genome | coverage, composition, graph, markers, genes, ESM-C, taxonomy)
```

Primary metric is not global AUROC. It is false-positive merge rate at very high precision, stratified by difficulty.

Only after pair calibration works should the score be used inside graph partitioning.

## Required ablations

Every new modality must justify itself:

- coverage + composition baseline
- + assembly/read links
- + marker conflict
- + taxonomy
- + gene profile / GENERanno
- + ESM-C
- + learned multimodal projection
- + global optimizer

Measure the marginal gain and the new failure modes at each step.

## Benchmark plan

Fast development:

- current deterministic Zymo reference-fragment benchmark
- current real Bridge assembly outputs

Generalization:

- CAMI2 strain-madness (essential for near-neighbor genomes)
- complex simulated mixtures across abundance/depth regimes
- real human gut
- real soil/ocean where taxonomy is less complete
- PacBio HiFi and ONT long-read benchmarks

Compare at minimum:

- MetaBAT2
- SemiBin2
- COMEBin
- VAMB / TaxVAMB where input mode is appropriate
- MetaCAT
- QuickBin

Metrics:

- weighted purity / contamination
- bp recall and F1
- ARI
- HQ/MQ MAG count
- CheckM2 and GUNC for real data
- pair-edge false-positive rate
- chimeric-bin bases
- unbinned bp
- runtime / peak RAM / optional GPU time

## Bridge-specific long-term advantage

The main opportunity is not to make another stand-alone k-mer + coverage binner. Bridge controls the assembler and therefore can preserve evidence most external binners never receive:

```
reads -> assembly graph -> conservative bins
      -> bin-conditioned graph/path resolution
      -> local reassembly / scaffold update
      -> re-bin
```

Binning can help resolve A-X-B versus A-X-D ambiguities; improved graph paths then provide stronger physical linkage back to binning. This assembly-binning feedback loop is the intended production direction after the v2 evidence model is validated.

## Selected references / implementations studied

- MetaBAT2: Kang et al., PeerJ 2019.
- VAMB: Nissen et al., Nature Biotechnology 2021; current RasmussenLab/vamb implementation.
- SemiBin2: Pan et al., Bioinformatics 2023; current BigDataBiology/SemiBin implementation.
- COMEBin: current ziyewang/COMEBin implementation and contrastive-learning pipeline.
- TaxVAMB: RasmussenLab/VAMB, Nature Biotechnology 2026.
- QuickBin: BBTools implementation and Communications Biology 2026 paper.
- MetaCAT: Nature Microbiology 2026 and liu-congcong/MetaCAT implementation.
- LorBin: LorMeBioAI/LorBin and Nature Communications 2025.
- GENERanno: GenerTeam/GENERanno, bioRxiv 2025 and current prokaryotic CDS annotator checkpoints.
- Constrained clustering overview: Artificial Intelligence Review 2025.
- Correlation clustering / signed graph literature for future global partitioning.
