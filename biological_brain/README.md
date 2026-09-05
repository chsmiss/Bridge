# Biological Brain: protein- and language-model-assisted graph decisions

This module adds biological priors to **existing** BridgeAsm graph edges.  It is designed
to improve branch choice and protein recovery without allowing a language model to
hallucinate nucleotide sequence.

## Decision hierarchy

1. **Hard nucleotide constraints**
   - an edge must already exist in the GFA;
   - overlap, direct-read, gapped-read, paired-end, topology, coverage, and cycle checks
     remain authoritative;
   - no protein or language model may synthesize a missing bridge in the default mode.
2. **Protein-assembly continuity**
   - Plass/Penguin-style protein assemblies provide long amino-acid paths;
   - each existing GFA edge is reconstructed locally, translated in all six frames, and
     tested for a stop-free peptide that crosses the exact nucleotide junction;
   - exact amino-acid k-mers on both sides of the breakpoint must map to one protein
     assembly on one diagonal;
   - duplicated/homologous matches are assigned high ambiguity and are not used as
     strong evidence.
3. **Protein language-model adjudication**
   - ESM compares local pseudo-log-likelihood with the peptide halves joined versus
     independent;
   - it is useful for `same ORF`, frameshift/stop, domain boundary, and random chimera
     discrimination;
   - it has a small weight and is ignored unless a nucleotide/protein-supported edge is
     already eligible.
4. **DNA language-model adjudication**
   - GENERator/CARBon-like models should score candidate junction contexts rather than
     generate sequence;
   - the interface is a TSV with `source`, `target`, and `dna_lm_delta`;
   - `kmer_junction_lm.py` is a transparent order-5 Markov baseline that validates the
     adapter and calibration path before large neural checkpoints are introduced.

## Why protein preassembly can help

Nucleotide assembly breaks at strain variation, synonymous substitutions, sequencing
errors, and low-coverage branches.  Protein-level assembly projects reads into a smaller
functional alphabet and can recover a long coding path even when the nucleotide graph
contains several locally plausible alternatives.  Back-projecting that path to an
**already existing** GFA edge can answer a narrow and safe question:

> Does this exact nucleotide junction preserve one protein path across the breakpoint?

It cannot establish intergenic sequence, chromosome-scale order, plasmid host, or taxon
identity by itself.  Conserved proteins and paralogs are therefore treated as ambiguous,
not as permission to join genomes.

## Implemented files

- `protein_bridge_evidence.py`
  - reads a sequence-bearing GFA and Plass/Penguin-style protein FASTA;
  - reconstructs each plus/plus junction;
  - translates six frames and reports `same_orf_supported`,
    `one_sided_protein_match`, `ambiguous_homology`,
    `stop_or_short_at_junction`, or `no_protein_assembly_match`;
  - emits edge scores and the junction peptide for optional ESM scoring.
- `esm_breakpoint_score.py`
  - optional small-ESM masked pseudo-likelihood scorer;
  - lazy-loads `torch`/`transformers`, so normal BridgeAsm builds stay lightweight.
- `kmer_junction_lm.py`
  - deterministic DNA-LM adapter baseline;
  - external neural DNA models should emit the same scalar-score schema.
- `src/bin/bridgeasm-evidence-path.rs`
  - fuses nucleotide, protein, ESM, and DNA-LM evidence;
  - keeps one-successor/one-predecessor path-cover and cycle prevention;
  - outputs a full per-edge audit report.

## Evidence TSV contract

The required protein evidence columns are:

```text
source  target  protein_score  unique_kmers  ambiguity
frame_consistency  protein_id  breakpoint_class
```

Optional ESM score file:

```text
source  target  esm_delta
```

Optional DNA model score file:

```text
source  target  dna_lm_delta
```

All scores are keyed by the exact GFA segment names.  Unknown rows are ignored.  Missing
scores are zero.

## Recommended experiment ladder

### A. Deterministic branch fixture

The committed unit test contains one true coding continuation and one coverage-matched
decoy.  It verifies that:

- only existing GFA links are scored;
- the true edge has matched amino-acid k-mers on both sides;
- the decoy is one-sided and receives zero usable protein score;
- an exact duplicate reference protein is marked ambiguous;
- the Rust path cover selects the true edge.

### B. Real Plass integration fixture

GitHub Actions generates overlapping nucleotide reads, runs the real Plass executable,
uses its assembled proteins as evidence, and checks that the nucleotide graph chooses
the true branch.  This tests the full translation/assembly/back-projection boundary.

### C. Zymo pilot

Run four arms on the same graph and same reads:

1. nucleotide-only baseline;
2. + Plass evidence;
3. + Plass + ESM;
4. + Plass + ESM + DNA-LM score.

Report both assembly and protein metrics:

- N50, NA50, largest contig, genome fraction, misassemblies, mismatch/indel rate;
- predicted ORFs, full-length ORF fraction, reference-protein recall at several identity
  and coverage thresholds;
- false cross-reference joins and duplicated single-copy markers;
- edge-level precision/recall for branches that can be labelled from the references;
- results stratified by abundance, coverage, protein length, conserved-domain content,
  and strain-shared sequence.

A higher N50 is accepted only when NA50 and misassembly/error metrics remain controlled.

## Default safeguards

- No model-created GFA edges.
- No generated nucleotide sequence.
- Protein support must cross the exact breakpoint and have evidence on both sides.
- Conserved/duplicated proteins are ambiguity-gated.
- ESM and DNA-LM scores are bounded with `tanh` and have much smaller weights than
  direct-read/protein evidence.
- Every selected edge remains auditable in the report.

## Next production step

For sensitivity beyond exact amino-acid k-mers, add a Plass/Penguin alignment adapter
that emits the same TSV from MMseqs2-style local alignments.  The path selector does not
need to change.  Model-specific GENERator/CARBon adapters should likewise be isolated
behind the scalar junction-score contract so checkpoint, context length, license, and
calibration changes cannot silently alter the assembler core.
