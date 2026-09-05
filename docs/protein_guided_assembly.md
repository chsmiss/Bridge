# Protein-guided assembly and language-model evidence

## Decision

BridgeAsm should keep nucleotide assembly as the source of sequence truth. Protein
assembly and language models may propose, rank, or veto candidate joins, but they
must not invent unsupported nucleotide sequence in the default mode.

The integration order is:

1. **PenguiN-guided nucleotide contigs** propose long-range ordering and short
   bridge sequences.
2. **PLASS protein assemblies** provide an independent protein-recovery benchmark
   and later an orthogonal edge-support signal.
3. **ESM-C** scores translated coding junctions only.
4. **GENERator or Carbon** compare already observed nucleotide alternatives. They
   are not used to free-generate gaps.

This order makes the first experiment deterministic, CPU-testable, auditable, and
possible to evaluate with reference genomes before adding expensive model inference.

## Why PLASS and PenguiN are useful

PLASS assembles six-frame-translated read fragments in amino-acid space. Protein
sequence is more conserved than DNA sequence and removes synonymous variation, so
protein-level overlaps can recover coding regions that are fragmented in a strict
nucleotide graph. Its output cannot recover non-coding sequence or uniquely resolve
synonymous nucleotide paths.

PenguiN is the closer fit for nucleotide assembly. Its `guided_nuclassemble`
workflow extracts and translates ORFs, finds amino-acid overlaps, maps those
alignments back to nucleotide coordinates, extends nucleotide and protein sequences,
then performs nucleotide assembly and redundancy reduction. BridgeAsm therefore uses
PenguiN nucleotide contigs as *guides*, while preserving its own contigs as immutable
backbone segments.

Upstream implementation inspected for this design:

- <https://github.com/soedinglab/plass>
- `data/guidedNuclAssemble.sh`
- `src/workflow/GuidedNuclassembler.cpp`
- `src/assembler/guidedassembleresult.cpp`
- `src/assembler/nuclassembleresult.cpp`

## Data flow

```text
paired reads
   |-------------------------> PLASS proteins
   |                                |
   |                                +--> protein recall benchmark
   |                                +--> future same-protein edge support
   |
   |-------------------------> PenguiN guided nucleotide contigs
   |                                |
   v                                v
BridgeAsm nucleotide contigs --> unique PAF anchors --> candidate joins
                                                    |
                                                    +--> hard sequence filters
                                                    +--> optional ESM score
                                                    +--> optional DNA-LM score
                                                    v
                                             acyclic path cover
                                                    v
                                  consensus FASTA + evidence TSV + JSON
```

For the first benchmark, `primary_contigs.fasta` is converted to a link-free GFA.
Every existing contig becomes one segment, so the protein-guided stage can join
contigs but cannot split or rewrite them.

## Implemented hard filters

`bridgeasm-proteinguide` currently requires:

- a primary, high-map-quality PAF anchor;
- sufficient aligned length, identity, and query coverage;
- no near-equal mapping to a different guide locus;
- monotonic ordering of anchors on a guide contig;
- orientation consistency;
- either a direct GFA edge, a validated suffix-prefix overlap, or a bounded guide
  gap;
- low ambiguity in any inserted guide bases;
- compatible segment coverage when coverage is available;
- consistent gap estimates when more than one guide supports an edge;
- an acyclic one-predecessor/one-successor path cover.

Three modes are exposed:

- `overlap-only`: never inserts guide bases;
- `conservative`: permits only short gaps between nearly full-length, very
  high-identity anchors;
- `protein`: permits larger PenguiN-supported gaps and remains experimental.

Every candidate is written to a TSV report with selection state, rejection reason,
identity, alignment fraction, map quality, projected gap, overlap, guide support,
coverage ratio, and optional model evidence.

## PLASS evidence beyond benchmarking

The next protein-specific edge scorer should operate on candidate junctions rather
than directly altering the graph:

1. call terminal ORFs on both sides of each candidate;
2. search those ORFs against the PLASS assembly;
3. support an edge when both sides align to the same PLASS protein in compatible
   order and frame;
4. veto an edge when the proposed join introduces an internal stop or contradicts
   a high-confidence PLASS protein;
5. do not apply this test to non-coding candidates.

This retains proteins recovered only by PLASS without forcing their amino-acid
consensus back into an arbitrary nucleotide sequence.

## ESM-C breakpoint scoring

ESM is useful only after a coding frame has been established. For each candidate:

1. predict ORFs on the joined nucleotide window;
2. retain ORFs that cross the candidate boundary;
3. translate the reference left and right context plus each competing join;
4. compute a junction-local pseudo-log-likelihood or a calibrated classifier over
   ESM-C representations;
5. compare candidate scores within the same locus;
6. emit a normalized score or `veto` in the model-evidence TSV.

A valid ESM scorer must be calibrated on held-out true and false breaks generated
from complete microbial genomes. It should include explicit features for internal
stops, frameshifts, low-complexity sequence, ORF length, and distance from the
junction. ESM must not rank candidates for which no coding ORF crosses the boundary.

## GENERator and Carbon

Both models are autoregressive genomic models and can be used for *relative
likelihood*, not unrestrained gap generation.

For a locus with candidate sequences `c_1 ... c_n`, score the same left context and
right context under each candidate:

```text
DNA_score(c) = mean_log_probability(left + c + right)
             - mean_log_probability(shuffled_or_null_control)
```

Required safeguards:

- score observed candidate paths only;
- score both forward and reverse-complement representations;
- use identical context lengths and token-boundary handling;
- normalize by candidate length and local GC/codon composition;
- reject any scorer with strong strand asymmetry on held-out controls;
- calibrate scores per organism abundance and coding/non-coding class;
- require independent read or guide evidence before an LM can affect selection.

GENERator uses a 6-mer tokenizer, so all compared contexts must receive identical
left padding or truncation. Carbon also uses 6-mer DNA tokens and requires its DNA
input tag. For Zymo bacterial data, the prokaryotic GENERator-v2 checkpoint is a
more natural first DNA-LM experiment than eukaryote-focused checkpoints. Carbon is
valuable as an independent model after the deterministic protein-guided result is
established.

## Model evidence interface

The Rust selector accepts an optional TSV:

```text
source  target  score  decision  scorer
u12+    u57-    0.83   neutral   esmc-600m-junction-v1
u91+    u33+   -1.20   veto      generator-prokaryote-v2
```

- `source` and `target` are oriented segment identifiers;
- scores are averaged across votes and added only after hard filters;
- `veto` or `reject` disqualifies a candidate;
- model evidence never creates a candidate or supplies nucleotide bases.

## Validation matrix

The Zymo Log smoke benchmark compares:

| Output | Purpose |
|---|---|
| BridgeAsm primary contigs | nucleotide baseline |
| PenguiN guided nucleotide assembly | external protein-guided baseline |
| PLASS proteins | protein-only recall ceiling/control |
| BridgeAsm + overlap-only guide | sequence-only safe join control |
| BridgeAsm + conservative guide | proposed default |
| BridgeAsm + protein guide | experimental upper-contiguity mode |
| MEGAHIT | established nucleotide comparator |

Nucleotide metrics:

- total aligned fraction / genome fraction;
- N50 and NA50;
- largest alignment;
- extensive and local misassemblies;
- mismatches and indels per 100 kbp;
- duplicated ratio and unaligned sequence.

Protein metrics:

- reference proteins recovered by one high-identity hit;
- reference proteins recovered by the union of fragments;
- reciprocal-complete protein recovery;
- fraction of reference amino acids covered;
- number and total length of predicted proteins;
- PLASS-only proteins absent from every nucleotide assembly.

## Promotion criteria

`conservative` should not become the default unless it is Pareto-improving on more
than one dataset. The initial acceptance gates are:

- no material loss of aligned genome fraction;
- extensive misassemblies no worse than the current production baseline by more
  than the predeclared tolerance;
- higher single-hit or reciprocal protein recall;
- a meaningful N50/NA50 increase that remains after reference correction;
- inserted guide bases remain a small, explicitly reported fraction of output;
- gains reproduce on a second mock community or simulated abundance profile.

The `protein` mode remains experimental even if it yields a high N50. A high N50
without NA50, misassembly, and protein-recall improvement is not considered a win.

## Repository components

- `src/bin/bridgeasm-proteinguide.rs`: auditable guide-anchor path cover;
- `scripts/fasta_to_gfa.py`: immutable-contig backbone converter;
- `scripts/protein_recall.py`: single-hit, union, and reciprocal protein recall;
- `.github/workflows/zymo-proteinguide-smoke.yml`: real-data comparison workflow.
