# Biological Brain matched Zymo Log 50k ablation

## Provenance

- Dataset: `ERR2935805`, deterministic first 50,000 paired reads.
- Source nucleotide graph: `bridgeasm-k31/assembly.gfa` from GitHub Actions run `33306726454`.
- Final matched Biological Brain run: `33317630673` (`success`).
- References: ZymoBIOMICS reference collection used by the existing reference-smoke workflow.
- Protein preassembly: PLASS, followed by exact duplicate and exact-contained protein collapse.
- Path mode: `conservative`; protein and DNA-LM signals only rerank physically eligible links that already exist in the GFA.

This is a matched ablation. It is intentionally separate from the default BridgeAsm primary emitter and must not be interpreted as a replacement production assembly benchmark.

## Nucleotide ablation

| output | GF (%) | N50 | NA50 | largest | misassemblies | mismatches / 100 kbp | indels / 100 kbp |
|---|---:|---:|---:|---:|---:|---:|---:|
| physical only | 1.501 | 338 | 334 | 3,704 | 50 | 280.69 | 53.96 |
| physical + PLASS | 1.501 | 338 | 334 | 3,704 | 48 | 280.63 | 53.95 |
| physical + PLASS + DNA-LM | 1.501 | 338 | 334 | 3,704 | 47 | 280.62 | 53.93 |

PLASS changed the selected path-cover edge set but did not change GF or N50. It reduced reported misassemblies by two in this matched path-cover experiment; the simple order-5 DNA-LM changed four additional edge decisions and reduced the count by one more. These are small correctness signals, not evidence of a completeness or contiguity breakthrough.

The high mismatch/indel rates in this table belong to the experimental `bridgeasm-evidence-path` path cover over the raw unitig graph. They are not directly comparable to the default BridgeAsm k31 primary output, which on the same 50k source run had GF 1.427%, N50/NA50 328 bp, zero reported misassemblies, 2.17 mismatches/100 kbp and 0.59 indels/100 kbp. MEGAHIT on that source run had GF 2.284%, N50 623 bp and 79 misassemblies; metaSPAdes had GF 3.169%, N50 579 bp and 47 misassemblies.

## Final matched protein-recall comparison

Reference set: 79,971 predicted reference proteins / 17,471,191 amino acids. Thresholds are 80% identity, 80% reference coverage, and 80% target coverage for reciprocal completeness.

| output | predicted proteins | >=100 aa | single recall | union recall | reciprocal recall | reference AA coverage | matched predicted |
|---|---:|---:|---:|---:|---:|---:|---:|
| BridgeAsm k31 primary | 3,392 | 998 | 0.1788% | 0.2839% | 0.1751% | 1.8185% | 96.93% |
| MEGAHIT | 3,830 | 2,420 | 0.6715% | 1.0541% | 0.6640% | 2.9488% | 97.34% |
| metaSPAdes | 5,823 | 2,857 | **0.8666%** | **1.8732%** | **0.8541%** | **4.1159%** | 96.43% |
| PLASS standalone | 1,292 | 163 | 0.0263% | 0.0300% | 0.0263% | 0.6030% | 89.86% |
| physical only | 3,706 | 1,090 | 0.1976% | 0.3089% | 0.1926% | 1.8684% | 96.74% |
| physical + PLASS | 3,706 | 1,092 | 0.1976% | 0.3089% | 0.1926% | 1.8686% | 96.76% |
| physical + PLASS + DNA-LM | 3,707 | 1,093 | 0.1976% | 0.3089% | 0.1926% | 1.8686% | 96.76% |

On this shallow 50k-pair subset, the hypothesized standalone PLASS recovery advantage is **not observed**. PLASS standalone recovers substantially fewer complete reference proteins and less reference amino-acid sequence than all three nucleotide-first assemblers. metaSPAdes is the strongest protein-recall baseline in this test.

The experimental physical path cover has a small protein-recall advantage over the default BridgeAsm k31 primary output, but PLASS reranking itself does not add measurable complete-protein recall. The gain is therefore attributable to the different nucleotide path cover, not to protein guidance.

This result should not be generalized to deep metagenomes yet: 50,000 pairs is an unusually shallow input for protein-level assembly, PLASS retained only 1,292 representatives after exact containment collapse, and only 163 were >=100 aa. A 500k/full-log protein benchmark is required before rejecting PLASS/PenguiN as a parallel low-abundance protein-recovery layer.

## Edge-level effect

PLASS produced 1,573 protein records; exact deduplication left 1,499 unique sequences and exact containment collapse retained 1,292 representatives.

Protein evidence over existing GFA links was classified as:

- `same_orf_supported`: 127 edges;
- `one_sided_protein_match`: 112;
- `ambiguous_homology`: 5;
- `no_protein_assembly_match`: 1,390;
- `stop_or_short_at_junction`: 2,562.

The conservative physical-only path cover selected 1,759 edges. PLASS retained the same number but replaced 8 selected edges with 8 alternatives (16-edge symmetric difference). Of the 127 `same_orf_supported` edges, 73 were selected without protein ranking and 79 with PLASS ranking. Adding the Markov DNA-LM changed four more selected edges.

## Production decision after the completed 50k benchmark

1. **Do not promote protein-guided path selection as a primary continuity/completeness mechanism.** It does not improve GF, N50, or complete-protein recall on this matched benchmark.
2. **Keep protein evidence as an ambiguous-edge adjudicator/veto.** It changed a small number of physically eligible decisions and reduced the experimental path-cover misassembly count from 50 to 48 without changing GF.
3. **Keep DNA-LM as a bounded secondary reranker only.** It changed four more edges and reduced the count from 48 to 47, again without GF/N50 gain.
4. **Do not claim PLASS standalone improves protein recall from this 50k experiment.** The opposite is observed here.
5. **Run the protein-recovery hypothesis again at 500k and Full Log.** That is the appropriate scale to test whether amino-acid-space assembly recovers low-abundance coding sequence that nucleotide contig assembly misses.
6. **ESM should be evaluated only on the small physically eligible ambiguous-edge set**, preferably as a veto/ranking signal against known wrong alternatives, not across every graph edge and never as a source of new nucleotide sequence.

The resulting architecture remains:

```text
nucleotide DBG / iterative multi-k / physical read evidence
                         |
                         +--> primary nucleotide assembly
                         |
                         +--> ambiguous candidate edges
                                  |
                                  +--> PLASS/PenguiN continuity evidence
                                  +--> optional ESM breakpoint score
                                  +--> bounded DNA-LM score
                                  |
                                  +--> rerank or veto existing candidates only

raw reads --------------------------> parallel protein assembly product
```

The next evidence-bearing experiment is the 500k matched protein benchmark, not another 50k parameter sweep.
