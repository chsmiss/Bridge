# Biological Brain matched Zymo Log 50k ablation

## Provenance

- Dataset: `ERR2935805`, deterministic first 50,000 paired reads.
- Source nucleotide graph: `bridgeasm-k31/assembly.gfa` from GitHub Actions run `33306726454`.
- Matched Biological Brain run: `33312388063`.
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

## Protein-recall ablation from the three path-cover outputs

Reference set: 79,971 predicted reference proteins / 17,471,191 amino acids. Thresholds are 80% identity, 80% reference coverage, and 80% target coverage for reciprocal completeness.

| output | predicted proteins | >=100 aa | single recall | union recall | reciprocal recall | reference AA coverage |
|---|---:|---:|---:|---:|---:|---:|
| physical only | 3,706 | 1,090 | 0.1976% | 0.3089% | 0.1926% | 1.8684% |
| physical + PLASS | 3,706 | 1,092 | 0.1976% | 0.3089% | 0.1926% | 1.8686% |
| physical + PLASS + DNA-LM | 3,707 | 1,093 | 0.1976% | 0.3089% | 0.1926% | 1.8686% |

Protein-guided reranking therefore did not materially increase protein recall in this experiment.

## Edge-level effect

PLASS produced 1,573 protein records; exact deduplication left 1,499 unique sequences and exact containment collapse retained 1,292 representatives.

Protein evidence over existing GFA links was classified as:

- `same_orf_supported`: 127 edges;
- `one_sided_protein_match`: 112;
- `ambiguous_homology`: 5;
- `no_protein_assembly_match`: 1,390;
- `stop_or_short_at_junction`: 2,562.

The conservative physical-only path cover selected 1,759 edges. PLASS retained the same number but replaced 8 selected edges with 8 alternatives (16-edge symmetric difference). Of the 127 `same_orf_supported` edges, 73 were selected without protein ranking and 79 with PLASS ranking. Adding the Markov DNA-LM changed four more selected edges.

## Decision

**Do not promote protein-guided path selection as a primary continuity/completeness mechanism.** On this matched test it changes a small number of ambiguous decisions but produces no measurable GF, N50, or complete-protein-recall gain.

Retain the Biological Brain architecture for two narrower roles:

1. **ambiguous-edge adjudication / veto** after physical evidence has produced candidate joins;
2. **parallel protein recovery**, where PLASS/PenguiN-style protein assemblies are reported as an independent protein product rather than forcing their paths back into nucleotide primary contigs.

The next matched workflow extends the protein evaluation to compare PLASS standalone directly against proteins predicted from BridgeAsm k31, MEGAHIT, metaSPAdes, and the Biological Brain path outputs. That result determines whether protein preassembly is valuable as an independent recovery layer even though it does not improve nucleotide continuity here.

ESM remains an optional bounded reranker with synthetic smoke validation. A Zymo ESM ablation should be restricted to the small set of physically eligible ambiguous edges; scoring every graph edge is not justified by the present result.
