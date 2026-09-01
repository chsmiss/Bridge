# BridgeAsm engineering optimization plan

This note records performance work that is reference-free and must preserve assembly decisions unless explicitly described otherwise.

## Validated 50k diagnostic

Dataset: ERR2935805, 50,000 paired-end spots selected as five seeded non-overlapping windows (seed 43), matching the existing random-50k Bridge benchmarks.

The `kmer_engineering_bench` diagnostic compares three independent whole-read scans using the current 5x-u64 `KmerKey` with one scan that updates compact counters for k=21,25,31.

| mode | wall seconds | peak RSS MiB | k-mer key bytes |
| --- | ---: | ---: | ---: |
| repeated wide keys | 6.93 | 592.7 | 40 |
| single-pass compact | 1.87 | 412.6 | 8 for k<=31 |

The observation and distinct-k-mer counts are exactly identical for every tested k:

- k21: 8,094,658 observations; 3,796,352 distinct
- k25: 7,694,799 observations; 3,745,292 distinct
- k31: 7,095,080 observations; 3,644,423 distinct

This is a 3.71x wall-time improvement and a 30.4% lower peak RSS even though the compact test keeps all three maps live at once. The result validates both one-pass multi-k extraction and width-specialized keys as high-priority engineering work.

## P0: typed compact k-mer data plane

Do not use an enum as the hash-map key because its size is determined by the largest variant. Dispatch once at the outer assembly level and monomorphize the hot path:

- k <= 31: `u64`
- k <= 63: `u128`
- larger supported k: fixed multiword fallback

The typed key must flow through evidence counting, retained/rescued membership, graph indexing, and read threading; otherwise a conversion back to the 40-byte key removes much of the memory benefit.

## P0: decode reads once

The current assembly pipeline scans/decompresses FASTQ repeatedly for counting/filtering, mercy/singleton recovery, graph edge construction, and read threading. Introduce a packed `ReadStore` built in the first pass:

- 2-bit bases plus an ambiguity bitmap
- compact read offsets/lengths
- mates implicit from pair order or compact pair metadata
- quality bytes retained only while low-depth evidence scoring needs them
- no read-name strings in hot phases

Later phases replay the binary store rather than reparsing gzip/FASTQ.

For large inputs, the same format should support memory mapping so the assembler does not require the entire store in anonymous RAM.

## P0: consume and release k-mer evidence before edge counting

Graph construction currently overlaps the full `KmerSet` (`evidence`, `retained`, `rescued`) with a sorted key vector, key-to-node index, and edge-count table. Convert retained keys to node IDs, copy only counts and rescued flags required by graph construction, then release the large evidence/hash sets before scanning transitions.

Represent rescued/solid state as node-ID bit vectors after node assignment rather than another k-mer hash set.

## P1: one rolling max-k state for several exact k values

For a small candidate set such as 21/25/31, a read should be decoded once. A further optimization over independent rolling states is to maintain the largest rolling forward/reverse-complement word and derive smaller suffix k-mers with masks/shifts at each position. Hash insertion remains O(number of k values), but base decoding and rolling updates are shared.

For broad k sweeps, do not keep every whole-data table resident simultaneously.

## P1: signature/minimizer bins plus sort/reduce

Replace large per-k `FxHashMap` counting with a two-stage data plane:

1. scan the packed reads and emit compact k-mer/super-k-mer records to balanced signature/minimizer bins;
2. process bins independently with radix/sort + run-length reduction.

Benefits:

- bounded memory
- sequential access instead of random hash-table access
- easy per-thread partitioning without lock contention
- deterministic output
- a natural place to emit several k sizes from one read traversal

## P1: streaming k-mer iterator

Avoid allocating a `Vec` of canonical k-mers for every read in graph construction and threading. Expose a rolling iterator/callback yielding `(position, key, orientation)` and consume adjacent pairs online. Reuse per-thread scratch buffers for fragment-level deduplication.

## P1: integer-ID graph as early as possible

After retained k-mers receive node IDs:

- adjacency should use packed `u32` IDs
- edge counts can use packed `u64` `(source,target)` keys or sorted edge records
- rescued/solid flags should be bit vectors
- raw k-mer keys should not be duplicated in multiple maps

Keep raw keys only for phases that genuinely need sequence reconstruction. Consider streaming/debug output so the full raw graph does not have to remain in `AssemblyProduct` after unitigs and threading are complete.

## P1: phase-owned memory and instrumentation

Make large structures phase-owned and explicitly release them at handoff boundaries. Record phase-level RSS, allocation estimates, read-store bytes, and temporary-bin bytes in `run_profile.json`, in addition to timing. Optimization acceptance should use both maximum RSS and the phase responsible for it.

## P2: parallel data plane

Once reads are in `ReadStore`, partition by read ranges. Each worker emits to private bins or private sorted buffers, followed by deterministic reduction. Avoid one shared mutable k-mer hash table. The current `--threads` setting should accelerate counting and edge extraction, not only downstream work.

## P2: local-k routing from the read index

Local adaptive-k should query a persistent minimizer/signature-to-read index rather than rescan FASTQ or materialize neighborhood FASTQs. Reads can be soft-assigned by read IDs and decoded directly from `ReadStore`.

## Acceptance order

1. Preserve exact k-mer counts in unit tests and the engineering benchmark.
2. Preserve 50k MetaQUAST GF/misassembly metrics for pure engineering changes.
3. Measure wall/RSS on deterministic first-50k and random seed-43 50k.
4. Validate 500k before promoting the data plane to default.
5. Only then combine the faster data plane with new low-depth evidence scoring and graph-level local k-lifting.
