#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path

SAMPLES = [f"sample{i}" for i in range(1, 7)]
PROFILES = {
    "Bacillus_subtilis": [20.0, 4.0, 15.0, 6.0, 28.0, 10.0],
    "Enterococcus_faecalis": [8.0, 25.0, 5.0, 18.0, 7.0, 30.0],
    "Escherichia_coli": [30.0, 12.0, 25.0, 7.0, 15.0, 20.0],
    # Deliberately identical to E. coli: this pair must be separated by sequence composition.
    "Salmonella_enterica": [30.0, 12.0, 25.0, 7.0, 15.0, 20.0],
    "Lactobacillus_fermentum": [12.0, 20.0, 8.0, 25.0, 5.0, 15.0],
    # Deliberately similar to Lactobacillus to make coverage less decisive.
    "Listeria_monocytogenes": [10.0, 18.0, 6.0, 23.0, 8.0, 16.0],
    "Pseudomonas_aeruginosa": [5.0, 8.0, 30.0, 12.0, 25.0, 6.0],
    "Staphylococcus_aureus": [25.0, 6.0, 10.0, 30.0, 12.0, 8.0],
}


def read_fasta(path):
    records = []
    name = None
    seq = []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq).upper()))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def discover_reference(ref_dir, species):
    root = Path(ref_dir)
    for suffix in (".fasta", ".fa", ".fna"):
        candidate = root / f"{species}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"missing reference for {species} in {root}")


def fragment_records(records, rng, min_len, max_len):
    fragments = []
    for record_name, seq in records:
        start = 0
        index = 0
        while start < len(seq):
            target = rng.randint(min_len, max_len)
            end = min(len(seq), start + target)
            fragment = seq[start:end]
            if fragment:
                fragments.append((record_name, index, start, end, fragment))
            start = end
            index += 1
    return fragments


def main():
    parser = argparse.ArgumentParser(
        description="Create a deterministic BridgeBin benchmark from real Zymo reference genomes"
    )
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=43)
    parser.add_argument("--min-fragment", type=int, default=2200)
    parser.add_argument("--max-fragment", type=int, default=8500)
    parser.add_argument("--min-contig", type=int, default=1500)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    master = random.Random(args.seed)
    truth_rows = []
    coverage_rows = []

    with (out / "contigs.fa").open("w") as fasta:
        for species in sorted(PROFILES):
            ref = discover_reference(args.reference_dir, species)
            rng = random.Random(master.getrandbits(64))
            fragments = fragment_records(
                read_fasta(ref), rng, args.min_fragment, args.max_fragment
            )
            for serial, (record, _, start, end, seq) in enumerate(fragments, start=1):
                contig = f"{species}__c{serial:05d}"
                fasta.write(f">{contig}\n")
                for pos in range(0, len(seq), 80):
                    fasta.write(seq[pos : pos + 80] + "\n")
                eligible = int(len(seq) >= args.min_contig)
                truth_rows.append((contig, species, len(seq), eligible, record, start, end))

                contig_scale = 1.0 + master.uniform(-0.10, 0.10)
                depths = []
                for base_depth in PROFILES[species]:
                    sample_noise = 1.0 + master.uniform(-0.025, 0.025)
                    depths.append(base_depth * contig_scale * sample_noise)
                coverage_rows.append((contig, depths))

    with (out / "coverage.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", *SAMPLES])
        for contig, depths in coverage_rows:
            writer.writerow([contig, *[f"{x:.6f}" for x in depths]])

    with (out / "truth.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", "genome", "length", "eligible", "record", "start", "end"])
        writer.writerows(truth_rows)

    with (out / "truth_abundance.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome", "sample", "depth", "relative_abundance"])
        for sample_index, sample in enumerate(SAMPLES):
            total = sum(profile[sample_index] for profile in PROFILES.values())
            for species in sorted(PROFILES):
                depth = PROFILES[species][sample_index]
                writer.writerow(
                    [species, sample, f"{depth:.6f}", f"{depth / total:.8f}"]
                )

    eligible = sum(row[3] for row in truth_rows)
    eligible_bp = sum(row[2] for row in truth_rows if row[3])
    print(
        f"generated {len(truth_rows)} Zymo reference fragments; "
        f"eligible={eligible}; eligible_bp={eligible_bp}; species={len(PROFILES)}"
    )


if __name__ == "__main__":
    main()
