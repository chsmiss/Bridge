#!/usr/bin/env python3
import argparse
import gzip
import math
from collections import defaultdict
from pathlib import Path


def open_text(path):
    path = Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, 'rt')
    return path.open()


def read_fasta_lengths(path):
    lengths = {}
    name = None
    n = 0
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
        if name is not None:
            lengths[name] = n
    return lengths


def read_coverage(path):
    names = []
    rows = {}
    first = True
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if first:
                first = False
                try:
                    float(fields[1])
                except ValueError:
                    names = fields[1:]
                    continue
            if not names:
                names = [f'sample{i+1}' for i in range(len(fields) - 1)]
            rows[fields[0]] = [float(x) for x in fields[1:]]
    return names, rows


def metabat_depth(args):
    lengths = read_fasta_lengths(args.fasta)
    samples, coverage = read_coverage(args.coverage)
    with open(args.output, 'w') as out:
        header = ['contigName', 'contigLen', 'totalAvgDepth']
        for sample in samples:
            header.extend([sample, f'{sample}-var'])
        print('\t'.join(header), file=out)
        for contig, length in lengths.items():
            row = coverage.get(contig, [0.0] * len(samples))
            avg = sum(row) / len(row) if row else 0.0
            fields = [contig, str(length), f'{avg:.8f}']
            for depth in row:
                # The benchmark has mean depth but no per-base variance. Use the
                # Poisson variance implied by that mean rather than leaking truth.
                fields.extend([f'{depth:.8f}', f'{max(depth, 1e-6):.8f}'])
            print('\t'.join(fields), file=out)


def iter_fasta_ids(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                yield line[1:].split()[0]


def bins_to_assignments(args):
    lengths = read_fasta_lengths(args.fasta)
    membership = {}
    paths = []
    root = Path(args.bins)
    for pattern in args.pattern:
        paths.extend(root.glob(pattern))
    paths = sorted(set(paths))
    for bin_index, path in enumerate(paths):
        label = f'bin_{bin_index + 1:04d}'
        for contig in iter_fasta_ids(path):
            if contig in membership:
                raise SystemExit(f'duplicate contig across bins: {contig}')
            membership[contig] = label
    with open(args.output, 'w') as out:
        print('contig\tbin\tlength\tscore', file=out)
        for contig, length in lengths.items():
            print(f'{contig}\t{membership.get(contig, "unbinned")}\t{length}\t{1.0 if contig in membership else 0.0:.6f}', file=out)


def read_assignments(path):
    assignments = {}
    lengths = {}
    with open(path) as fh:
        header = next(fh, None)
        for line in fh:
            fields = line.rstrip().split('\t')
            if len(fields) < 3:
                continue
            contig, bin_name, length = fields[:3]
            assignments[contig] = None if bin_name == 'unbinned' else bin_name
            lengths[contig] = int(length)
    return assignments, lengths


def weighted_median(observations):
    if not observations:
        return 0.0
    observations = sorted(observations)
    total = sum(w for _, w in observations)
    target = (total + 1) // 2
    acc = 0
    for value, weight in observations:
        acc += weight
        if acc >= target:
            return value
    return observations[-1][0]


def quantify(args):
    samples, coverage = read_coverage(args.coverage)
    assignments, lengths = read_assignments(args.assignments)
    bins = sorted({b for b in assignments.values() if b is not None})
    obs = {b: [[] for _ in samples] for b in bins}
    for contig, bin_name in assignments.items():
        if bin_name is None or contig not in coverage:
            continue
        length = lengths[contig]
        for i, depth in enumerate(coverage[contig]):
            obs[bin_name][i].append((depth, length))
    robust = {b: [weighted_median(x) for x in obs[b]] for b in bins}
    totals = [sum(robust[b][i] for b in bins) for i in range(len(samples))]
    with open(args.output, 'w') as out:
        print('bin\tsample\trobust_depth\tmean_depth\trelative_abundance\tcovered_bp\tcovered_contigs', file=out)
        for b in bins:
            for i, sample in enumerate(samples):
                vals = obs[b][i]
                bp = sum(length for _, length in vals)
                mean = sum(depth * length for depth, length in vals) / bp if bp else 0.0
                rel = robust[b][i] / totals[i] if totals[i] else 0.0
                print(f'{b}\t{sample}\t{robust[b][i]:.6f}\t{mean:.6f}\t{rel:.8f}\t{bp}\t{len(vals)}', file=out)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('metabat-depth')
    p.add_argument('--fasta', required=True)
    p.add_argument('--coverage', required=True)
    p.add_argument('--output', required=True)
    p.set_defaults(func=metabat_depth)

    p = sub.add_parser('bins-to-assignments')
    p.add_argument('--fasta', required=True)
    p.add_argument('--bins', required=True)
    p.add_argument('--pattern', action='append', default=['*.fa', '*.fasta', '*.fa.gz', '*.fna', '*.fna.gz'])
    p.add_argument('--output', required=True)
    p.set_defaults(func=bins_to_assignments)

    p = sub.add_parser('quantify')
    p.add_argument('--coverage', required=True)
    p.add_argument('--assignments', required=True)
    p.add_argument('--output', required=True)
    p.set_defaults(func=quantify)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
