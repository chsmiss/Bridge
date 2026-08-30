#!/usr/bin/env python3
import argparse
import gzip


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def fasta_records(path):
    name = None
    seq = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq).upper()


def main():
    parser = argparse.ArgumentParser(description='Convert FASTA contigs to deterministic high-quality single-end FASTQ.')
    parser.add_argument('fasta')
    parser.add_argument('fastq')
    parser.add_argument('--min-length', type=int, default=1)
    parser.add_argument('--copies', type=int, default=2, help='Repeat each contig to provide solid k-mer multiplicity.')
    args = parser.parse_args()
    if args.copies < 1:
        raise SystemExit('--copies must be >= 1')
    written = 0
    with open(args.fastq, 'w') as out:
        for name, seq in fasta_records(args.fasta):
            if len(seq) < args.min_length or any(base not in 'ACGT' for base in seq):
                continue
            qual = 'I' * len(seq)
            for copy in range(args.copies):
                out.write(f'@iter_{written:08d}_{name}_c{copy + 1}\n{seq}\n+\n{qual}\n')
            written += 1
    print(f'converted {written} FASTA records to {written * args.copies} FASTQ records')


if __name__ == '__main__':
    main()
