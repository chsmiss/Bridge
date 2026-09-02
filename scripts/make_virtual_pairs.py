#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip
from pathlib import Path

_COMP = str.maketrans('ACGTNacgtn','TGCANtgcan')

def rc(s:str)->str: return s.translate(_COMP)[::-1].upper()

def fasta_records(path:Path):
    h=None; chunks=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: yield h,''.join(chunks).upper()
                h=line[1:]; chunks=[]
            else: chunks.append(line)
    if h is not None: yield h,''.join(chunks).upper()

def emit_pair(o1,o2,idx,seq,pos,read_len,insert):
    left=seq[pos:pos+read_len]
    rstart=pos+insert-read_len
    right=rc(seq[rstart:rstart+read_len])
    if len(left)!=read_len or len(right)!=read_len or 'N' in left or 'N' in right: return False
    q='I'*read_len
    ident=f'virtual_{idx:09d}_{pos}'
    o1.write(f'@{ident}/1\n{left}\n+\n{q}\n')
    o2.write(f'@{ident}/2\n{right}\n+\n{q}\n')
    return True

def main():
    ap=argparse.ArgumentParser(description='Project trusted contigs into the next k as low-copy virtual paired reads.')
    ap.add_argument('fasta', type=Path, nargs='+')
    ap.add_argument('--read1', type=Path, required=True)
    ap.add_argument('--read2', type=Path, required=True)
    ap.add_argument('--read-length', type=int, default=101)
    ap.add_argument('--insert-size', type=int, default=250)
    ap.add_argument('--stride', type=int, default=180)
    ap.add_argument('--min-length', type=int, default=500)
    ap.add_argument('--max-pairs-per-record', type=int, default=1000)
    args=ap.parse_args()
    if args.insert_size < args.read_length*2: raise SystemExit('insert-size must be >= 2*read-length')
    args.read1.parent.mkdir(parents=True,exist_ok=True); args.read2.parent.mkdir(parents=True,exist_ok=True)
    seen=set(); idx=0; pairs=0; bases=0
    with gzip.open(args.read1,'wt',compresslevel=3) as o1, gzip.open(args.read2,'wt',compresslevel=3) as o2:
        for fp in args.fasta:
            if not fp.exists(): continue
            for _h,seq in fasta_records(fp):
                if len(seq)<max(args.min_length,args.insert_size) or seq in seen: continue
                seen.add(seq); bases += len(seq); idx += 1
                emitted=0
                last=max(0,len(seq)-args.insert_size)
                positions=list(range(0,last+1,args.stride))
                if positions and positions[-1]!=last: positions.append(last)
                for pos in positions:
                    if emitted>=args.max_pairs_per_record: break
                    if emit_pair(o1,o2,idx,seq,pos,args.read_length,args.insert_size):
                        emitted += 1; pairs += 1
    print(f'records\t{len(seen)}')
    print(f'source_bases\t{bases}')
    print(f'virtual_pairs\t{pairs}')

if __name__=='__main__': main()
