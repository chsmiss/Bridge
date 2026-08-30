#!/usr/bin/env python3
from __future__ import annotations
import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

BASE = {65:0,67:1,71:2,84:3}
TRANS = bytes.maketrans(b'ACGTN', b'TGCAN')

def records(path: Path) -> Iterator[tuple[str, bytes]]:
    h=None; chunks=[]
    with path.open('rb') as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith(b'>'):
                if h is not None: yield h, b''.join(chunks).upper()
                h=line[1:].decode('utf-8','replace'); chunks=[]
            else:
                if h is None: raise ValueError(f'sequence before FASTA header in {path}')
                chunks.append(line)
    if h is not None: yield h,b''.join(chunks).upper()

def rc(seq: bytes)->bytes:
    return seq.translate(TRANS)[::-1]

def canonical(seq:bytes)->bytes:
    r=rc(seq); return r if r<seq else seq

def seed_key(seq:bytes, start:int, k:int)->int|None:
    x=0
    for b in seq[start:start+k]:
        v=BASE.get(b)
        if v is None:return None
        x=(x<<2)|v
    return x

def rolling_positions(seq:bytes,k:int,max_start:int):
    mask=(1<<(2*k))-1
    x=0; valid=0
    for i,b in enumerate(seq):
        v=BASE.get(b)
        if v is None:
            x=0; valid=0; continue
        x=((x<<2)|v)&mask; valid+=1
        if valid>=k:
            p=i+1-k
            if p>max_start: break
            yield p,x

@dataclass(frozen=True)
class Match:
    source:int; target:int; overlap:int

class DSU:
    def __init__(self,n): self.p=list(range(n)); self.sz=[1]*n
    def find(self,x):
        while self.p[x]!=x:
            self.p[x]=self.p[self.p[x]]; x=self.p[x]
        return x
    def union(self,a,b):
        a=self.find(a); b=self.find(b)
        if a==b:return False
        if self.sz[a]<self.sz[b]:a,b=b,a
        self.p[b]=a; self.sz[a]+=self.sz[b]; return True

def top_unique(values:dict[int,int], margin:int):
    if not values:return None
    ranked=sorted(values.items(), key=lambda kv:(-kv[1],kv[0]))
    if len(ranked)>1 and (ranked[0][1]==ranked[1][1] or ranked[0][1]-ranked[1][1]<margin):
        return None
    return ranked[0]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('output',type=Path)
    ap.add_argument('inputs',nargs='+',type=Path)
    ap.add_argument('--min-overlap',type=int,default=100)
    ap.add_argument('--overlap-margin',type=int,default=20)
    ap.add_argument('--seed-length',type=int,default=31)
    ap.add_argument('--max-seed-occurrences',type=int,default=128)
    ap.add_argument('--min-length',type=int,default=200)
    args=ap.parse_args()
    if args.seed_length>args.min_overlap: raise SystemExit('seed length > min overlap')
    seq_source={}
    for path in args.inputs:
        for h,s in records(path):
            if len(s)<args.min_length:continue
            c=canonical(s)
            seq_source.setdefault(c,set()).add(path.name)
    seqs=sorted(seq_source,key=lambda s:(-len(s),s))
    n=len(seqs)
    oriented=[]
    for s in seqs: oriented.extend((s,rc(s)))
    prefix=defaultdict(list)
    for state,s in enumerate(oriented):
        key=seed_key(s,0,args.seed_length)
        if key is not None: prefix[key].append(state)
    # Mark repetitive prefix seeds as unusable.
    prefix={k:v for k,v in prefix.items() if len(v)<=args.max_seed_occurrences}
    out_maps=[{} for _ in range(2*n)]
    in_maps=[{} for _ in range(2*n)]
    checks=0; matches=0
    for src,a in enumerate(oriented):
        max_start=len(a)-args.min_overlap
        if max_start<0:continue
        local={}
        for pos,key in rolling_positions(a,args.seed_length,max_start):
            candidates=prefix.get(key)
            if not candidates:continue
            ov=len(a)-pos
            suffix=None
            for tgt in candidates:
                if tgt//2==src//2:continue
                b=oriented[tgt]
                if len(b)<=ov:continue
                checks+=1
                if suffix is None:suffix=a[pos:]
                if b.startswith(suffix):
                    old=local.get(tgt,0)
                    if ov>old:local[tgt]=ov
        for tgt,ov in local.items():
            out_maps[src][tgt]=ov
            in_maps[tgt][src]=ov
            matches+=1
    best_out=[top_unique(m,args.overlap_margin) for m in out_maps]
    best_in=[top_unique(m,args.overlap_margin) for m in in_maps]
    reciprocal=[]
    seen_physical=set()
    for src,item in enumerate(best_out):
        if item is None:continue
        tgt,ov=item
        if best_in[tgt] != (src,ov):continue
        rev_src=src^1; rev_tgt=tgt^1
        physical=min((src,tgt),(rev_tgt,rev_src))
        if physical in seen_physical:continue
        seen_physical.add(physical)
        reciprocal.append(Match(src,tgt,ov))
    reciprocal.sort(key=lambda m:(-m.overlap,m.source,m.target))
    succ=[None]*(2*n); pred=[None]*(2*n); edge_ov={}; dsu=DSU(n); accepted=[]
    for m in reciprocal:
        s,t,ov=m.source,m.target,m.overlap
        rs,rt=s^1,t^1
        if succ[s] is not None or pred[t] is not None or succ[rt] is not None or pred[rs] is not None:continue
        if not dsu.union(s//2,t//2):continue
        succ[s]=t; pred[t]=s; edge_ov[(s,t)]=ov
        succ[rt]=rs; pred[rs]=rt; edge_ov[(rt,rs)]=ov
        accepted.append(m)
    used_records=set(); outputs=[]; chains=[]
    for state in range(2*n):
        rid=state//2
        if rid in used_records or pred[state] is not None:continue
        path=[]; cur=state
        while cur is not None and cur//2 not in used_records:
            path.append(cur); used_records.add(cur//2); cur=succ[cur]
        if not path:continue
        merged=oriented[path[0]]
        for s,t in zip(path,path[1:]): merged += oriented[t][edge_ov[(s,t)]:]
        merged=canonical(merged)
        outputs.append(merged); chains.append(path)
    # Cycles should have been excluded by DSU; emit any unvisited singleton defensively.
    for rid,s in enumerate(seqs):
        if rid not in used_records: outputs.append(s); chains.append([2*rid])
    outputs=sorted(set(outputs),key=lambda s:(-len(s),s))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('wb') as f:
        for i,s in enumerate(outputs,1):
            f.write(f'>stitched_{i:08d} len={len(s)}\n'.encode())
            for p in range(0,len(s),80):f.write(s[p:p+80]+b'\n')
    lengths=sorted((len(s) for s in outputs),reverse=True); total=sum(lengths); half=(total+1)//2; acc=0;n50=0
    for x in lengths:
        acc+=x
        if acc>=half:n50=x;break
    chain_lens=[len(p) for p in chains]
    print(f'inputs={n} outputs={len(outputs)} candidate_checks={checks} exact_matches={matches} reciprocal={len(reciprocal)} accepted={len(accepted)}')
    print(f'total_bp={total} n50={n50} largest={lengths[0] if lengths else 0} max_chain={max(chain_lens,default=0)} multi_record_chains={sum(x>1 for x in chain_lens)}')
if __name__=='__main__': main()
