#!/usr/bin/env python3
"""Project flanked strain/bulge paths onto a primary backbone.

A candidate is retained as evidence when both sequence flanks anchor to the
same backbone contig and it carries k-mers absent from the backbone. This is
closer to bulge projection than whole-sequence similarity: the variable middle
is expected to differ, while the physical flanks identify where it belongs.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(_COMP)[::-1]
def canon(s): return min(s,rc(s))
def records(path):
    h=None;c=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: yield h,''.join(c).upper()
                h=line[1:];c=[]
            else:c.append(line)
    if h is not None: yield h,''.join(c).upper()
def kmers(s,k,stride=1):
    if len(s)<k:return
    for i in range(0,len(s)-k+1,stride):
        q=s[i:i+k]
        if 'N' not in q: yield canon(q)
def flank_hits(segment,index,k,stride):
    counts=defaultdict(int); seen=set()
    for q in kmers(segment,k,stride):
        if q in seen: continue
        seen.add(q)
        for bid in index.get(q,()): counts[bid]+=1
    return counts,len(seen)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('backbone',type=Path);ap.add_argument('candidates',type=Path,nargs='+');ap.add_argument('-o','--output',type=Path,required=True);ap.add_argument('--map',type=Path,required=True)
    ap.add_argument('--k',type=int,default=31);ap.add_argument('--projection-k',type=int,default=21);ap.add_argument('--flank-length',type=int,default=90);ap.add_argument('--projection-stride',type=int,default=2);ap.add_argument('--min-flank-hits',type=int,default=2);ap.add_argument('--min-novel-kmers',type=int,default=3);ap.add_argument('--min-novel-fraction',type=float,default=0.01);a=ap.parse_args()
    backs=list(records(a.backbone)); global_k=set(); index=defaultdict(set)
    for bid,(_h,s) in enumerate(backs):
        global_k.update(kmers(s,a.k,1))
        for q in set(kmers(s,a.projection_k,a.projection_stride)): index[q].add(bid)
    selected=[]; rows=[]; seen=set()
    for fp in a.candidates:
        if not fp.exists(): continue
        for h,s in records(fp):
            c=canon(s)
            if c in seen or len(s)<max(a.k,a.projection_k): continue
            seen.add(c); ks=list(kmers(s,a.k,1)); novel_set={q for q in ks if q not in global_k}; novel=len(novel_set); novel_frac=novel/max(1,len(set(ks)))
            flank=min(a.flank_length,max(a.projection_k,(len(s)-a.projection_k)//3)); left=s[:flank]; right=s[-flank:]
            lh,ln=flank_hits(left,index,a.projection_k,a.projection_stride); rh,rn=flank_hits(right,index,a.projection_k,a.projection_stride); common=set(lh)&set(rh)
            if common:
                bid=max(common,key=lambda b:(min(lh[b],rh[b]),lh[b]+rh[b],-b)); left_hits=lh[bid]; right_hits=rh[bid]; target=backs[bid][0]
            else: bid=None; left_hits=max(lh.values(),default=0); right_hits=max(rh.values(),default=0); target='.'
            keep=bid is not None and left_hits>=a.min_flank_hits and right_hits>=a.min_flank_hits and novel>=a.min_novel_kmers and novel_frac>=a.min_novel_fraction
            rows.append((h,len(s),target,left_hits,right_hits,ln,rn,novel_frac,novel,keep))
            if keep:selected.append((h,c))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,(h,s) in enumerate(selected,1):
            f.write(f'>projected_strain_{i:07d} len={len(s)} source={h}\n')
            for j in range(0,len(s),80):f.write(s[j:j+80]+'\n')
    with a.map.open('w') as f:
        f.write('candidate\tlength\tbackbone_target\tleft_hits\tright_hits\tleft_flank_kmers\tright_flank_kmers\tnovel_fraction\tnovel_kmers\tselected\n')
        for row in rows:f.write('\t'.join(map(str,row))+'\n')
    print(f'projected_candidates\t{len(rows)}');print(f'selected_strain_paths\t{len(selected)}');print(f'selected_bases\t{sum(len(s) for _,s in selected)}')
if __name__=='__main__':main()
