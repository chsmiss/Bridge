#!/usr/bin/env python3
"""Project flanked strain paths onto a backbone and retain novel sequence evidence."""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s):return s.translate(_COMP)[::-1]
def canon(s):return min(s,rc(s))
def records(path):
    h=None;c=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                if h is not None:yield h,''.join(c).upper()
                h=line[1:];c=[]
            else:c.append(line)
    if h is not None:yield h,''.join(c).upper()
def kmers(s,k,stride=1):
    for i in range(0,len(s)-k+1,stride):
        q=s[i:i+k]
        if 'N' not in q:yield canon(q)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('backbone',type=Path);ap.add_argument('candidates',type=Path,nargs='+');ap.add_argument('-o','--output',type=Path,required=True);ap.add_argument('--map',type=Path,required=True);ap.add_argument('--k',type=int,default=31);ap.add_argument('--projection-k',type=int,default=21);ap.add_argument('--min-novel-fraction',type=float,default=0.03);ap.add_argument('--min-projection-fraction',type=float,default=0.10);a=ap.parse_args()
    backs=list(records(a.backbone)); global_k=set(); index=defaultdict(set)
    for bid,(_h,s) in enumerate(backs):
        global_k.update(kmers(s,a.k))
        for q in set(kmers(s,a.projection_k,3)):index[q].add(bid)
    selected=[]; rows=[]; seen=set()
    for fp in a.candidates:
        if not fp.exists():continue
        for h,s in records(fp):
            c=canon(s)
            if c in seen:continue
            seen.add(c)
            ks=list(kmers(s,a.k)); novel=sum(q not in global_k for q in ks); novel_frac=novel/max(1,len(ks))
            counts=defaultdict(int); pks=list(kmers(s,a.projection_k,3))
            for q in pks:
                for bid in index.get(q,()):counts[bid]+=1
            if counts:
                bid,hits=max(counts.items(),key=lambda x:(x[1],-x[0])); proj_frac=hits/max(1,len(pks)); target=backs[bid][0]
            else:proj_frac=0.0;target='.'
            keep=novel_frac>=a.min_novel_fraction and proj_frac>=a.min_projection_fraction
            rows.append((h,len(s),target,proj_frac,novel_frac,novel,keep))
            if keep:selected.append((h,c))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,(h,s) in enumerate(selected,1):
            f.write(f'>projected_strain_{i:07d} len={len(s)} source={h}\n')
            for j in range(0,len(s),80):f.write(s[j:j+80]+'\n')
    with a.map.open('w') as f:
        f.write('candidate\tlength\tbackbone_target\tprojection_fraction\tnovel_fraction\tnovel_kmers\tselected\n')
        for row in rows:f.write('\t'.join(map(str,row))+'\n')
    print(f'projected_candidates\t{len(rows)}');print(f'selected_strain_paths\t{len(selected)}');print(f'selected_bases\t{sum(len(s) for _,s in selected)}')
if __name__=='__main__':main()
