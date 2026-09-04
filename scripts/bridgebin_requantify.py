#!/usr/bin/env python3
"""Recompute BridgeBin abundance after an external/refinement assignment step.

Matches the Rust quantifier: contig-length weighted median depth (robust_depth),
length-weighted mean depth, then sample-wise relative abundance from robust depths.
This helper is also the stable metric boundary for target-local Biological Brain tests.
"""

from __future__ import annotations

import argparse, csv, math
from collections import defaultdict
from pathlib import Path


def weighted_median(values):
    if not values:
        return 0.0
    values=sorted(values,key=lambda x:x[0])
    total=sum(w for _,w in values)
    if total<=0: return 0.0
    target=(total+1)//2; acc=0
    for value,weight in values:
        acc+=weight
        if acc>=target: return value
    return values[-1][0]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',type=Path,required=True)
    p.add_argument('--coverage',type=Path,required=True)
    p.add_argument('--truth-lengths',type=Path,required=True,help='TSV containing contig,length (benchmark)')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    lengths={r['contig']:int(r['length']) for r in csv.DictReader(a.truth_lengths.open(),delimiter='\t')}
    assignment={r['contig']:r['bin'] for r in csv.DictReader(a.assignments.open(),delimiter='\t')}
    with a.coverage.open() as h:
        reader=csv.DictReader(h,delimiter='\t'); samples=[x for x in reader.fieldnames if x!='contig']; cov={r['contig']:[float(r[s]) for s in samples] for r in reader}
    obs=defaultdict(lambda:[[] for _ in samples])
    for contig,bin_name in assignment.items():
        if bin_name in {'','.','NA','unbinned'} or contig not in cov or contig not in lengths: continue
        for i,depth in enumerate(cov[contig]): obs[bin_name][i].append((depth,lengths[contig]))
    robust={}; means={}; bp={}; counts={}
    for b,per_sample in obs.items():
        robust[b]=[]; means[b]=[]; bp[b]=[]; counts[b]=[]
        for values in per_sample:
            total_bp=sum(w for _,w in values); weighted=sum(v*w for v,w in values)
            robust[b].append(weighted_median(values)); means[b].append(weighted/total_bp if total_bp else 0.0); bp[b].append(total_bp); counts[b].append(len(values))
    totals=[sum(robust[b][i] for b in robust) for i in range(len(samples))]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['bin','sample','robust_depth','mean_depth','relative_abundance','covered_bp','covered_contigs'])
        for b in sorted(robust):
            for i,s in enumerate(samples):
                rel=robust[b][i]/totals[i] if totals[i]>0 else 0.0
                w.writerow([b,s,f'{robust[b][i]:.6f}',f'{means[b][i]:.6f}',f'{rel:.8f}',bp[b][i],counts[b][i]])
    print(f'bridgebin-requantify: bins={len(robust)} samples={len(samples)} output={a.output}')

if __name__=='__main__': main()
