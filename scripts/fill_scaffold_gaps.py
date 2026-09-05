#!/usr/bin/env python3
"""Conservatively fill N-gaps with local read-supported de Bruijn paths."""
from __future__ import annotations
import argparse, gzip, re
from collections import defaultdict
from pathlib import Path

_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s:str)->str: return s.translate(_COMP)[::-1]
def canon(s:str)->str: return min(s,rc(s))

def fasta(path:Path):
    out=[]; h=None; c=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: out.append((h,''.join(c).upper()))
                h=line[1:].split()[0]; c=[]
            else:c.append(line)
    if h is not None: out.append((h,''.join(c).upper()))
    return out

def fastq_pairs(p1:Path,p2:Path):
    def reader(p):
        op=gzip.open if p.suffix=='.gz' else open
        with op(p,'rt') as f:
            while True:
                h=f.readline()
                if not h: return
                s=f.readline().strip(); f.readline(); f.readline()
                yield h.strip().split()[0].lstrip('@'),s.upper()
    for aa,bb in zip(reader(p1),reader(p2)):
        yield aa,bb

def anchor_kmers(s,k,stride=4):
    vals=[]
    if len(s)<k:return vals
    for i in range(0,len(s)-k+1,stride):
        q=s[i:i+k]
        if 'N' not in q: vals.append(canon(q))
    if (len(s)-k)%stride:
        q=s[-k:]
        if 'N' not in q: vals.append(canon(q))
    return vals

def local_path(left,right,reads,k,expected,max_steps,dominance):
    seqs=[left,right]
    for s in reads: seqs.extend((s,rc(s)))
    out=defaultdict(lambda:defaultdict(int))
    for idx,s in enumerate(seqs):
        weight=3 if idx<2 else 1
        for i in range(len(s)-k):
            aa=s[i:i+k]; bb=s[i+1:i+k+1]
            if 'N' in aa or 'N' in bb: continue
            out[aa][bb]+=weight
    start=left[-k:]; goal=right[:k]
    if start not in out:return None
    path=[start]; cur=start; seen={start}; max_walk=min(max_steps,max(k+1,expected+2*k+250))
    for _ in range(max_walk):
        if cur==goal: break
        choices=[(cnt,nxt) for nxt,cnt in out.get(cur,{}).items() if nxt not in seen or nxt==goal]
        if not choices:return None
        choices.sort(reverse=True)
        if len(choices)>1 and choices[0][0]/sum(c for c,_ in choices)<dominance:return None
        _,nxt=choices[0]; path.append(nxt); cur=nxt
        if cur!=goal: seen.add(cur)
    if cur!=goal:return None
    assembled=path[0]+''.join(x[-1] for x in path[1:])
    if len(assembled)<2*k:return None
    fill=assembled[k:-k]
    tolerance=max(80,int(expected*0.75)+30)
    if abs(len(fill)-expected)>tolerance:return None
    return fill

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('scaffolds',type=Path); ap.add_argument('-1','--read1',type=Path,required=True); ap.add_argument('-2','--read2',type=Path,required=True)
    ap.add_argument('-o','--output',type=Path,required=True); ap.add_argument('--report',type=Path,required=True); ap.add_argument('--anchor-k',type=int,default=31); ap.add_argument('--local-k',type=int,default=21); ap.add_argument('--flank',type=int,default=180); ap.add_argument('--dominance',type=float,default=0.65); ap.add_argument('--max-reads-per-gap',type=int,default=1000)
    a=ap.parse_args(); recs=fasta(a.scaffolds); gaps=[]; anchor_index=defaultdict(list)
    for rid,(h,s) in enumerate(recs):
        for m in re.finditer(r'N+',s):
            ll=max(0,m.start()-a.flank); rr=min(len(s),m.end()+a.flank); left=s[ll:m.start()]; right=s[m.end():rr]
            if len(left)<a.anchor_k or len(right)<a.anchor_k: continue
            gid=len(gaps); gaps.append({'rid':rid,'start':m.start(),'end':m.end(),'left':left,'right':right,'reads':[]})
            for q in anchor_kmers(left[-min(len(left),100):],a.anchor_k): anchor_index[q].append(gid)
            for q in anchor_kmers(right[:min(len(right),100)],a.anchor_k): anchor_index[q].append(gid)
    for (_id1,s1),(_id2,s2) in fastq_pairs(a.read1,a.read2):
        hit=set()
        for s in (s1,s2):
            for i in range(len(s)-a.anchor_k+1):
                q=canon(s[i:i+a.anchor_k])
                for gid in anchor_index.get(q,()): hit.add(gid)
        for gid in hit:
            buf=gaps[gid]['reads']
            if len(buf)<a.max_reads_per_gap: buf.extend((s1,s2))
    replacements=defaultdict(list); report=[]
    for gid,g in enumerate(gaps):
        expected=g['end']-g['start']
        fill=local_path(g['left'],g['right'],g['reads'],a.local_k,expected,expected+500,a.dominance)
        if fill is not None:
            replacements[g['rid']].append((g['start'],g['end'],fill)); status='filled'
        else: status='unresolved'
        report.append((gid,recs[g['rid']][0],g['start'],g['end'],expected,len(g['reads']),status,0 if fill is None else len(fill)))
    output=[]
    for rid,(h,s) in enumerate(recs):
        for st,en,fill in sorted(replacements[rid],reverse=True): s=s[:st]+fill+s[en:]
        output.append((h,s))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,(_h,s) in enumerate(output,1):
            f.write(f'>gap_refined_{i:07d} len={len(s)}\n')
            for j in range(0,len(s),80):f.write(s[j:j+80]+'\n')
    with a.report.open('w') as f:
        f.write('gap_id\tscaffold\tstart\tend\testimated_gap\tlocal_reads\tstatus\tfilled_bases\n')
        for row in report:f.write('\t'.join(map(str,row))+'\n')
    print(f'gaps\t{len(gaps)}'); print(f'filled\t{sum(r[6]=="filled" for r in report)}')
if __name__=='__main__':main()
