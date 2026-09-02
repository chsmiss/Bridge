#!/usr/bin/env python3
"""Extract evidence-supported alternate paths from BridgeAsm GFA.

The extractor deliberately allows high-coverage backbone edges/nodes to be reused
by multiple paths. Alternate branch seeds consume residual edge capacity, while
extension through another ambiguous branch requires a dominant physical signal.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path

_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s:str)->str: return s.translate(_COMP)[::-1]
def canonical(s:str)->str: return min(s,rc(s))

def parse_tags(fields):
    out={}
    for item in fields:
        p=item.split(':',2)
        if len(p)==3: out[p[0]]=p[2]
    return out

def parse_gfa(path:Path):
    seq={}; cov={}; edges={}; k=None
    with path.open() as f:
        for raw in f:
            fields=raw.rstrip('\n').split('\t')
            if not fields: continue
            if fields[0]=='S':
                tags=parse_tags(fields[3:])
                seq[fields[1]]=fields[2].upper()
                cov[fields[1]]=float(tags.get('KC','1'))
            elif fields[0]=='L' and fields[2]=='+' and fields[4]=='+':
                overlap=fields[5]
                if overlap.endswith('M'):
                    ov=int(overlap[:-1]); k=ov if k is None else k
                tags=parse_tags(fields[6:])
                edges[(fields[1],fields[3])]=(int(tags.get('DR','0')),int(tags.get('GR','0')),int(tags.get('PE','0')))
    return seq,cov,edges,k or 0

def physical(ev):
    dr,gr,pe=ev
    return max(dr,gr,pe)

def score(ev):
    dr,gr,pe=ev
    return 100.0*dr + 40.0*gr + 20.0*pe

def choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen):
    valid=[]
    for edge,ev in cands:
        other=edge[0] if edge[1] in seen else edge[1]
        if other in seen or physical(ev)<min_support or edge_use[edge]>=capacity[edge]: continue
        valid.append((score(ev),edge,ev))
    if not valid: return None
    valid.sort(reverse=True,key=lambda x:(x[0],physical(x[2]),x[1]))
    best=valid[0]
    if len(valid)>1 and best[0]/max(total,1e-9)<dominance: return None
    return best[1]

def assemble(path,seq,k):
    if not path: return ''
    out=seq[path[0]]
    for n in path[1:]: out += seq[n][min(k,len(seq[n])):]
    return out

def extract_one(path:Path, secondary:float, dominance:float, min_support:int, min_length:int, max_copy:int):
    seq,cov,edges,k=parse_gfa(path)
    outgoing=defaultdict(list); incoming=defaultdict(list)
    for edge,ev in edges.items():
        outgoing[edge[0]].append((edge,ev)); incoming[edge[1]].append((edge,ev))
    out_total={n:sum(score(ev) for _,ev in es) for n,es in outgoing.items()}
    in_total={n:sum(score(ev) for _,ev in es) for n,es in incoming.items()}
    capacity={}
    for e,ev in edges.items():
        c=max(1,min(max_copy,int(round(min(cov.get(e[0],1.0),cov.get(e[1],1.0))))))
        capacity[e]=c if physical(ev)>=min_support else 1
    edge_use=defaultdict(int)
    seeds=[]
    for e,ev in edges.items():
        if physical(ev)<min_support: continue
        so=score(ev); fo=so/max(out_total.get(e[0],so),1e-9); fi=so/max(in_total.get(e[1],so),1e-9)
        branching=len(outgoing[e[0]])>1 or len(incoming[e[1]])>1
        if not branching: continue
        if max(fo,fi)>=dominance and min(fo,fi)>=secondary:
            seeds.append((so*min(fo,fi),so,e,ev,fo,fi))
        elif fo>=secondary and fi>=secondary:
            seeds.append((so*min(fo,fi)*0.75,so,e,ev,fo,fi))
    seeds.sort(reverse=True,key=lambda x:(x[0],x[1],x[2]))
    emitted=[]; seen_seq=set(); meta=[]
    for _rank,_so,seed,ev,fo,fi in seeds:
        if edge_use[seed]>=capacity[seed]: continue
        left,right=seed; seen={left,right}; prefix=[left]; suffix=[right]
        current=left
        while True:
            cands=incoming.get(current,[])
            total=sum(score(xev) for _,xev in cands if physical(xev)>=min_support)
            chosen=choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen)
            if chosen is None or chosen[1]!=current: break
            prev=chosen[0]; prefix.append(prev); seen.add(prev); current=prev
        prefix.reverse()
        current=right
        while True:
            cands=outgoing.get(current,[])
            total=sum(score(xev) for _,xev in cands if physical(xev)>=min_support)
            chosen=choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen)
            if chosen is None or chosen[0]!=current: break
            nxt=chosen[1]; suffix.append(nxt); seen.add(nxt); current=nxt
        nodes=prefix+suffix
        s=assemble(nodes,seq,k)
        if len(s)<min_length: continue
        can=canonical(s)
        if can in seen_seq: continue
        seen_seq.add(can); emitted.append(can)
        for aa,bb in zip(nodes,nodes[1:]):
            if (aa,bb) in edges: edge_use[(aa,bb)] += 1
        meta.append((len(emitted),path.name,len(nodes),len(can),physical(ev),fo,fi,capacity[seed]))
    return emitted,meta

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('gfa',type=Path,nargs='+')
    ap.add_argument('-o','--output',type=Path,required=True)
    ap.add_argument('--metadata',type=Path)
    ap.add_argument('--secondary-dominance',type=float,default=0.18)
    ap.add_argument('--extension-dominance',type=float,default=0.72)
    ap.add_argument('--min-support',type=int,default=2)
    ap.add_argument('--min-length',type=int,default=300)
    ap.add_argument('--max-copy',type=int,default=6)
    a=ap.parse_args()
    allseq=[]; allmeta=[]; seen=set()
    for g in a.gfa:
        if not g.exists(): continue
        seqs,meta=extract_one(g,a.secondary_dominance,a.extension_dominance,a.min_support,a.min_length,a.max_copy)
        for s,m in zip(seqs,meta):
            c=canonical(s)
            if c not in seen: seen.add(c); allseq.append(c); allmeta.append(m)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,s in enumerate(allseq,1):
            f.write(f'>residual_path_{i:07d} len={len(s)}\n')
            for j in range(0,len(s),80): f.write(s[j:j+80]+'\n')
    if a.metadata:
        with a.metadata.open('w') as f:
            f.write('id\tgfa\tnodes\tlength\tseed_support\tsource_fraction\ttarget_fraction\tseed_capacity\n')
            for row in allmeta: f.write('\t'.join(map(str,row))+'\n')
    print(f'residual_paths\t{len(allseq)}')
    print(f'residual_bases\t{sum(map(len,allseq))}')
if __name__=='__main__': main()
