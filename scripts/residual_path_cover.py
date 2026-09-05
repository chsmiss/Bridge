#!/usr/bin/env python3
"""Extract non-redundant residual branch patches from BridgeAsm GFA.

Shared backbone may carry residual capacity for several strain paths, but shared
sequence is never emitted repeatedly. Full residual paths are first extracted
from coverage + DR/GR/PE evidence, then reduced to sequence patches containing
k-mers absent from the current primary backbone. This preserves alternate
branches without inflating duplication ratio by copying the backbone.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path

_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s:str)->str: return s.translate(_COMP)[::-1]
def canonical(s:str)->str: return min(s,rc(s))
def ckmer(s:str)->str: return canonical(s)

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

def kmers(s:str,k:int):
    if len(s)<k: return
    for i in range(len(s)-k+1):
        q=s[i:i+k]
        if 'N' not in q: yield i,ckmer(q)

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
                tags=parse_tags(fields[3:]); seq[fields[1]]=fields[2].upper(); cov[fields[1]]=float(tags.get('KC','1'))
            elif fields[0]=='L' and fields[2]=='+' and fields[4]=='+':
                overlap=fields[5]
                if overlap.endswith('M'): k=int(overlap[:-1]) if k is None else k
                tags=parse_tags(fields[6:]); edges[(fields[1],fields[3])]=(int(tags.get('DR','0')),int(tags.get('GR','0')),int(tags.get('PE','0')))
    return seq,cov,edges,k or 0

def physical(ev):
    dr,gr,pe=ev; return max(dr,gr,pe)

def score(ev):
    dr,gr,pe=ev; return 100.0*dr + 40.0*gr + 20.0*pe

def choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen):
    valid=[]
    for edge,ev in cands:
        other=edge[0] if edge[1] in seen else edge[1]
        if other in seen or physical(ev)<min_support or edge_use[edge]>=capacity[edge]: continue
        valid.append((score(ev),edge,ev))
    if not valid: return None
    valid.sort(reverse=True,key=lambda x:(x[0],physical(x[2]),x[1])); best=valid[0]
    if len(valid)>1 and best[0]/max(total,1e-9)<dominance: return None
    return best[1]

def assemble(path,seq,k):
    if not path: return ''
    out=seq[path[0]]
    for n in path[1:]: out += seq[n][min(k,len(seq[n])):]
    return out

def extract_one(path:Path,secondary:float,dominance:float,min_support:int,max_copy:int):
    seq,cov,edges,k=parse_gfa(path); outgoing=defaultdict(list); incoming=defaultdict(list)
    for edge,ev in edges.items(): outgoing[edge[0]].append((edge,ev)); incoming[edge[1]].append((edge,ev))
    out_total={n:sum(score(ev) for _,ev in es) for n,es in outgoing.items()}; in_total={n:sum(score(ev) for _,ev in es) for n,es in incoming.items()}
    capacity={}
    for e,ev in edges.items():
        c=max(1,min(max_copy,int(round(min(cov.get(e[0],1.0),cov.get(e[1],1.0)))))); capacity[e]=c if physical(ev)>=min_support else 1
    edge_use=defaultdict(int); seeds=[]
    for e,ev in edges.items():
        if physical(ev)<min_support: continue
        so=score(ev); fo=so/max(out_total.get(e[0],so),1e-9); fi=so/max(in_total.get(e[1],so),1e-9)
        if len(outgoing[e[0]])<=1 and len(incoming[e[1]])<=1: continue
        if max(fo,fi)>=dominance and min(fo,fi)>=secondary: seeds.append((so*min(fo,fi),so,e,ev,fo,fi))
        elif fo>=secondary and fi>=secondary: seeds.append((so*min(fo,fi)*0.75,so,e,ev,fo,fi))
    seeds.sort(reverse=True,key=lambda x:(x[0],x[1],x[2])); paths=[]; seen_seq=set()
    for rank,so,seed,ev,fo,fi in seeds:
        if edge_use[seed]>=capacity[seed]: continue
        left,right=seed; seen={left,right}; prefix=[left]; suffix=[right]; current=left
        while True:
            cands=incoming.get(current,[]); total=sum(score(xev) for _,xev in cands if physical(xev)>=min_support)
            chosen=choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen)
            if chosen is None or chosen[1]!=current: break
            prev=chosen[0]; prefix.append(prev); seen.add(prev); current=prev
        prefix.reverse(); current=right
        while True:
            cands=outgoing.get(current,[]); total=sum(score(xev) for _,xev in cands if physical(xev)>=min_support)
            chosen=choose_extension(cands,total,dominance,min_support,edge_use,capacity,seen)
            if chosen is None or chosen[0]!=current: break
            nxt=chosen[1]; suffix.append(nxt); seen.add(nxt); current=nxt
        nodes=prefix+suffix; s=assemble(nodes,seq,k); can=canonical(s)
        if can in seen_seq: continue
        seen_seq.add(can)
        for aa,bb in zip(nodes,nodes[1:]):
            if (aa,bb) in edges: edge_use[(aa,bb)]+=1
        paths.append((rank,so,can,path.name,len(nodes),physical(ev),fo,fi,capacity[seed]))
    return paths

def novel_patches(seq:str,backbone:set[str],k:int,flank:int,max_novel_gap:int,min_novel_kmers:int,min_novel_fraction:float,min_length:int,max_patch_length:int):
    positions=[i for i,q in kmers(seq,k) if q not in backbone]
    if len(positions)<min_novel_kmers: return []
    clusters=[]; cur=[positions[0]]
    for p in positions[1:]:
        if p-cur[-1] <= max_novel_gap: cur.append(p)
        else: clusters.append(cur); cur=[p]
    clusters.append(cur); out=[]
    for cluster in clusters:
        if len(cluster)<min_novel_kmers: continue
        start=max(0,cluster[0]-flank); end=min(len(seq),cluster[-1]+k+flank)
        if end-start>max_patch_length:
            center=(cluster[0]+cluster[-1]+k)//2; half=max_patch_length//2; start=max(0,center-half); end=min(len(seq),start+max_patch_length); start=max(0,end-max_patch_length)
        patch=seq[start:end]; pkm=[q for _,q in kmers(patch,k)]; novel={q for q in pkm if q not in backbone}; frac=len(novel)/max(1,len(set(pkm)))
        if len(patch)>=min_length and len(novel)>=min_novel_kmers and frac>=min_novel_fraction: out.append((canonical(patch),len(novel),frac,start,end))
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('gfa',type=Path,nargs='+'); ap.add_argument('--backbone',type=Path,required=True); ap.add_argument('-o','--output',type=Path,required=True); ap.add_argument('--metadata',type=Path)
    ap.add_argument('--secondary-dominance',type=float,default=0.20); ap.add_argument('--extension-dominance',type=float,default=0.75); ap.add_argument('--min-support',type=int,default=3); ap.add_argument('--max-copy',type=int,default=3)
    ap.add_argument('--novel-k',type=int,default=31); ap.add_argument('--flank',type=int,default=120); ap.add_argument('--max-novel-gap',type=int,default=96); ap.add_argument('--min-novel-kmers',type=int,default=4); ap.add_argument('--min-novel-fraction',type=float,default=0.05); ap.add_argument('--min-length',type=int,default=200); ap.add_argument('--max-patch-length',type=int,default=1500); ap.add_argument('--max-patches',type=int,default=1500); ap.add_argument('--max-total-fraction',type=float,default=0.20)
    a=ap.parse_args(); backbone=set(); backbone_bases=0
    for _h,s in fasta_records(a.backbone): backbone_bases+=len(s); backbone.update(q for _,q in kmers(s,a.novel_k))
    raw=[]
    for g in a.gfa:
        if g.exists(): raw.extend(extract_one(g,a.secondary_dominance,a.extension_dominance,a.min_support,a.max_copy))
    candidates=[]
    for rank,so,s,gfa,nodes,sup,fo,fi,cap in raw:
        for patch,novel,frac,start,end in novel_patches(s,backbone,a.novel_k,a.flank,a.max_novel_gap,a.min_novel_kmers,a.min_novel_fraction,a.min_length,a.max_patch_length): candidates.append((sup*novel*frac,so,patch,gfa,nodes,sup,fo,fi,cap,novel,frac,start,end))
    candidates.sort(reverse=True,key=lambda x:(x[0],x[1],len(x[2]),x[2])); selected=[]; selected_novel=set(); seen=set(); total=0; max_total=max(a.min_length,int(backbone_bases*a.max_total_fraction))
    for item in candidates:
        patch=item[2]
        if patch in seen: continue
        novel_set={q for _,q in kmers(patch,a.novel_k) if q not in backbone}; new=novel_set-selected_novel
        if len(new)<a.min_novel_kmers: continue
        if len(selected)>=a.max_patches or total+len(patch)>max_total: continue
        seen.add(patch); selected_novel.update(novel_set); total+=len(patch); selected.append(item)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,item in enumerate(selected,1):
            s=item[2]; f.write(f'>residual_patch_{i:07d} len={len(s)} novel_kmers={item[9]} novel_fraction={item[10]:.4f}\n')
            for j in range(0,len(s),80): f.write(s[j:j+80]+'\n')
    if a.metadata:
        with a.metadata.open('w') as f:
            f.write('id\tgfa\tnodes\tlength\tseed_support\tsource_fraction\ttarget_fraction\tseed_capacity\tnovel_kmers\tnovel_fraction\tpatch_start\tpatch_end\n')
            for i,item in enumerate(selected,1): f.write('\t'.join(map(str,(i,item[3],item[4],len(item[2]),item[5],item[6],item[7],item[8],item[9],item[10],item[11],item[12])))+'\n')
    print(f'raw_residual_paths\t{len(raw)}'); print(f'candidate_patches\t{len(candidates)}'); print(f'residual_patches\t{len(selected)}'); print(f'residual_bases\t{total}'); print(f'novel_kmers_retained\t{len(selected_novel)}'); print(f'backbone_bases\t{backbone_bases}')
if __name__=='__main__': main()
