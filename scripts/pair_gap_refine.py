#!/usr/bin/env python3
"""Paired-end contig linker with empirical insert-size estimation.

Maps reads to preliminary contigs with minimap2, estimates the library insert
size from unique same-contig FR pairs, then joins reciprocal-dominant contig
ends. Negative/zero estimated gaps are closed only by exact suffix/prefix
overlap; positive gaps are represented with Ns for later local assembly.
"""
from __future__ import annotations
import argparse, re, statistics, subprocess
from collections import defaultdict
from pathlib import Path

_COMP=str.maketrans('ACGTN','TGCAN')
def rc(s:str)->str: return s.translate(_COMP)[::-1]
def canonical(s:str)->str: return min(s,rc(s))

def fasta(path:Path):
    out=[]; h=None; chunks=[]
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: out.append((h,''.join(chunks).upper()))
                h=line[1:].split()[0]; chunks=[]
            else: chunks.append(line)
    if h is not None: out.append((h,''.join(chunks).upper()))
    return out

def ref_span(cigar:str)->int:
    if cigar=='*': return 0
    return sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar) if op in 'MDN=X')

def parse_rec(fields):
    flag=int(fields[1])
    if flag & (0x4|0x100|0x800): return None
    start=int(fields[3])-1; span=ref_span(fields[5])
    return {'q':fields[0],'flag':flag,'r':fields[2],'start':start,'end':start+span,
            'mapq':int(fields[4]),'rev':bool(flag&16),'tlen':int(fields[8]),'qlen':0 if fields[9]=='*' else len(fields[9])}

def run_mapping(contigs:Path,r1:Path,r2:Path,threads:int):
    cmd=['minimap2','-ax','sr','--secondary=no','-t',str(threads),str(contigs),str(r1),str(r2)]
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True,bufsize=1)
    pending={}
    assert proc.stdout is not None
    for line in proc.stdout:
        if line.startswith('@'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<11: continue
        rec=parse_rec(f)
        if rec is None: continue
        old=pending.pop(rec['q'],None)
        if old is None: pending[rec['q']]=rec
        else: yield old,rec
    code=proc.wait()
    if code: raise subprocess.CalledProcessError(code,cmd)

def robust_insert(vals):
    if not vals: return 300.0,50.0,0
    med=float(statistics.median(vals)); dev=[abs(x-med) for x in vals]
    mad=float(statistics.median(dev)) if dev else 0.0
    return med,max(mad,10.0),len(vals)

def source_roles(rec,length,end_window):
    roles=[]
    if not rec['rev'] and length-rec['end']<=end_window:
        roles.append(('+',float(length-rec['start']),'R'))
    if rec['rev'] and rec['start']<=end_window:
        roles.append(('-',float(rec['end']),'L'))
    return roles

def target_roles(rec,length,end_window):
    roles=[]
    if rec['rev'] and rec['start']<=end_window:
        roles.append(('+',float(rec['end']),'L'))
    if not rec['rev'] and length-rec['end']<=end_window:
        roles.append(('-',float(length-rec['start']),'R'))
    return roles

def canonical_end_edge(s,sside,t,tside):
    aa=(s,sside); bb=(t,tside)
    return (aa,bb) if aa<=bb else (bb,aa)

def longest_overlap(a,b,min_overlap,max_overlap=None):
    limit=min(len(a),len(b),max_overlap or max(len(a),len(b)))
    for n in range(limit,min_overlap-1,-1):
        if a[-n:]==b[:n]: return n
    return 0

def orient(seq,linked_left_side):
    return seq if linked_left_side=='L' else rc(seq)

def build_scaffolds(records,selected,min_overlap,max_gap):
    seq={h:s for h,s in records}; nbr={}; edge_data={}
    for (aa,bb),data in selected.items():
        nbr[aa]=bb; nbr[bb]=aa; edge_data[frozenset((aa,bb))]=data
    used=set(); outputs=[]; links=[]
    for cid in seq:
        if cid in used: continue
        left=(cid,'L'); right=(cid,'R')
        if left in nbr and right in nbr: continue
        free_side='L' if left not in nbr else 'R'
        current_seq=seq[cid] if free_side=='L' else rc(seq[cid]); used.add(cid)
        outgoing=(cid,'R' if free_side=='L' else 'L')
        while outgoing in nbr:
            other=nbr[outgoing]; oid,oside=other
            if oid in used: break
            nxt=orient(seq[oid],oside); data=edge_data[frozenset((outgoing,other))]
            gap=int(round(data['gap'])); overlap=0
            if gap<=0:
                overlap=longest_overlap(current_seq,nxt,min_overlap,min(len(current_seq),len(nxt),max(4*min_overlap,abs(gap)+150)))
            if overlap:
                current_seq += nxt[overlap:]; emitted_gap=0
            else:
                emitted_gap=max(1,min(max_gap,max(0,gap))); current_seq += 'N'*emitted_gap + nxt
            links.append((outgoing[0],outgoing[1],oid,oside,data['support'],data['gap'],emitted_gap,overlap))
            used.add(oid); outgoing=(oid,'R' if oside=='L' else 'L')
        outputs.append(canonical(current_seq))
    for cid,s in seq.items():
        if cid not in used: outputs.append(canonical(s)); used.add(cid)
    uniq=[]; seen=set()
    for s in sorted(outputs,key=lambda x:(-len(x),x)):
        c=canonical(s)
        if c not in seen: seen.add(c); uniq.append(c)
    return uniq,links

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('contigs',type=Path); ap.add_argument('-1','--read1',type=Path,required=True); ap.add_argument('-2','--read2',type=Path,required=True)
    ap.add_argument('-o','--output',type=Path,required=True); ap.add_argument('--links',type=Path,required=True)
    ap.add_argument('--threads',type=int,default=2); ap.add_argument('--min-mapq',type=int,default=20); ap.add_argument('--min-support',type=int,default=3)
    ap.add_argument('--dominance',type=float,default=0.75); ap.add_argument('--end-window',type=int,default=500); ap.add_argument('--min-overlap',type=int,default=31); ap.add_argument('--max-gap',type=int,default=1000)
    a=ap.parse_args(); records=fasta(a.contigs); lengths={h:len(s) for h,s in records}
    insert_values=[]
    for x,y in run_mapping(a.contigs,a.read1,a.read2,a.threads):
        if x['mapq']<a.min_mapq or y['mapq']<a.min_mapq or x['r']!=y['r'] or x['rev']==y['rev']:
            continue
        t=abs(x['tlen']) or abs(y['tlen'])
        if 100<=t<=2000: insert_values.append(t)
    insert,mad,n_insert=robust_insert(insert_values)
    support=defaultdict(set); gaps=defaultdict(list)
    for x,y in run_mapping(a.contigs,a.read1,a.read2,a.threads):
        if x['mapq']<a.min_mapq or y['mapq']<a.min_mapq or x['r']==y['r'] or x['r'] not in lengths or y['r'] not in lengths: continue
        generated=set()
        for src,tgt in ((x,y),(y,x)):
            for _ori,outer,side_s in source_roles(src,lengths[src['r']],a.end_window):
                for _tori,touter,side_t in target_roles(tgt,lengths[tgt['r']],a.end_window):
                    if src['r']==tgt['r']: continue
                    edge=canonical_end_edge(src['r'],side_s,tgt['r'],side_t)
                    if edge in generated: continue
                    generated.add(edge); gap=insert-outer-touter
                    if gap < -1000 or gap > a.max_gap*2: continue
                    support[edge].add(src['q']); gaps[edge].append(gap)
    endpoint_total=defaultdict(int); candidates={}
    for edge,names in support.items():
        n=len(names)
        if n<a.min_support: continue
        med=float(statistics.median(gaps[edge])); spread=float(statistics.median([abs(g-med) for g in gaps[edge]])) if gaps[edge] else 0.0
        if spread>max(3*mad,80.0): continue
        candidates[edge]={'support':n,'gap':med,'spread':spread}; endpoint_total[edge[0]] += n; endpoint_total[edge[1]] += n
    selected={}; best={}
    for edge,data in candidates.items():
        for ep in edge:
            key=(data['support'], -abs(data['gap']), tuple(edge))
            if ep not in best or key>best[ep][0]: best[ep]=(key,edge)
    for edge,data in candidates.items():
        if best.get(edge[0],(None,None))[1]!=edge or best.get(edge[1],(None,None))[1]!=edge: continue
        if data['support']/max(1,endpoint_total[edge[0]])<a.dominance or data['support']/max(1,endpoint_total[edge[1]])<a.dominance: continue
        selected[edge]=data
    seqs,links=build_scaffolds(records,selected,a.min_overlap,a.max_gap)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w') as f:
        for i,s in enumerate(seqs,1):
            f.write(f'>pair_refined_{i:07d} len={len(s)}\n')
            for j in range(0,len(s),80): f.write(s[j:j+80]+'\n')
    with a.links.open('w') as f:
        f.write(f'# insert_median={insert:.3f}\tinsert_mad={mad:.3f}\tinsert_pairs={n_insert}\n')
        f.write('source\tsource_side\ttarget\ttarget_side\tsupport\testimated_gap\temitted_gap\toverlap\n')
        for row in links: f.write('\t'.join(map(str,row))+'\n')
    print(f'insert_median\t{insert:.3f}'); print(f'insert_mad\t{mad:.3f}'); print(f'insert_pairs\t{n_insert}')
    print(f'candidate_links\t{len(candidates)}'); print(f'selected_links\t{len(selected)}'); print(f'output_records\t{len(seqs)}')
if __name__=='__main__': main()
