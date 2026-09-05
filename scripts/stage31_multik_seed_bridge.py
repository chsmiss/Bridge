#!/usr/bin/env python3
"""Stage31: multi-k evidence-gated seed-to-seed bridging for 200k Zymo.

The Stage24 gap-free contigs are immutable. k31-aggressive is discovery only;
independent k55/k77/k99 assemblies (same seed + raw reads, mercy disabled)
validate topology. Only exact seed-end bridges that win at both physical ends
are used, cycles are forbidden, and no N gaps are emitted. A final conservative
k31 one-ended extension is applied after bridge topology is fixed.
"""
from __future__ import annotations
import argparse, json, sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator
import bridge_carry_forward_pipeline as carry
import stage26_stage24_carry_forward as s26

TR=bytes.maketrans(b'ACGTN',b'TGCAN')
W={31:1,55:2,77:3,99:4}
VARIANTS={
 'consensus2':(2,77,2), 'k99_single':(1,99,2), 'k77_single':(1,77,1),
 'k55_single':(1,55,1), 'k31_only':(0,31,1),
}

def rc(s:bytes)->bytes:return s.translate(TR)[::-1]
def fasta(p:Path)->Iterator[tuple[str,bytes]]:
    n=None; c=[]
    with p.open('rb') as h:
        for raw in h:
            x=raw.strip()
            if not x: continue
            if x.startswith(b'>'):
                if n is not None: yield n,b''.join(c).upper()
                n=x[1:].decode(errors='replace'); c=[]
            else:c.append(x)
    if n is not None: yield n,b''.join(c).upper()

def write_fa(p:Path,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('wb') as h:
        for n,s in rows:
            h.write(b'>'+n.encode()+b'\n')
            for i in range(0,len(s),80): h.write(s[i:i+80]+b'\n')

def n50(rows):
    a=sorted((len(s) for s in rows),reverse=True); half=(sum(a)+1)//2; z=0
    for x in a:
        z+=x
        if z>=half:return x
    return 0

def enc31(s:bytes):
    if len(s)!=31:return None
    z=0
    for b in s:
        if b==65:v=0
        elif b==67:v=1
        elif b==71:v=2
        elif b==84:v=3
        else:return None
        z=(z<<2)|v
    return z

def roll31(s:bytes):
    mask=(1<<62)-1; z=0; valid=0
    for i,b in enumerate(s):
        if b==65:v=0
        elif b==67:v=1
        elif b==71:v=2
        elif b==84:v=3
        else:z=0;valid=0;continue
        z=((z<<2)|v)&mask;valid+=1
        if valid>=31:yield i-30,z

def endpoint(state:int,source:bool):
    sid=state//2; rev=state%2
    return (sid,('L' if rev else 'R') if source else ('R' if rev else 'L'))

@dataclass(frozen=True)
class P:
    k:int; ls:int; rs:int; lo:int; ro:int; mid:bytes; cand:bytes; name:str
    @property
    def le(self):return endpoint(self.ls,True)
    @property
    def re(self):return endpoint(self.rs,False)
    @property
    def pair(self):return tuple(sorted((self.le,self.re)))
    @property
    def midkey(self):return min(self.mid,rc(self.mid))
    def flip(self):return P(self.k,self.rs^1,self.ls^1,self.ro,self.lo,rc(self.mid),rc(self.cand),self.name+'/rc')

def best(items,margin):
    if not items:return None
    a=sorted(items,key=lambda x:(-x[1],x[0]))
    if len(a)>1 and (a[0][1]==a[1][1] or a[0][1]-a[1][1]<margin):return None
    return a[0]

def discover(seedrows,candrows,k,minov=120,margin=30):
    seeds=[]
    for _,s in seedrows: seeds += [s,rc(s)]
    cands=[]
    for j,(n,s) in enumerate(candrows): cands += [(f'{n}|{j}|+',s),(f'{n}|{j}|-',rc(s))]
    cp=defaultdict(list)
    for j,(_,s) in enumerate(cands):
        if len(s)>=minov:
            q=enc31(s[:31])
            if q is not None: cp[q].append(j)
    left=defaultdict(list)
    for si,s in enumerate(seeds):
        for pos,q in roll31(s):
            ov=len(s)-pos
            if ov<minov:break
            for cj in cp.get(q,[]):
                c=cands[cj][1]
                if ov<=len(c) and c.startswith(s[pos:]): left[cj].append((si,ov))
    sp=defaultdict(list)
    for si,s in enumerate(seeds):
        if len(s)>=minov:
            q=enc31(s[:31])
            if q is not None: sp[q].append(si)
    right=defaultdict(list)
    for cj,(_,c) in enumerate(cands):
        for pos,q in roll31(c):
            ov=len(c)-pos
            if ov<minov:break
            for si in sp.get(q,[]):
                s=seeds[si]
                if ov<=len(s) and s.startswith(c[pos:]): right[cj].append((si,ov))
    out=[]; seen=set(); both=amb=same=short=0
    for cj,(name,c) in enumerate(cands):
        if not left.get(cj) or not right.get(cj):continue
        both+=1; L=best(left[cj],margin); R=best(right[cj],margin)
        if L is None or R is None:amb+=1;continue
        ls,lo=L;rs,ro=R
        if ls//2==rs//2:same+=1;continue
        m=len(c)-lo-ro
        if m<1:short+=1;continue
        p=P(k,ls,rs,lo,ro,c[lo:len(c)-ro],c,name); key=(p.pair,k,p.midkey)
        if key not in seen:seen.add(key);out.append(p)
    return out,{'records':len(candrows),'both_end':both,'ambiguous':amb,'same_seed':same,'short_internal':short,'proposals':len(out)}

def select(proposals,min_high,min_k,end_margin,nseed):
    groups=defaultdict(lambda:defaultdict(list))
    for p in proposals: groups[p.pair][p.midkey].append(p)
    edges=[]; ambiguous=0
    for pair,alts in groups.items():
        ranked=[]
        for mid,g in alts.items():
            ks=tuple(sorted({p.k for p in g}));score=sum(W[x] for x in ks)
            ranked.append((score,mid,g,ks))
        ranked.sort(key=lambda x:(-x[0],x[1]))
        if len(ranked)>1 and ranked[0][0]<=ranked[1][0]:ambiguous+=1;continue
        score,_,g,ks=ranked[0]
        if sum(x>=55 for x in ks)<min_high or max(ks)<min_k:continue
        p=max(g,key=lambda x:(x.k,x.lo+x.ro,len(x.cand)))
        edges.append((p,ks,score))
    inc=defaultdict(list)
    for e in edges:
        p=e[0];inc[p.le].append(e);inc[p.re].append(e)
    win={}
    for ep,a in inc.items():
        a=sorted(a,key=lambda e:(-e[2],-max(e[1]),-(e[0].lo+e[0].ro),e[0].pair))
        win[ep]=None if len(a)>1 and a[0][2]-a[1][2]<end_margin else a[0]
    rec=[e for e in edges if win.get(e[0].le) is e and win.get(e[0].re) is e]
    parent=list(range(nseed))
    def f(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    chosen=[];cycles=0
    for e in sorted(rec,key=lambda e:(-e[2],-max(e[1]),e[0].pair)):
        a,b=e[0].le[0],e[0].re[0];ra,rb=f(a),f(b)
        if ra==rb:cycles+=1;continue
        parent[rb]=ra;chosen.append(e)
    return chosen,{'endpoint_pairs':len(groups),'ambiguous_sequences':ambiguous,'consensus_edges':len(edges),'reciprocal_edges':len(rec),'cycle_rejected':cycles}

def assemble(seedrows,chosen):
    inc={}
    for e in chosen:
        for ep in (e[0].le,e[0].re):
            if ep in inc:raise RuntimeError('endpoint reused')
            inc[ep]=e
    used=set();out=[]
    def orient(e,sid,end):
        p=e[0]
        return p if p.le==(sid,end) else p.flip()
    def chain(sid,rev):
        seq=rc(seedrows[sid][1]) if rev else seedrows[sid][1]; ids=[sid];used.add(sid)
        while True:
            rend='L' if rev else 'R';e=inc.get((sid,rend))
            if e is None:break
            p=orient(e,sid,rend); nid=p.re[0]; nrev=p.rs%2
            if nid in used:raise RuntimeError('cycle')
            ns=rc(seedrows[nid][1]) if nrev else seedrows[nid][1]
            if not p.cand.startswith(seq[-p.lo:]) or not ns.startswith(p.cand[-p.ro:]):raise RuntimeError('overlap mismatch')
            seq += p.cand[p.lo:] + ns[p.ro:]
            sid,rev=nid,nrev;ids.append(sid);used.add(sid)
        return ids,seq
    with_edges={x[0] for x in inc}
    for sid in sorted(with_edges):
        if sid in used:continue
        L=(sid,'L') in inc;R=(sid,'R') in inc
        if L and R:continue
        ids,s=chain(sid,L and not R);out.append((f'stage31_chain seeds={",".join(str(i+1) for i in ids)}',s))
    for sid in sorted(with_edges):
        if sid not in used:
            ids,s=chain(sid,False);out.append((f'stage31_chain seeds={",".join(str(i+1) for i in ids)}',s))
    for sid,(n,s) in enumerate(seedrows):
        if sid not in used:out.append((f'stage31_single seed={sid+1} source={n}',s))
    return out

def repl(cmd,flag,val):cmd[cmd.index(flag)+1]=str(val)
def build(bridgeasm,r1,r2,seed,root,k,threads,mode):
    vr1=root/f'v{k}_{mode}_1.fastq.gz';vr2=root/f'v{k}_{mode}_2.fastq.gz'
    nr,nb=carry.append_virtual_pairs(r1,r2,seed,vr1,vr2,200)
    out=root/mode
    cmd=s26.seeded_assemble_command(bridgeasm,vr1,vr2,out,k,16 if k==31 else 0,threads,major_path_cover=(mode=='k31_discovery'))
    if mode=='k31_discovery': repl(cmd,'--min-primary-support',2);repl(cmd,'--primary-dominance',.55);repl(cmd,'--path-cover-secondary-dominance',.10)
    elif mode=='k31_strict': pass
    else: repl(cmd,'--mercy-max-kmers',0);repl(cmd,'--min-primary-support',3);repl(cmd,'--primary-dominance',.70);repl(cmd,'--path-cover-secondary-dominance',.20)
    sec=carry.run(cmd,env=s26.clean_recovery_env());vr1.unlink(missing_ok=True);vr2.unlink(missing_ok=True)
    return out/'primary_contigs.fasta',{'k':k,'virtual_records':nr,'virtual_bases':nb,'seconds':sec}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bridgeasm',type=Path,required=True);ap.add_argument('--read1',type=Path,required=True);ap.add_argument('--read2',type=Path,required=True);ap.add_argument('--seed',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--threads',type=int,default=2);a=ap.parse_args()
    seedrows=list(fasta(a.seed));assert seedrows and all(not(set(s)-set(b'ACGT')) for _,s in seedrows)
    a.output.mkdir(parents=True,exist_ok=True); builds={}; paths={}
    for k,mode in [(31,'k31_discovery'),(31,'k31_strict'),(55,'k55_validator'),(77,'k77_validator'),(99,'k99_validator')]:
        p,info=build(a.bridgeasm,a.read1,a.read2,a.seed,a.output,k,a.threads,mode);paths[mode]=p;builds[mode]=info
    props=[];src={}
    for k,mode in [(31,'k31_discovery'),(55,'k55_validator'),(77,'k77_validator'),(99,'k99_validator')]:
        q,st=discover(seedrows,list(fasta(paths[mode])),k);props+=q;src[f'k{k}']=st
    variants={};scripts=Path(__file__).resolve().parent
    for name,(mh,mk,em) in VARIANTS.items():
        chosen,st=select(props,mh,mk,em,len(seedrows)); bridged=assemble(seedrows,chosen);vr=a.output/name;vr.mkdir(exist_ok=True)
        bfa=vr/'bridged_seed.fasta';write_fa(bfa,bridged)
        rows=['left_seed\tleft_end\tright_seed\tright_end\tsources\tscore\tinternal_bp']
        for p,ks,score in chosen:rows.append(f'{p.le[0]+1}\t{p.le[1]}\t{p.re[0]+1}\t{p.re[1]}\t{",".join(map(str,ks))}\t{score}\t{len(p.mid)}')
        (vr/'bridges.tsv').write_text('\n'.join(rows)+'\n')
        bst={**st,'selected_bridges':len(chosen),'seed_records':len(seedrows),'output_records':len(bridged),'seed_n50':n50(s for _,s in seedrows),'bridge_n50':n50(s for _,s in bridged),'max_chain_seeds':max((n.count(',')+1 for n,_ in bridged if n.startswith('stage31_chain')),default=1)}
        (vr/'bridge_stats.json').write_text(json.dumps(bst,indent=2,sort_keys=True)+'\n')
        final=vr/'primary_contigs.fasta';ext=vr/'extension_stats.json'
        carry.run([sys.executable,str(scripts/'seed_locked_extensions.py'),str(final),str(bfa),str(paths['k31_strict']),'--min-overlap','500','--overlap-margin','30','--seed-length','31','--min-extension','20','--max-seed-occurrences','64','--stats-json',str(ext)])
        variants[name]={'bridge':bst,'extension':json.loads(ext.read_text()),'final':str(final)}
    m={'pipeline':'stage31-multik-seed-bridge-v1','builds':builds,'source_stats':src,'total_proposals':len(props),'variants':variants}
    (a.output/'stage31_manifest.json').write_text(json.dumps(m,indent=2,sort_keys=True)+'\n');print(json.dumps(m,indent=2,sort_keys=True))
if __name__=='__main__':main()
