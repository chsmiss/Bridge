#!/usr/bin/env python3
"""Embed sparse BridgeBin DNA candidates with GenomeOcean.

GenomeOcean is trained directly on large-scale metagenomic assemblies and its official
embedding path uses a length-masked mean of the final hidden state.  This adapter mirrors
that pooling while adding forward/reverse-complement averaging and candidate gating so it
can act as an independent DNA expert beside DNABERT-S.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Set, Tuple

import torch
from transformers import AutoModel, PreTrainedTokenizerFast


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contigs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--pairs',type=Path,help='optional candidate table; embed only endpoints')
    p.add_argument('--model',default='DOEJGI/GenomeOcean-100M')
    p.add_argument('--device',default='cpu')
    p.add_argument('--window-bp',type=int,default=3000)
    p.add_argument('--max-windows-per-contig',type=int,default=2)
    p.add_argument('--max-tokens',type=int,default=1024)
    p.add_argument('--batch-size',type=int,default=2)
    p.add_argument('--no-reverse-complement',action='store_true')
    return p.parse_args(argv)


def fasta(path:Path)->Iterator[Tuple[str,str]]:
    name=None; chunks=[]
    with path.open() as h:
        for raw in h:
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name,''.join(chunks).upper()
                name=line[1:].split()[0]; chunks=[]
            else:
                if name is None: raise ValueError(f'{path}: sequence before header')
                chunks.append(line)
    if name is not None: yield name,''.join(chunks).upper()


def endpoints(path:Optional[Path])->Optional[Set[str]]:
    if path is None: return None
    out=set()
    with path.open(newline='') as h:
        r=csv.DictReader(h,delimiter='\t')
        for row in r:
            for key in ('left','right','source','target','contig_a','contig_b'):
                v=(row.get(key) or '').strip()
                if v: out.add(v)
    return out


def revcomp(seq:str)->str:
    return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


def positions(length:int,window:int,count:int)->List[int]:
    if length<=window or count<=1: return [0]
    count=min(count,max(1,math.ceil(length/window)))
    span=length-window
    return sorted({round(i*span/(count-1)) for i in range(count)})


def normalize(x:torch.Tensor)->torch.Tensor:
    return x/x.norm(p=2,dim=-1,keepdim=True).clamp_min(1e-12)


def pooled_hidden(model,tokenizer,sequences:List[str],device:torch.device,max_tokens:int)->torch.Tensor:
    enc=tokenizer(sequences,return_tensors='pt',padding=True,truncation=True,max_length=max_tokens)
    # GenomeOcean is Mistral-based. Some fast-tokenizer versions emit token_type_ids,
    # but MistralModel.forward does not accept them.
    enc.pop('token_type_ids',None)
    enc={k:v.to(device) for k,v in enc.items()}
    with torch.inference_mode():
        out=model(**enc,output_hidden_states=False,return_dict=True)
    hidden=out.last_hidden_state.float()
    mask=enc.get('attention_mask')
    if mask is None: mask=torch.ones(hidden.shape[:2],device=hidden.device,dtype=torch.long)
    mf=mask.to(hidden.dtype).unsqueeze(-1)
    pooled=(hidden*mf).sum(1)/mf.sum(1).clamp_min(1)
    return normalize(pooled).cpu()


def main(argv:Optional[Sequence[str]]=None)->int:
    a=parse_args(argv)
    if a.window_bp<256 or a.max_windows_per_contig<1 or a.max_tokens<32 or a.batch_size<1:
        raise SystemExit('invalid window/token/batch arguments')
    selected=endpoints(a.pairs)
    device=torch.device(a.device)
    tokenizer=PreTrainedTokenizerFast.from_pretrained(a.model)
    tokenizer.model_max_length=a.max_tokens
    tokenizer.padding_side='left'
    if tokenizer.pad_token_id is None:
        tokenizer.pad_token=tokenizer.eos_token or tokenizer.unk_token
    dtype=torch.float32 if device.type=='cpu' else torch.bfloat16
    model=AutoModel.from_pretrained(a.model,torch_dtype=dtype,attn_implementation='sdpa').to(device).eval()
    records=[]
    for contig,seq in fasta(a.contigs):
        if selected is not None and contig not in selected: continue
        for rank,start in enumerate(positions(len(seq),a.window_bp,a.max_windows_per_contig)):
            window=seq[start:start+a.window_bp]
            records.append((contig,rank,start,window,'fwd'))
            if not a.no_reverse_complement:
                records.append((contig,rank,start,revcomp(window),'rc'))
    if not records: raise SystemExit('no sequences selected')
    embeddings=[]
    for start in range(0,len(records),a.batch_size):
        batch=records[start:start+a.batch_size]
        vectors=pooled_hidden(model,tokenizer,[x[3] for x in batch],device,a.max_tokens)
        embeddings.extend(vectors.tolist())
    by={}
    for rec,vec in zip(records,embeddings):
        contig=rec[0]; by.setdefault(contig,[]).append(vec)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['contig','embedding','windows','model'])
        for contig in sorted(by):
            x=torch.tensor(by[contig],dtype=torch.float32).mean(0,keepdim=True); x=normalize(x)[0]
            w.writerow([contig,','.join(f'{float(v):.7g}' for v in x.tolist()),len(by[contig]),a.model])
    missing=0 if selected is None else len(selected-set(by))
    print(f'bridgebin-genomeocean: contigs={len(by)} windows={len(records)} selected={len(selected) if selected is not None else len(by)} missing_pair_endpoints={missing} model={a.model}')
    return 0

if __name__=='__main__': raise SystemExit(main())
