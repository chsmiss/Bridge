#!/usr/bin/env python3
import random
import stage31_multik_seed_bridge as s

def dna(seed,n):
    r=random.Random(seed)
    return bytes(r.choice(b'ACGT') for _ in range(n))

def main():
    a,b,c=dna(1,900),dna(2,900),dna(3,900)
    x,y=dna(4,150),dna(5,170)
    seeds=[('a',a),('b',b),('c',c)]
    props=[]
    for k in (31,55,77):
        props += [
            s.P(k,0,2,160,160,x,a[-160:]+x+b[:160],f'ab{k}'),
            s.P(k,2,4,160,160,y,b[-160:]+y+c[:160],f'bc{k}'),
        ]
    chosen,_=s.select(props,2,77,2,3)
    assert len(chosen)==2
    out=s.assemble(seeds,chosen)
    assert len(out)==1
    assert all(z in out[0][1] for z in (a,b,c))
    assert b'N' not in out[0][1]
    cand=a[-170:]+x+b[:190]
    q,stats=s.discover([('a',a),('b',b)],[('ab',cand)],77)
    assert len(q)==1,(q,stats)
    print('stage31 tests: passed')

if __name__=='__main__':
    main()
