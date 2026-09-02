#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path

SPECIES={
 'Bacillus_subtilis':['bacillus subtilis','bacillus_subtilis','b_subtilis','bsubtilis'],
 'Enterococcus_faecalis':['enterococcus faecalis','enterococcus_faecalis','e_faecalis','efaecalis'],
 'Escherichia_coli':['escherichia coli','escherichia_coli','e_coli','ecoli'],
 'Lactobacillus_fermentum':['lactobacillus fermentum','lactobacillus_fermentum','l_fermentum','lfermentum'],
 'Listeria_monocytogenes':['listeria monocytogenes','listeria_monocytogenes','l_monocytogenes','lmonocytogenes'],
 'Pseudomonas_aeruginosa':['pseudomonas aeruginosa','pseudomonas_aeruginosa','p_aeruginosa','paeruginosa'],
 'Salmonella_enterica':['salmonella enterica','salmonella_enterica','s_enterica','senterica'],
 'Staphylococcus_aureus':['staphylococcus aureus','staphylococcus_aureus','s_aureus','saureus'],
 'Saccharomyces_cerevisiae':['saccharomyces cerevisiae','saccharomyces_cerevisiae'],
 'Cryptococcus_neoformans':['cryptococcus neoformans','cryptococcus_neoformans'],
}

def classify(path:Path):
    probe=str(path).lower().replace('-','_')
    try:
        with path.open(errors='ignore') as f:
            heads=' '.join(line.strip() for line in f if line.startswith('>'))[:20000].lower().replace('-','_')
    except Exception: heads=''
    probe+=' '+heads
    hits=[]
    for sp,terms in SPECIES.items():
        for t in terms:
            if t.replace(' ','_') in probe or t in probe:
                hits.append(sp); break
    return hits

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('root',type=Path); ap.add_argument('-o','--output',type=Path,required=True); a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    files=[p for p in a.root.rglob('*') if p.is_file() and p.suffix.lower() in {'.fa','.fasta','.fna','.fas'}]
    grouped={sp:[] for sp in SPECIES}; unknown=[]
    for p in files:
        hits=classify(p)
        if len(hits)==1: grouped[hits[0]].append(p)
        elif not hits: unknown.append(p)
    for sp,paths in grouped.items():
        if not paths: continue
        with (a.output/f'{sp}.fasta').open('w') as out:
            for p in sorted(paths): out.write(p.read_text())
        print(f'{sp}\t{len(paths)}\t'+','.join(str(p) for p in paths))
    if unknown:
        print('UNCLASSIFIED',*(str(p) for p in unknown),sep='\n')
    missing=[sp for sp in list(SPECIES)[:8] if not grouped[sp]]
    if missing: raise SystemExit('missing bacterial references: '+','.join(missing))
if __name__=='__main__':main()
