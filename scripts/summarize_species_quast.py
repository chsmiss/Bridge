#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
BACTERIA=['Bacillus_subtilis','Enterococcus_faecalis','Escherichia_coli','Lactobacillus_fermentum','Listeria_monocytogenes','Pseudomonas_aeruginosa','Salmonella_enterica','Staphylococcus_aureus']

def parse_report(path:Path):
    rows=list(csv.reader(path.open(),delimiter='\t'))
    if not rows:return {}
    assemblies=rows[0][1:]; out={a:{} for a in assemblies}
    for row in rows[1:]:
        if not row:continue
        metric=row[0]
        for a,v in zip(assemblies,row[1:]):out[a][metric]=v
    return out

def num(v):
    try:return float(v)
    except:return 0.0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('root',type=Path); ap.add_argument('--json',type=Path); a=ap.parse_args()
    species={}
    for sp in BACTERIA:
        report=a.root/sp/'report.tsv'
        if not report.exists(): raise SystemExit(f'missing {report}')
        species[sp]=parse_report(report)
    assemblies=sorted(next(iter(species.values())).keys())
    result={}; rows=['assembly\tbacteria_mean_gf\tbacteria_mean_na50\tlisteria_gf\tlisteria_na50\tpseudomonas_gf\tpseudomonas_na50']
    for asm in assemblies:
        gfs=[num(species[sp][asm].get('Genome fraction (%)',0)) for sp in BACTERIA]
        na=[num(species[sp][asm].get('NA50',0)) for sp in BACTERIA]
        rec={
          'bacteria_mean_gf':sum(gfs)/len(gfs),'bacteria_mean_na50':sum(na)/len(na),
          'listeria_gf':num(species['Listeria_monocytogenes'][asm].get('Genome fraction (%)',0)),
          'listeria_na50':num(species['Listeria_monocytogenes'][asm].get('NA50',0)),
          'pseudomonas_gf':num(species['Pseudomonas_aeruginosa'][asm].get('Genome fraction (%)',0)),
          'pseudomonas_na50':num(species['Pseudomonas_aeruginosa'][asm].get('NA50',0)),
        }
        result[asm]=rec
        rows.append(f"{asm}\t{rec['bacteria_mean_gf']:.4f}\t{rec['bacteria_mean_na50']:.1f}\t{rec['listeria_gf']:.4f}\t{rec['listeria_na50']:.1f}\t{rec['pseudomonas_gf']:.4f}\t{rec['pseudomonas_na50']:.1f}")
    print('\n'.join(rows))
    if a.json:a.json.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
