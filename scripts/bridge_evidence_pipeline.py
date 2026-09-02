#!/usr/bin/env python3
"""Cumulative six-step BridgeAsm evidence recovery experiment.

step0 baseline                 existing recovery-v3
step1 iterative                low-k trusted paths -> virtual pairs -> next k graph
step2 residual                 shared-backbone residual path cover from GFA evidence
step3 pairs                    empirical insert-size paired-end link graph
step4 second_pass              remap reads to preliminary scaffolds and resolve again
step5 gapfill                  local read-supported DBG fill of N gaps
step6 strain_projection        project flanked bulges, feed them back as strain-path evidence
"""
from __future__ import annotations
import argparse, json, os, resource, shutil, subprocess, sys, time
from pathlib import Path

def run(cmd,*,env=None,stdout=None):
    print('+',' '.join(map(str,cmd)),flush=True); st=time.monotonic()
    subprocess.run(list(map(str,cmd)),check=True,env=env,stdout=stdout)
    return time.monotonic()-st

def bridge_cmd(exe,r1,r2,out,k,mercy,threads):
    return [exe,'assemble','-1',r1,'-2',r2,'-o',out,'-k',k,'--min-count',2,'--mercy-max-kmers',mercy,
            '--mercy-min-support',1,'--mercy-min-quality',25,'--min-read-support',2,'--min-pair-support',2,
            '--min-primary-support',5,'--primary-dominance',0.75,'--threaded-path-cover','--major-path-cover',
            '--path-cover-secondary-dominance',0.25,'--min-contig-length',200,'--threads',threads]

def concat_gzip(a:Path,b:Path,out:Path):
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('wb') as w:
        for p in (a,b):
            with p.open('rb') as r: shutil.copyfileobj(r,w,1<<20)

def postprocess(scripts:Path,inputs:list[Path],outdir:Path,min_overlap=31,margin=10):
    outdir.mkdir(parents=True,exist_ok=True)
    merged=outdir/'union.fasta'; filtered=outdir/'noncontained.fasta'; final=outdir/'primary_contigs.fasta'
    run([sys.executable,scripts/'merge_fasta_unique.py',merged,*inputs,'--min-length',200])
    run([sys.executable,scripts/'filter_contained_fasta.py',merged,filtered,'--min-length',200,'--seed-k',21,'--window',12,
         '--candidate-minimizers',16,'--removed-tsv',outdir/'contained_removed.tsv','--stats-json',outdir/'containment_stats.json'])
    run([sys.executable,scripts/'stitch_exact_overlaps.py',final,filtered,'--min-overlap',min_overlap,'--overlap-margin',margin,
         '--seed-length',31,'--max-seed-occurrences',64,'--min-length',200])
    return final

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--bridgeasm',type=Path,required=True); ap.add_argument('--read1',type=Path,required=True); ap.add_argument('--read2',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--threads',type=int,default=2)
    a=ap.parse_args(); out=a.output; out.mkdir(parents=True,exist_ok=True); scripts=Path(__file__).resolve().parent; times={}; started=time.monotonic()

    step0=out/'step0_baseline'
    times['step0_baseline']=run([sys.executable,scripts/'bridge_recovery_pipeline.py','--bridgeasm',a.bridgeasm,'--read1',a.read1,'--read2',a.read2,'--output',step0,'--threads',a.threads,'--singleton-fraction',0.50,'--singleton-quality',35,'--mate-terminal-mercy',96,'--stitch-min-overlap',31,'--stitch-overlap-margin',10])

    iterative=out/'iterative'; iterative.mkdir(exist_ok=True)
    stages=[('k21_recall',21,24),('k31_resolve',31,16),('k41_resolve',41,12),('k55_resolve',55,8)]
    candidates=[]; cur_r1=a.read1; cur_r2=a.read2
    for idx,(name,k,mercy) in enumerate(stages):
        sd=iterative/name; env=os.environ.copy()
        if k==21:
            env['BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION']='0.50'; env['BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY']='35'; env['BRIDGEASM_MATE_TERMINAL_MERCY_KMERS']='96'
        else:
            env.pop('BRIDGEASM_SINGLETON_ISLAND_MIN_FRACTION',None); env.pop('BRIDGEASM_SINGLETON_ISLAND_MIN_QUALITY',None); env.pop('BRIDGEASM_MATE_TERMINAL_MERCY_KMERS',None)
        times[f'step1_{name}']=run(bridge_cmd(a.bridgeasm,cur_r1,cur_r2,sd,k,mercy,a.threads),env=env)
        candidates += [sd/'primary_contigs.fasta',sd/'haplotigs.fasta']
        if idx+1<len(stages):
            vr1=iterative/f'{name}.virtual_R1.fastq.gz'; vr2=iterative/f'{name}.virtual_R2.fastq.gz'
            times[f'step1_{name}_virtualize']=run([sys.executable,scripts/'make_virtual_pairs.py',sd/'primary_contigs.fasta',sd/'haplotigs.fasta','--read1',vr1,'--read2',vr2,'--read-length',101,'--insert-size',250,'--stride',180,'--min-length',500])
            ar1=iterative/f'{name}.aug_R1.fastq.gz'; ar2=iterative/f'{name}.aug_R2.fastq.gz'; concat_gzip(a.read1,vr1,ar1); concat_gzip(a.read2,vr2,ar2); cur_r1,cur_r2=ar1,ar2
    s1dir=out/'step1_iterative'; s1=postprocess(scripts,candidates,s1dir); shutil.copy2(s1,out/'step1_iterative.fasta')

    residual=out/'step2_residual_paths.fasta'; residual_meta=out/'step2_residual_paths.tsv'
    times['step2_residual_extract']=run([sys.executable,scripts/'residual_path_cover.py',iterative/'k31_resolve'/'assembly.gfa',iterative/'k41_resolve'/'assembly.gfa',iterative/'k55_resolve'/'assembly.gfa','-o',residual,'--metadata',residual_meta,'--secondary-dominance',0.18,'--extension-dominance',0.72,'--min-support',2,'--min-length',300,'--max-copy',6])
    s2=postprocess(scripts,[s1,residual],out/'step2_residual'); shutil.copy2(s2,out/'step2_residual.fasta')

    s3=out/'step3_pairs.fasta'
    times['step3_pair_graph']=run([sys.executable,scripts/'pair_gap_refine.py',s2,'-1',a.read1,'-2',a.read2,'-o',s3,'--links',out/'step3_pair_links.tsv','--threads',a.threads,'--min-mapq',20,'--min-support',3,'--dominance',0.75,'--end-window',500,'--min-overlap',31,'--max-gap',1000])

    s4=out/'step4_second_pass.fasta'
    times['step4_second_pass']=run([sys.executable,scripts/'pair_gap_refine.py',s3,'-1',a.read1,'-2',a.read2,'-o',s4,'--links',out/'step4_pair_links.tsv','--threads',a.threads,'--min-mapq',20,'--min-support',3,'--dominance',0.75,'--end-window',600,'--min-overlap',31,'--max-gap',1000])

    s5=out/'step5_gapfill.fasta'
    times['step5_gapfill']=run([sys.executable,scripts/'fill_scaffold_gaps.py',s4,'-1',a.read1,'-2',a.read2,'-o',s5,'--report',out/'step5_gapfill.tsv','--anchor-k',31,'--local-k',21,'--flank',180,'--dominance',0.65])

    projected=out/'step6_projected_strain_paths.fasta'; projection_map=out/'step6_projection_map.tsv'
    haplotigs=[iterative/name/'haplotigs.fasta' for name,_,_ in stages]
    times['step6_project']=run([sys.executable,scripts/'strain_projection.py',s5,*haplotigs,'-o',projected,'--map',projection_map,'--k',31,'--projection-k',21,'--min-novel-fraction',0.03,'--min-projection-fraction',0.10])
    pvr1=out/'step6.virtual_R1.fastq.gz'; pvr2=out/'step6.virtual_R2.fastq.gz'
    times['step6_virtualize']=run([sys.executable,scripts/'make_virtual_pairs.py',projected,'--read1',pvr1,'--read2',pvr2,'--read-length',101,'--insert-size',250,'--stride',120,'--min-length',300])
    par1=out/'step6.aug_R1.fastq.gz'; par2=out/'step6.aug_R2.fastq.gz'; concat_gzip(a.read1,pvr1,par1); concat_gzip(a.read2,pvr2,par2)
    projasm=out/'step6_projection_k31'
    times['step6_projection_k31']=run(bridge_cmd(a.bridgeasm,par1,par2,projasm,31,16,a.threads))
    novel_proj=out/'step6_projection_novel.fasta'
    times['step6_select_novel']=run([sys.executable,scripts/'strain_projection.py',s5,projasm/'primary_contigs.fasta',projasm/'haplotigs.fasta','-o',novel_proj,'--map',out/'step6_projection_assembly_map.tsv','--k',31,'--projection-k',21,'--min-novel-fraction',0.02,'--min-projection-fraction',0.05])
    s6=postprocess(scripts,[s5,projected,novel_proj],out/'step6_strain_projection'); shutil.copy2(s6,out/'step6_strain_projection.fasta')

    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest={'pipeline':'bridge-evidence-six-step-v1','wall_seconds':time.monotonic()-started,'peak_child_rss_mib':usage.ru_maxrss/1024.0,'timings_seconds':times,
              'outputs':{'step0_baseline':str(step0/'primary_contigs.fasta'),'step1_iterative':str(out/'step1_iterative.fasta'),'step2_residual':str(out/'step2_residual.fasta'),'step3_pairs':str(s3),'step4_second_pass':str(s4),'step5_gapfill':str(s5),'step6_strain_projection':str(out/'step6_strain_projection.fasta')}}
    (out/'pipeline_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n'); print(json.dumps(manifest,indent=2,sort_keys=True))
if __name__=='__main__':main()
