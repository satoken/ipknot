#!/usr/bin/env python3
"""Compare actual ILP objectives against an immutable pre-pruning binary."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import random
import re
import subprocess

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before',type=Path,default=Path('/tmp/ipknot-noe-optimized-build/ipknot'))
    parser.add_argument('--after',type=Path,default=Path('/tmp/ipknot-noe-pruned-build/ipknot'))
    parser.add_argument('--output',type=Path,default=Path('new_nmr/noe_pruning_results/random_regression'))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    rng=random.Random(20261004)
    cases=[]
    for seed in range(4):
        fasta=args.output/f'random_{seed}.fa'
        fasta.write_text(f'>random_{seed}\n'+''.join(rng.choice('ACGU') for _ in range(12+4*seed))+'\n')
        for level in (1,2):
            for bulge in ('none','all'):
                for shared in (False,True):
                    for isolated in (False,True):
                        name=f's{seed}_L{level}_{bulge}_shared{int(shared)}_isolated{int(isolated)}'
                        command=['-n','1','-e','lpc','-r','0','-t',','.join(['.2']*level),
                                 '--nmr-soft','--nmr-bulge-mode',bulge,'--coaxial-stacking',
                                 '--loglevel','debug','--base-pairs','AA=1,UU=1,GU=2']
                        for pattern in ('GC AU','GU GC','GC GC','GC AU AU'):
                            command += ['--stack-constraint',pattern]
                        if shared: command += ['--nmr-allow-shared-stack-pairs']
                        if isolated: command += ['--allow-isolated']
                        command.append(str(fasta))
                        cases.append((name,command))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    # The integrated executable defaults to DD; this compares ILP objectives.
    decoder_args = {}
    for binary in (args.before, args.after):
        help_result = subprocess.run([str(binary), '--help'], capture_output=True,
                                     text=True, check=True, env=env)
        decoder_args[binary] = ['--decoder', 'ilp'] if '--decoder' in help_result.stdout else []
    def compare(case):
        name,command=case
        values=[]
        structures=[]
        removed=0
        for label,binary in (('before',args.before),('after',args.after)):
            result=subprocess.run([str(binary),*decoder_args[binary],*command],capture_output=True,text=True,timeout=20,env=env)
            (args.output/f'{name}_{label}.log').write_text(result.stdout+result.stderr)
            if result.returncode: raise RuntimeError((name,label,result.stderr))
            objectives=re.findall(r'IP objective: ([^\n]+)',result.stdout)
            values.append([float(x) for x in objectives])
            structures.append(re.findall(r'^[.()\[\]{}<>]+$',result.stdout,re.M))
            if label=='after':
                removed=sum(int(n) for n in re.findall(r'and (\d+) coaxial proven infeasible',result.stdout))
        equal=len(values[0])==len(values[1]) and all(abs(a-b)<1e-7 for a,b in zip(*values))
        return dict(case=name,objectives=values,objective_equal=equal,
                    structure_equal=structures[0]==structures[1],coaxial_removed=removed)
    with ThreadPoolExecutor(max_workers=4) as pool: rows=list(pool.map(compare,cases))
    summary=dict(cases=len(rows),matching_objectives=sum(r['objective_equal'] for r in rows),
                 matching_structures=sum(r['structure_equal'] for r in rows),
                 cases_actually_pruned=sum(r['coaxial_removed']>0 for r in rows),results=rows)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='results'}))
    if not all(r['objective_equal'] for r in rows): raise SystemExit(1)

if __name__=='__main__': main()
