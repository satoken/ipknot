"""Graph parity, numeric scoring, zero identities and native rank selection."""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory() as directory:
    root=Path(directory);bpp=root/'numeric.bpp';out=root/'out.bpseq';model=root/'model.txt'
    evidence={(1,14):.6,(2,13):.6,(3,12):.6,(4,11):.6,
              (7,20):.2,(8,19):.2,(9,18):.2,(10,17):.2}
    bpp.write_text('\n'.join(f'{i} A'+''.join(f' {j}:{p}' for (a,j),p in evidence.items() if a==i) for i in range(1,21))+'\n')
    def run(flags,threshold='auto,auto'):
        r=subprocess.run([sys.argv[1],'-x','-r','0','-t',threshold,'-B',str(out),
            '--decoder','dd','--dd-dp','improved-beam','--dd-beam','100','--dd-max-iter','50',
            '--dd-patience','0','--dd-crossing-beam','100','--dd-witnesses','16',
            '--loglevel','info',*flags,str(bpp)],capture_output=True,text=True,timeout=15)
        return r,r.stdout+r.stderr
    r,log=run([]);assert r.returncode==0,log;baseline=out.read_bytes()
    model.write_text('IPKNOT_PK_EXCLUSION_V1\nintercept 0\nenergy_scale 0\ntemperature 37\n')
    common=['--pk-learned-model',str(model),'--pk-h-max-stem','200','--pk-h-max-loop','200']
    for core in (0,2,3):
        for form in ('crossing','projected'):
            flags=[*common,'--pk-h-formulation',form,'--pk-core-width',str(core)]
            r,log=run([*flags,'--pk-learned-scale','0']);assert r.returncode==0 and out.read_bytes()==baseline,log
            r,log=run(flags);assert r.returncode==0,log
            scores=[float(x) for x in re.findall(r'PK_score=([^,]+),',log)]
            assert scores and max(scores)>0,log
            if form=='crossing':
                r,log=run([*flags,'--pk-best-partner']);assert r.returncode==0,log
                r,log=run([*flags,'--pk-best-partner','--pk-learned-scale','0']);assert r.returncode==0 and out.read_bytes()==baseline,log
    pool=root/'pool.jsonl';rank=root/'rank.txt'
    rank.write_text('IPKNOT_PK_RANK_V1\n'+''.join(f'f{k} {(-10 if k==3 else 0)}\n' for k in range(8)))
    flags=[*common,'--pk-h-formulation','projected','--pk-core-width','2']
    r,base_log=run(flags);assert r.returncode==0,base_log;unranked=out.read_bytes()
    r,log=run([*flags,'--pk-rank-model',str(rank),'--pk-rank-scale','0']);assert r.returncode==0 and out.read_bytes()==unranked,log
    r,log=run([*flags,'--pk-rank-model',str(rank),'--pk-rank-scale','1','--pk-rank-output',str(pool)])
    assert r.returncode==0,log
    assert re.findall(r'DD: iterations=.*',log)==re.findall(r'DD: iterations=.*',base_log)
    record=json.loads(pool.read_text());candidates=record['candidates']
    scores=[c['base']-10*c['features'][3] for c in candidates]
    assert record['selected']==max(range(len(scores)),key=lambda k:scores[k])
    physical={(int(f[0])-1,int(f[2])-1) for line in out.read_text().splitlines()
              if len(f:=line.split())==3 and f[0].isdigit() and int(f[2])>int(f[0])}
    assert physical=={tuple(p) for p in candidates[record['selected']]['pairs']}
    for bad in [['--pk-core-width','1'],['--pk-best-partner','--pk-h-formulation','projected'],
                ['--pk-rank-scale','1'],['--pk-rank-model',str(rank),'--pk-rank-scale','-1']]:
        assert run(bad)[0].returncode!=0,bad
    assert run(['--pk-rank-model',str(rank)],threshold='.2,.2')[0].returncode!=0
print('alternative CLI: short cores, both forms, best-partner, zero parity, rank/DD graph reuse passed')
