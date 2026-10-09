"""Native DP conversions and candidate reuse, with a numeric synthetic BPP."""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

binary=sys.argv[1]
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory);bpp=root/'synthetic.bpp';output=root/'out.bpseq';model=root/'model.txt'
    evidence={(1,10):.8,(2,9):.8,(5,14):.45,(6,13):.45}
    bpp.write_text('\n'.join(f'{i} A'+''.join(f' {j}:{p}' for (left,j),p in evidence.items() if left==i) for i in range(1,15))+'\n')
    def run(extra,threshold='0.4,0.4'):
        r=subprocess.run([binary,'--decoder','dd','--dd-dp','improved-beam','--dd-beam','100',
            '--dd-max-iter','50','--dd-patience','0','-x','-r','0','-t',threshold,'-B',str(output),
            '--loglevel','info',*extra,str(bpp)],capture_output=True,text=True,timeout=15)
        return r,r.stdout+r.stderr
    def pairs():
        return {(int(f[0]),int(f[2])) for line in output.read_text().splitlines() if len(f:=line.split())==3 and f[0].isdigit() and int(f[2])>int(f[0])}
    r,log=run([]);assert r.returncode==0,log;baseline=output.read_bytes()
    for kind in ('LINEAR','LOCAL'):
        model.write_text(f'IPKNOT_PK_DP_{kind}_V1\nintercept 0\nenergy_scale 1\ntemperature 37\n')
        for form in ('crossing','projected'):
            flags=['--pk-learned-model',str(model),'--pk-h-formulation',form]
            r,log=run(flags);assert r.returncode==0 and pairs()=={(1,10),(2,9)},log
            r,log=run([*flags,'--pk-learned-scale','0']);assert r.returncode==0 and output.read_bytes()==baseline,log
    r,log=run(['--pk-ensemble','--pk-ensemble-scale','1']);assert r.returncode==0 and output.read_bytes()==baseline,log # singleton pool
    r,baseline_log=run([],threshold='auto,auto');assert r.returncode==0,baseline_log
    dump=root/'pool.jsonl'
    r,log=run(['--pk-ensemble','--pk-ensemble-scale','.1','--pk-ensemble-output',str(dump)],threshold='auto,auto')
    assert r.returncode==0 and 'extra solves=0' in log,log
    pattern=r'DD: iterations=.*'
    assert re.findall(pattern,log)==re.findall(pattern,baseline_log),'ensemble changed DD problems'
    pool=json.loads(dump.read_text());assert pool['visits']>=len(pool['candidates'])
    assert abs(sum(c['probability'] for c in pool['candidates'])-1)<1e-12
    selected={(i+1,j+1) for i,j in pool['candidates'][pool['selected']]['pairs']}
    assert selected==pairs()
    marginal={}
    for c in pool['candidates']:
        for ij in c['pairs']:marginal[tuple(ij)]=marginal.get(tuple(ij),0)+c['probability']
    assert all(sum(p for ij,p in marginal.items() if i in ij)<=1+1e-12 for i in range(14))
    for flags in [['--pk-ensemble-scale','.1'],['--pk-ensemble','--pk-ensemble-temperature','0'],
        ['--pk-ensemble','--pk-ensemble-scale','nan'],['--pk-ensemble','--pk-level-penalty','1'],
        ['--pk-ensemble','--pk-ensemble-threshold','2']]:
        assert run(flags)[0].returncode!=0,flags
print('DP CLI signed/zero scores, singleton identity, original DD graph reuse, finite pool output verified')
