from pathlib import Path
import json, subprocess, sys, tempfile

binary=sys.argv[1]; linked=sys.argv[2]=='ON'
root=Path(__file__).resolve().parents[1]
def run(args, success=True, message=''):
    p=subprocess.run([binary,'--decoder','dd','-r','0','--loglevel','info',*map(str,args)],text=True,capture_output=True,timeout=40)
    text=p.stdout+p.stderr
    assert (p.returncode==0)==success,(args,text)
    assert message in text,(message,text)
    return text
with tempfile.TemporaryDirectory() as td:
    d=Path(td); bpp=d/'input.bpp';trace=d/'trace.jsonl'
    bpp.write_text('1 G 8:0.9\n2 G 7:0.9\n3 A\n4 C 11:0.9\n5 C 10:0.9\n6 A\n7 C\n8 C\n9 A\n10 G\n11 G\n')
    args=['-x','-t','.01,.01','--dd-dp','nussinov','--stack-constraint','GC GC','--nmr-soft','--dd-noe-solver','ilp','--dd-trace',trace,'--dd-trace-state',bpp]
    if not linked:
        run(args,False,'needs a linked ILP solver')
    else:
        for mode in ('linear','full'):
            run([*args[:-1],'--dd-constraints',mode,bpp])
            events=[json.loads(x)for x in trace.read_text().splitlines()]
            problem=next(x for x in events if x['event']=='problem')
            summary=next(x for x in events if x['event']=='summary')
            assert problem['noe_solver']=='ilp' and problem['noe_ilp_variables']>0
            assert summary['noe_ilp_calls']>0 and summary['noe_ilp_variables']==problem['noe_ilp_variables']
            assert summary['lower_bound']<=summary['upper_bound']+1e-7
            paircols={p[3] for p in problem['pairs']}
            assert paircols.isdisjoint(problem['noe_columns'])
            for e in events:
                if e['event']!='iteration':continue
                for bound,lower,upper,terms in problem['noe_rows']:
                    activity=sum(a*e['selected'][c] for c,a in terms)
                    if bound in (1,3,4):assert activity>=lower-1e-7
                    if bound in (2,3):assert activity<=upper+1e-7
                    if bound==4:assert activity<=lower+1e-7
        # Real coaxial witnesses, continuation pairs, shared faces, hard/soft.
        for kind in ('multibranch','last_child','children'):
            fixture=root/f'tests/data/nmr_coaxial_bulge_{kind}.fa'
            run(['--dd-noe-solver','ilp','--coaxial-stacking','--nmr-soft',
                 '--nmr-allow-shared-stack-pairs','--stack-constraint','UU CC',
                 '--stack-constraint','GC AU','--nmr-bulge-mode','all',
                 '-c',fixture.with_suffix('.constraint'),'-t','0,0',fixture])
            run(['--dd-noe-solver','ilp','--coaxial-stacking',
                 '--nmr-allow-shared-stack-pairs','--stack-constraint','UU CC',
                 '--stack-constraint','UU CC','--stack-constraint','GC AU',
                 '--nmr-bulge-mode','all','-c',fixture.with_suffix('.constraint'),
                 '-t','0,0',fixture],False)
        # No witnesses: a single forced violation is a valid exact NOE factor.
        run(['-x','-t','.01,.01','--dd-noe-solver','ilp','--nmr-soft',
             '--stack-constraint','UU UU',bpp])
    run(['-x','--dd-noe-solver','wrong',bpp],False,'must be relaxed or ilp')
