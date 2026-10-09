"""Length-doubling audit for bounded fixed/NMR decoding on sparse BPP input.

Usage: python3 tests/dd_constraint_scaling_test.py /path/to/ipknot \
    --lengths 512,1024,2048,4096,8192 --repeat 3 --report scaling.json
Counters, rather than wall-time ratios, are the deterministic regression check.
"""
from pathlib import Path
import argparse
import json
import math
import platform
import re
import statistics
import subprocess
import tempfile
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('binary')
parser.add_argument('--lengths', default='128,256,512,1024')
parser.add_argument('--repeat', type=int, default=1)
parser.add_argument('--report', type=Path)
args = parser.parse_args()
lengths = [int(value) for value in args.lengths.split(',')]
assert lengths and min(lengths) >= 16 and args.repeat > 0
root = Path(__file__).resolve().parents[1]
budgets = dict(iterations=3, states=32, passes=4, dp=8, crossing=8,
               contacts=8, pairs=8, witnesses=8, pattern=8, levels=2)
flags = ['--decoder','dd','-x','-r','0','-t','.01,.01', '--loglevel','info',
         '--dd-max-iter','3','--dd-constraint-states','32','--dd-constraint-passes','4',
         '--dd-nmr-pair-beam','8','--dd-nmr-witnesses','8','--dd-nmr-pattern-beam','8',
         '--dd-beam','8','--dd-crossing-beam','8','--dd-witnesses','8']
records = []
with tempfile.TemporaryDirectory() as directory:
    tmp = Path(directory)
    for family in ('nested_fixed', 'noncanonical_soft', 'coaxial_fixed'):
        previous = None
        for requested_length in lengths:
            fixed_pairs = {}
            fixed_file = tmp/'fixed.constraint'
            if family == 'nested_fixed':
                half = (requested_length-4)//2
                seq = 'G'*half+'AAAA'+'C'*half
                fixed_pairs = {i:len(seq)-1-i for i in range(half)}
                observations = ['-c',str(fixed_file),'--base-pairs',f'GC={half}',
                                '--stack-constraint','GC GC']
            elif family == 'noncanonical_soft':
                seq = 'U'*requested_length
                observations = ['--base-pairs','UU=1','--without-canonical-neighbor',
                                '--stack-constraint','UU UU','--nmr-soft',
                                '--nmr-count-penalty','.0001','--nmr-stack-penalty','.0001']
            else:
                fixture = root/'tests/data/nmr_coaxial_bulge_multibranch'
                unit = fixture.with_suffix('.fa').read_text().splitlines()[1]
                copies = math.ceil(requested_length/len(unit))
                seq = unit*copies
                local = [(int(i)-1,int(j)-1) for i,j in
                         (line.split() for line in fixture.with_suffix('.constraint').read_text().splitlines())
                         if j != 'x' and int(i) < int(j)]
                fixed_pairs = {i+k*len(unit):j+k*len(unit)
                               for k in range(copies) for i,j in local}
                observations = ['-c',str(fixed_file),'--coaxial-stacking',
                                '--stack-constraint','UU CC','--stack-constraint','GC AU']
            n = len(seq)
            partners = {**fixed_pairs, **{j:i for i,j in fixed_pairs.items()}}
            fixed_file.write_text(''.join(f'{i+1} {partners[i]+1 if i in partners else "x"}\n'
                                          for i in range(n)))
            bpp = tmp/'input.bpp'
            bpp.write_text(''.join(f'{i+1} {base}' +
                                   (f' {fixed_pairs[i]+1}:0.9' if i in fixed_pairs else '')+'\n'
                                   for i,base in enumerate(seq)))
            trace = tmp/'trace.jsonl'
            output = tmp/'output.bpseq'
            durations = []
            for repetition in range(args.repeat):
                trace.unlink(missing_ok=True)
                start = time.perf_counter()
                result = subprocess.run([args.binary,*flags,'--dd-trace',str(trace),
                                         '-B',str(output),*observations,str(bpp)],
                                        capture_output=True,text=True,timeout=90)
                durations.append(time.perf_counter()-start)
                text = result.stdout+result.stderr
                assert result.returncode == 0, (family,n,text)
                events = [json.loads(line) for line in trace.read_text().splitlines()]
                problem = next(e for e in events if e['event']=='problem')
                summary = next(e for e in events if e['event']=='summary')
                logs = re.search(r'DD: iterations=.*',text).group(0)
                counters = {key:float(value) for key,value in
                            re.findall(r'(\w+)=(-?[\d.]+(?:e[+-]?\d+)?)',logs)}
                assert counters['linear'] == 1
                assert summary['upper_bound'] >= summary['lower_bound']-1e-7
                assert summary['iterations'] <= budgets['iterations']
                assert summary['repair_states'] <= budgets['states']
                assert summary['propagation_work'] <= (2*budgets['passes']*(budgets['states']+1)*summary['nonzeros'])
                planar_unit = 3*n*budgets['levels']+counters['pairs']
                assert summary['structural_work'] <= ((budgets['passes']*(budgets['states']+1)+
                                                       budgets['states']+budgets['iterations']+3)*planar_unit)
                prediction = {int(i)-1:int(j)-1 for line in output.read_text().splitlines()
                              if line and line[0].isdigit() for i,base,j in [line.split()]}
                assert prediction == {i:partners.get(i,-1) for i in range(n)}, (family,n)
                if family == 'noncanonical_soft':
                    # Any pair costs >=.01, more than both soft penalties.
                    assert abs(summary['lower_bound']+.0002)<1e-8
                    assert counters['nmr_pair_drops'] > 0
                else:
                    assert 'All stack constraints are satisfied' in text
            # These ceilings distinguish retained input graphs from a dense
            # pair or pair-pair matrix on the same sequences.
            assert problem['variables'] <= 24*n
            assert problem['row_count'] <= 20*n
            assert summary['nonzeros'] <= 320*n
            assert counters['nmr_seed_visits'] <= 32*n
            record = dict(family=family,length=n,input_pairs=len(fixed_pairs),
                          seconds_median=round(statistics.median(durations),6),
                          variables=problem['variables'],rows=problem['row_count'],
                          **{key:summary[key] for key in ('iterations','repair_states','nonzeros',
                              'propagation_work','structural_work','lower_bound','upper_bound')},
                          nmr_seed_visits=int(counters['nmr_seed_visits']))
            if previous:
                growth = n/previous['length']
                for field in ('variables','rows','nonzeros','nmr_seed_visits'):
                    assert record[field] <= 1.16*growth*max(previous[field],1)+32, (family,field,previous,record)
            records.append(record)
            previous = record
report = dict(machine=platform.machine(),budgets=budgets,repeat=args.repeat,records=records)
if args.report:
    args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
