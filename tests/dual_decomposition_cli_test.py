"""Exercise solver-free decoding, configuration errors and native refinement."""
from pathlib import Path
import json
import math
import subprocess
import sys
import tempfile

binary = sys.argv[1]
fixture = Path(__file__).resolve().parents[1] / 'examples/drz_Ppac_1_1.fa'
def run(*args):
    return subprocess.run([binary, *map(str,args)], capture_output=True, text=True, timeout=20)
def require_success(result):
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'DD: iterations=' in result.stdout + result.stderr

def check_recovered(problem, state):
    candidates = problem['pairs']
    selected = set(state['recovered'])
    assert len(selected) == len(state['recovered'])
    used = set()
    coordinates = set()
    for identifier in selected:
        left, right, level, _, allowed = candidates[identifier]
        assert allowed and left not in used and right not in used
        used.update((left, right))
        coordinates.add((left, right, level))
    for a in selected:
        left, right, level, _, _ = candidates[a]
        for b in selected:
            other_left, other_right, other_level, _, _ = candidates[b]
            assert level != other_level or not (
                left < other_left < right < other_right or
                other_left < left < other_right < right)
        if problem['no_lonely']:
            assert ((left + 1, right - 1, level) in coordinates or
                    (left - 1, right + 1, level) in coordinates)
    objective = sum(candidates[x][3] for x in selected)
    required = {identifier: set() for identifier in selected}
    for upper, contacts in problem['rows']:
        if upper not in selected:
            continue
        chosen = [(lower, score) for lower, score in contacts if lower in selected]
        assert chosen, 'Recovered upper pair has no required lower witness'
        required[upper].add(candidates[chosen[0][0]][2])
        objective += sum(score for _, score in chosen)
    for identifier in selected:
        assert required[identifier] == set(range(candidates[identifier][2]))
    assert math.isclose(objective, state['primal'], rel_tol=1e-10, abs_tol=1e-9)
    assert math.isfinite(state['upper_bound']) and math.isfinite(state['lower_bound'])
    assert state['upper_bound'] >= state['lower_bound'] - 1e-8

with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bpp = root / 'small.bpp'
    bpp.write_text('1 A 8:0.9\n2 A 7:0.9\n3 A\n4 A 11:0.9\n5 A 10:0.9\n6 A\n7 A\n8 A\n9 A\n10 A\n11 A\n')
    common = ['--decoder','dd','--loglevel','info','-x','-r','0','-t','.2,.1']
    direct=run(*common,bpp); require_success(direct)
    assert '((.[[.)).]]' in direct.stdout, direct.stdout
    # Omitting decoder and DD controls must use the requested DD defaults,
    # even when the executable is linked to an ILP solver.
    default_trace=root/'default.jsonl'; explicit_trace=root/'explicit.jsonl'
    default_output=root/'default.bpseq'; explicit_output=root/'explicit.bpseq'
    implicit=['--loglevel','info','-x','-r','0','-t','.2,.1']
    default_result=run(*implicit,'--dd-trace',default_trace,'-B',default_output,bpp)
    require_success(default_result)
    assert 'Decoder: dd' in default_result.stdout+default_result.stderr
    explicit_result=run(*common,'--dd-beam','100','--dd-max-iter','50',
                        '--dd-patience','0','--dd-crossing-beam','100','--dd-witnesses','16',
                        '--dd-trace',explicit_trace,'-B',explicit_output,bpp)
    require_success(explicit_result)
    def problem(path):
        return next(event for event in map(json.loads,path.read_text().splitlines())
                    if event['event']=='problem')
    assert problem(default_trace)==problem(explicit_trace), 'Implicit DD settings differ from explicit defaults'
    def numeric_lines(path):
        return [line for line in path.read_text().splitlines() if line and line[0].isdigit()]
    assert numeric_lines(default_output)==numeric_lines(explicit_output)
    # Exact Nussinov is a DP selection, independent of the configured beam,
    # and must also reach the unconstrained decoder and primal proposals.
    exact_outputs=[]
    for flags in (['--dd-dp','nussinov','--dd-beam','1'], ['--dd-beam','0']):
        exact_trace=root/'exact.jsonl'; exact_output=root/'exact.bpseq'
        result=run(*common,*flags,'--dd-recovery-every','1',
                   '--dd-trace',exact_trace,'-B',exact_output,bpp)
        require_success(result)
        assert problem(exact_trace)['beam']==0
        exact_outputs.append(numeric_lines(exact_output))
    assert exact_outputs[0]==exact_outputs[1], 'DP selector used the beam in exact Nussinov mode'
    trace=root/'trace.jsonl';trace.write_text('stale data\n')
    traced=run(*common,'--dd-trace',trace,'--dd-trace-state','--dd-unpruned-bound',bpp)
    require_success(traced)
    events=[json.loads(line) for line in trace.read_text().splitlines()]
    problem=next(x for x in events if x['event']=='problem')
    states=[x for x in events if x['event']=='iteration']
    summary=next(x for x in events if x['event']=='summary')
    assert len(states)==summary['iterations'] and states[-1]['stop']==summary['stop']
    assert problem['pairs'] and problem['rows'] and all('weights' in x for x in states)
    assert states[-1]['upper_bound']>=states[-1]['lower_bound']-1e-9
    assert summary['stop']=='bound_gap', summary
    for flags in (['--dd-global-bound'], ['--dd-bound-block','8','--dd-bound-every','1'],
                  ['--dd-recovery-every','10'],['--dd-recovery-every','10','--dd-recovery-target','best'], ['--dd-global-bound','--dd-bound-block','16','--dd-recovery-every','5']):
        traced=run(*common,'--dd-trace',trace,'--dd-trace-state',*flags,bpp)
        require_success(traced)
        events=[json.loads(line) for line in trace.read_text().splitlines()]
        summary=next(x for x in events if x['event']=='summary')
        assert summary['upper_bound']>=summary['lower_bound']-1e-9
        assert summary['lower_bound']>=1.5-1e-6
    require_success(run(*common,'--dd-schedule','diminishing','--dd-patience','0',bpp))
    # Several threshold subproblems must have distinct trace IDs.
    require_success(run('--decoder','dd','--dd-trace',trace,'--loglevel','info','-r','0','-x',bpp))
    ids=[x['solve_id'] for x in map(json.loads,trace.read_text().splitlines()) if x['event']=='summary']
    assert len(ids)>1 and len(ids)==len(set(ids))
    # Integer certificates may lie below the DD LP: check original feasible
    # structures and finite certificates, not certificate>=current dual value.
    for flags in (
            ['--dd-joint-bound','12','--dd-joint-states','0'],
            ['--dd-joint-bound','8','--dd-joint-clusters'],
            ['--dd-joint-bound','12','--dd-joint-clusters','--dd-exchange','12'],
            ['--dd-joint-bound','8','--dd-joint-shift','--dd-joint-matching'],
            ['--dd-joint-bound','8','--dd-joint-states','1'],
            ['--dd-exchange','12','--dd-exchange-passes','2'],
            ['--dd-exchange','8','--dd-exchange-passes','4','--dd-exchange-states','1'],
            ['--dd-joint-bound','8','--dd-joint-shift','--dd-exchange','8',
             '--dd-bound-block','8','--dd-bound-shift','--dd-bound-strict-stack',
             '--dd-exchange-every','2','--dd-recovery-every','2']):
        result = run(*common,
                     '--dd-trace',trace,'--dd-trace-state',*flags,bpp)
        require_success(result)
        events = [json.loads(line) for line in trace.read_text().splitlines()]
        problem = next(x for x in events if x['event'] == 'problem')
        summary = next(x for x in events if x['event'] == 'summary')
        states = [x for x in events if x['event'] == 'iteration']
        assert len(states) == summary['iterations'] and states[-1]['stop'] == summary['stop']
        assert math.isfinite(summary['upper_bound']) and math.isfinite(summary['lower_bound'])
        assert summary['upper_bound'] >= summary['lower_bound'] - 1e-8
        for state in states:
            check_recovered(problem, state)
        if '--dd-joint-bound' in flags:
            assert problem['joint_windows'] > 0 and math.isfinite(problem['static_certificate'])
            if '--dd-joint-states' in flags and flags[flags.index('--dd-joint-states') + 1] == '1':
                assert problem['joint_fallbacks'] > 0
        if '--dd-exchange' in flags:
            assert problem['exchange_width'] > 0
            assert all(summary[name] >= 0 for name in (
                'exchange_windows','exchange_improvements','exchange_budget_windows','exchange_states'))
    # Cached static proposals and shared decoder buffers preserve the original
    # recovery trajectory under both ordinary and improved beam decoding.
    for dp in ('beam', 'improved-beam'):
        recovered_outputs = []
        recovered_states = []
        for optimization in ('true', 'false'):
            output = root / f'recovery-{dp}-{optimization}.bpseq'
            result = run(*common,'--dd-dp',dp,'--dd-recovery-every','1',
                         '--dd-max-iter','40','--dd-patience','0',
                         f'--dd-recovery-cache={optimization}',f'--dd-recovery-share={optimization}',
                         '--dd-trace',trace,'--dd-trace-state','-B',output,bpp)
            require_success(result)
            recovered_outputs.append([line for line in output.read_text().splitlines()
                                      if line and line[0].isdigit()])
            events = [json.loads(line) for line in trace.read_text().splitlines()]
            recovered_states.append([(state['oracle_value'], state['baseline_lower_bound'],
                                      state['lower_bound'], state['eta'], state['selected'])
                                     for state in events if state['event'] == 'iteration'])
        assert recovered_outputs[0] == recovered_outputs[1], 'Recovery caching/sharing changed structure'
        assert recovered_states[0] == recovered_states[1], 'Recovery caching/sharing changed DD trajectory'
    for extra in ([], ['--dd-recovery-every', '1', '--dd-exchange', '6',
                       '--dd-unpruned-bound', '--dd-global-bound']):
        traced = run(*common, '--dd-dp', 'improved-beam', '--dd-beam', '100',
                     '--dd-max-iter', '50', '--dd-trace', trace,
                     '--dd-trace-state', *extra, bpp)
        require_success(traced)
        events = [json.loads(line) for line in trace.read_text().splitlines()]
        graph = next(x for x in events if x['event'] == 'problem')
        assert graph['improved_beam'] == 1 and graph['beam'] == 100
        assert all(score == 0 for _, contacts in graph['rows'] for _, score in contacts)
        summary = next(x for x in events if x['event'] == 'summary')
        assert summary['iterations'] <= 50
        for state in events:
            if state['event'] == 'iteration':
                check_recovered(graph, state)
    require_success(run('--decoder', 'dd', '--loglevel', 'info',
                        '--dd-dp', 'improved-beam', '-r', '1', '-t', 'auto,auto', fixture))
    # Empty candidates and one-level prediction remain valid.
    require_success(run(*common[:-2],'-t','1',bpp))
    for extra in (['--dd-max-iter','0'],['--dd-beam','-1'],['--dd-crossing-beam','-1'],
                  ['--dd-witnesses','-1'],['--dd-bound-block','65'],['--dd-bound-every','0'],
                  ['--dd-joint-bound','8','--dd-joint-clusters','--dd-joint-shift'],['--dd-joint-bound','-1'],['--dd-joint-bound','13'],['--dd-joint-states','-1'],
                  ['--dd-exchange','-1'],['--dd-exchange','13'],['--dd-exchange-states','-1'],
                  ['--dd-exchange-passes','0'],['--dd-exchange-passes','5'],['--dd-exchange-every','-1'],
                  ['--dd-recovery-every','-1'],['--dd-recovery-mode','wrong'],['--dd-recovery-target','wrong'],['--dd-schedule','wrong'],['--dd-trace-state'],['--dd-step','nan'],['--dd-step','2'],
                  ['--no-levelwise'],['--dd-constraint-states','-1'],['--decoder','wrong']):
        result=run(*common,*extra,bpp)
        assert result.returncode != 0, extra
    raw=root/'native.bpp'; first=root/'first.bpseq'; second=root/'second.bpseq'
    native=['--decoder','dd','--loglevel','info','-r','0','-t','.2,.1']
    require_success(run(*native,'--bpp',raw,'-B',first,fixture))
    raw.write_text('\n'.join(x for x in raw.read_text().splitlines() if not x.startswith('#'))+'\n')
    require_success(run(*native,'-x','-B',second,raw))
    def pairs(path):
        return [line for line in path.read_text().splitlines() if line and line[0].isdigit()]
    assert pairs(first)==pairs(second), 'Native and imported BPP changed DD decoding'
    require_success(run('--decoder','dd','--loglevel','info','-r','1','-t','auto,auto',fixture))
print('DD CLI, native sparse BPP roundtrip and automatic threshold/refinement passed')
