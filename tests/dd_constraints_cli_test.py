"""Verify ordinary fixed-structure constraints with a solver-free DD build."""
from pathlib import Path
import json
import math
import subprocess
import sys
import tempfile

binary = sys.argv[1]

def run(*args, success=True, message=None):
    thresholds = [] if '-t' in args else ['-t', '.01,.01']
    result = subprocess.run([binary, '--decoder', 'dd', '--loglevel', 'info',
                             '-r', '0', *thresholds, *map(str, args)],
                            capture_output=True, text=True, timeout=30)
    text = result.stdout + result.stderr
    assert (result.returncode == 0) == success, (args, text)
    if success:
        assert 'DD: iterations=' in text, text
    if message:
        assert message in text, text
    return text

def read_pairs(path):
    return {int(line.split()[0]) - 1: int(line.split()[2]) - 1
            for line in path.read_text().splitlines() if line and line[0].isdigit()}

with tempfile.TemporaryDirectory() as temporary:
    tmp = Path(temporary)
    output = tmp / 'prediction.bpseq'
    constraints = tmp / 'fixed.constraint'
    bpp = tmp / 'small.bpp'
    # Four stacked GC pairs forming a pseudoknot; all positions are independent
    # of posterior generation, so negative/fixed candidate behavior is visible.
    bpp.write_text('1 G 8:0.9\n2 G 7:0.9\n3 A\n4 C 11:0.9\n5 C 10:0.9\n6 A\n7 C\n8 C\n9 A\n10 G\n11 G\n')
    constraints.write_text('1 8\n2 7\n3 x\n4 11\n5 10\n6 x\n9 x\n')
    trace = tmp / 'trace.jsonl'
    text = run('-x', '-c', constraints, '-B', output,
               '--dd-trace', trace, '--dd-trace-state', bpp)
    pairs = read_pairs(output)
    assert {i:j for i,j in pairs.items() if i<j} == {0:7, 1:6, 3:10, 4:9}, pairs
    events = [json.loads(line) for line in trace.read_text().splitlines()]
    summary = next(event for event in events if event['event'] == 'summary')
    assert summary['constrained'] and math.isfinite(summary['upper_bound'])
    assert summary['upper_bound'] >= summary['lower_bound'] - 1e-8
    problem = next(event for event in events if event['event'] == 'problem')
    def audit(values, expected):
        assert len(values) == len(problem['variable_domains'])
        actual = 0
        for value, (coefficient, lower, upper, integer) in zip(values, problem['variable_domains']):
            assert lower - 1e-8 <= value <= upper + 1e-8
            assert not integer or abs(value - round(value)) < 1e-8
            actual += coefficient * value
        for bound, equality, terms in problem['rows']:
            total = sum(coefficient * values[column] for column, coefficient in terms)
            assert abs(total - bound) < 1e-8 if equality else total <= bound + 1e-8
        assert math.isclose(actual, expected, abs_tol=1e-8)
    for event in events:
        if event['event'] == 'iteration' and event['lower_bound'] is not None:
            audit(event['recovered'], event['lower_bound'])
        if event['event'] == 'repair':
            audit(event['selected'], summary['lower_bound'])
    run('-x', '-c', constraints, '--dd-exchange', '4', bpp, success=False,
        message='extra unconstrained bound/recovery options are unsupported')
    # Improved beam preserves fixed pseudoknots in both constraint modes,
    # including when a forced pair is missing from the posterior matrix.
    missing = tmp / 'missing.bpp'
    missing.write_text(bpp.read_text().replace('1 G 8:0.9', '1 G'))
    for mode in ('linear', 'full'):
        for evidence in (bpp, missing):
            run('-x', '-c', constraints, '-B', output,
                '--dd-dp', 'improved-beam', '--dd-constraints', mode, evidence)
            assert read_pairs(output) == pairs
    # Force one below-threshold pair with no stacking requirement. Other
    # candidates must stay free; a negative feasible objective remains valid.
    constraints.write_text('1 8\n3 x\n')
    run('-x', '-i', '-t', '1', '-c', constraints, '-B', output, bpp)
    assert read_pairs(output)[0] == 7 and read_pairs(output)[2] == -1
    # One-sided and either-sided forced pairing plus an unpaired position.
    constraints.write_text('1 <\n8 >\n3 x\n')
    run('-x', '-i', '-c', constraints, '-B', output, bpp)
    values = read_pairs(output)
    assert values[0] > 0 and 0 <= values[7] < 7 and values[2] == -1
    constraints.write_text('2 |\n')
    run('-x', '-i', '-c', constraints, '-B', output, bpp)
    assert read_pairs(output)[1] >= 0
    constraints.write_text('3 |\n')
    run('-x', '-i', '-c', constraints, bpp, success=False,
        message='no candidate partner')
    constraints.write_text('2 |\n')
    # The original all-candidate model remains an explicit diagnostic mode.
    run('-x', '-c', constraints, '--dd-constraints', 'full',
        '--dd-beam', '0', '--dd-constraint-states', '0', bpp)
    # Exact per-level Nussinov also supports fixed structures.
    run('-x', '-c', constraints, '--dd-dp', 'nussinov', '--dd-beam', '1', '--dd-trace', trace, bpp)
    summaries = [json.loads(line) for line in trace.read_text().splitlines()
                 if json.loads(line)['event'] == 'summary']
    assert summaries and all(s['exact_dp'] and s['dp_pruned_states'] == 0 for s in summaries)
    problems = [json.loads(line) for line in trace.read_text().splitlines()
                if json.loads(line)['event'] == 'problem']
    assert all(p['dp'] == 'nussinov' and p['beam'] == 0 for p in problems)
    run('-x', '-c', constraints, '--dd-dp', 'improved-beam', '--dd-beam', '100',
        '--dd-trace', trace, bpp)
    improved = [json.loads(line) for line in trace.read_text().splitlines()
                if json.loads(line)['event'] == 'problem']
    assert improved and all(p['dp'] == 'improved-beam' and p['beam'] == 100 for p in improved)
    run('-x', '-c', constraints, '--dd-dp', 'wrong', bpp,
        success=False, message='DD DP must be beam, improved-beam, nussinov or wavefront')
    run('-x', '-c', constraints, '--dd-constraint-recovery-every', '-1', bpp,
        success=False, message='recovery interval must be nonnegative')
    for option in ('--dd-crossing-beam', '--dd-witnesses', '--dd-constraint-states'):
        run('-x', '-c', constraints, option, '0', bpp, success=False,
            message='use --dd-constraints full')
    # Native refinement must preserve the original partial fixed structure.
    native = tmp / 'fixed.fa'
    native.write_text('>fixed\nAUAUAUAUGGCCAUAUAUAU\n')
    constraints.write_text('1 20\n3 x\n')
    run('-c', constraints, '-r', '1', '-B', output, '--dd-trace', trace, native)
    assert sum(json.loads(line)['event'] == 'summary' for line in trace.read_text().splitlines()) == 2
    values = read_pairs(output)
    assert values[0] == 19 and values[2] == -1
    run('-c', constraints, '-t', 'auto,auto', '-B', output, native)
    values = read_pairs(output)
    assert values[0] == 19 and values[2] == -1
print('Fixed structures, negative/missing candidates, exact DP, trace feasibility, refinement and validation passed')
