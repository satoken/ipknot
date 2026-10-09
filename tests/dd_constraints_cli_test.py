"""Verify fixed/NMR constraints with a solver-free DD build and real witnesses."""
from pathlib import Path
import json
import math
import random
import subprocess
import sys
import tempfile

binary = sys.argv[1]
root = Path(__file__).resolve().parents[1]
data = root / 'tests/data'

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

def count(pairs, seq, kind):
    return sum(1 for left, right in pairs.items()
               if left < right and sorted((seq[left], seq[right])) == sorted(kind))

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
    text = run('-x', '-c', constraints, '-B', output, '--base-pairs', 'GC=4',
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
    # Crossing/blocks positive and negative PK terms, including normalized
    # continuous score columns, preserve fixed pseudoknot structure.
    for intercept in ('8', '-8'):
        run('-x', '-c', constraints, '-B', output, '--pk-h-intercept', intercept,
            '--pk-crossing-normalize', bpp)
        assert read_pairs(output) == pairs
        for mode in ('linear','full'):
            text=run('-x','-c',constraints,'-B',output,'--dd-dp','improved-beam',
                     '--dd-constraints',mode,'--pk-h-formulation','projected',
                     '--pk-h-intercept',intercept,bpp)
            assert read_pairs(output) == pairs
    # Learned projection abstains from confidence for a block containing a
    # forced pair absent from BPP; its optional shape term remains defined.
    missing = tmp/'missing.bpp'
    missing.write_text(bpp.read_text().replace('1 G 8:0.9', '1 G'))
    learned = tmp/'model.txt'
    names = ('bias anchor_support weak_support weak_min_support support_product '
             'support_gap weak_specificity anchor_specificity_product min_length '
             'weak_length outer_loop middle_loop').split()
    learned.write_text('IPKNOT_PK_LINEAR_V1\n' + ''.join(
        f'{name} {1 if name == "bias" else 0}\n' for name in names))
    for mode in ('linear','full'):
        for shape in ([], ['--pk-hybrid-shape','--pk-h-intercept','-8']):
            run('-x','-c',constraints,'-B',output,'--dd-constraints',mode,
                '--pk-h-formulation','projected','--pk-learned-model',learned,
                '--pk-learned-scale','.05',*shape,missing)
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
    # Auxiliary import must use stack observations as well as base-pair counts.
    run('-x', '--stack-constraint', 'GC GC', bpp)
    constraints.write_text('2 |\n')
    # The original all-candidate model remains an explicit diagnostic mode.
    run('-x', '-c', constraints, '--dd-constraints', 'full',
        '--dd-beam', '0', '--dd-constraint-states', '0', bpp)
    # Exact per-level Nussinov may be combined with bounded NMR candidates.
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
        success=False, message='DD DP must be beam, improved-beam or nussinov')
    run('-x', '-c', constraints, '--dd-constraint-recovery-every', '-1', bpp,
        success=False, message='recovery interval must be nonnegative')
    for option in ('--dd-crossing-beam', '--dd-witnesses', '--dd-constraint-states'):
        run('-x', '-c', constraints, option, '0', bpp, success=False,
            message='use --dd-constraints full')
    # Empty retained candidates still pay soft count penalties.
    run('-x', '-t', '1', '--base-pairs', 'GC=1', '--nmr-soft', bpp,
        message='Soft NMR base-pair constraint GC violated')
    # Canonical exact versus lower-bound totals, and both soft penalty outcomes.
    run('--base-pairs', 'AU=1', '-B', output, data/'nmr_constraints.fa')
    native_seq = 'AUAUAUAUGGCCAUAUAUAU'
    assert count(read_pairs(output), native_seq, 'AU') == 1
    run('--base-pairs', 'AU=1', '--nmr-count-mode', 'lower-bound', '-B', output,
        data/'nmr_constraints.fa')
    assert count(read_pairs(output), native_seq, 'AU') > 1
    run('--base-pairs', 'UU=1', '--without-canonical-neighbor', '-B', output,
        data/'nmr_constraints.fa')
    assert count(read_pairs(output), native_seq, 'UU') == 1
    run('--base-pairs', 'GC=99', data/'nmr_constraints.fa', success=False,
        message='infeasible')
    run('--base-pairs', 'GC=99', '--nmr-soft', data/'nmr_constraints.fa',
        message='Soft NMR base-pair constraint GC violated: expected 99')
    run('--base-pairs', 'UU=1', '--without-canonical-neighbor', '--nmr-soft',
        data/'nmr_constraints.fa', message='Soft NMR base-pair constraint UU=1 satisfied')
    run('--base-pairs', 'UU=1', '--without-canonical-neighbor', '--nmr-soft',
        '--nmr-count-penalty', '.0001', data/'nmr_constraints.fa',
        message='Soft NMR base-pair constraint UU violated')
    run('--base-pairs', 'AU=0', '--nmr-count-mode', 'lower-bound', '--nmr-soft',
        data/'nmr_constraints.fa', message='Soft NMR base-pair constraint AU=0 satisfied')
    # Direct, compact and left/right one-base bulge observations.
    run('--stack-constraint', 'AU AU', data/'nmr_constraints.fa', message='All stack constraints are satisfied')
    run('--stack-constraint', '(GU)-G-U', data/'nmr_compact_stack.fa', message='All stack constraints are satisfied')
    for name, expected in [('nmr_single_bulge', {0:10, 2:9}),
                           ('nmr_single_bulge_right', {0:10, 1:8})]:
        run('--stack-constraint', 'GC AU', '-B', output, data/f'{name}.fa')
        assert {i:j for i,j in read_pairs(output).items() if i<j} == expected
    run('--stack-constraint', 'GC AU', '--nmr-bulge-mode', 'none', data/'nmr_single_bulge.fa', success=False)
    run('--stack-constraint', 'GC AU', '--nmr-bulge-mode', 'fallback',
        data/'nmr_single_bulge.fa', message='All stack constraints are satisfied')
    duplicate = ['--stack-constraint','GC AU','--stack-constraint','GC AU']
    run(*duplicate, data/'nmr_single_bulge.fa', success=False)
    run(*duplicate, '--nmr-allow-shared-stack-pairs', data/'nmr_single_bulge.fa')
    run('--stack-constraint', 'UU CC', '-t', '0,0', '--nmr-soft',
        data/'nmr_coaxial_multibranch.fa', message='Soft NMR stack constraint 1 violated')
    # All three coaxial contexts, bulges, pair sharing and helix-face capacity.
    for kind in ('multibranch', 'children', 'last_child'):
        fixture = data/f'nmr_coaxial_bulge_{kind}'
        common = ['-t','0,0','--nmr-bulge-mode','all','--coaxial-stacking',
                  '--stack-constraint','UU CC','--stack-constraint','GC AU',
                  '-c',fixture.with_suffix('.constraint')]
        for flags in ([], ['--nmr-soft'], ['--nmr-allow-shared-stack-pairs'],
                      ['--nmr-soft','--nmr-allow-shared-stack-pairs']):
            run(*common, *flags, fixture.with_suffix('.fa'), message='All stack constraints are satisfied')
        for mode in ("linear", "full"):
            run(*common, "--dd-constraints", mode,
                "--nmr-coaxial-no-adjacent-bulge", fixture.with_suffix(".fa"),
                message="All stack constraints are satisfied")
        faces = [*common, '--stack-constraint','UU CC','--nmr-allow-shared-stack-pairs']
        run(*faces, fixture.with_suffix('.fa'), success=False, message='infeasible')
        text = run(*faces, '--nmr-soft', fixture.with_suffix('.fa'))
        assert ('Soft NMR stack constraint 1 violated' in text or
                'Soft NMR stack constraint 3 violated' in text), text
    run("--nmr-coaxial-no-adjacent-bulge", data/"nmr_constraints.fa",
        success=False, message="requires --coaxial-stacking")
    # A fixed bulge at an observed coaxial terminus is permitted by default
    # and incurs a soft violation only when the local straight-stem prior is on.
    adjacent = ["-t", "0,0", "--coaxial-stacking", "--nmr-bulge-mode", "all",
                "--nmr-allow-shared-stack-pairs", "--stack-constraint", "UU CC",
                "--stack-constraint", "AC UU", "--stack-constraint", "GC AU",
                "-c", data/"nmr_coaxial_adjacent_bulge.constraint"]
    for mode in ("linear", "full"):
        run(*adjacent, "--dd-constraints", mode, data/"nmr_coaxial_bulge_multibranch.fa",
            message="All stack constraints are satisfied")
        run(*adjacent, "--dd-constraints", mode, "--nmr-coaxial-no-adjacent-bulge",
            data/"nmr_coaxial_bulge_multibranch.fa", success=False)
        run(*adjacent, "--dd-constraints", mode, "--nmr-coaxial-no-adjacent-bulge", "--nmr-soft",
            data/"nmr_coaxial_bulge_multibranch.fa", message="Soft NMR stack constraint 1 violated")
    # Fallback semantics use global direct-pattern existence, even when the
    # bounded candidate pool is empty. Compare the word index with a separate
    # exhaustive scan over both sequence strands and pattern orientations.
    rng = random.Random(91)
    for trial in range(40):
        seq = ''.join(rng.choice('ACGU') for _ in range(28))
        width = 2 + trial % 3
        kinds = [rng.choice(('AU','GC','GU','UU','AA','CC','GG','AC','AG','CU')) for _ in range(width)]
        expected = any(all(sorted((seq[i+d],seq[j-d])) == sorted(kind)
                           for d,kind in enumerate(pattern))
                       for pattern in (kinds,kinds[::-1])
                       for i in range(len(seq))
                       for j in range(i+4+2*(width-1),len(seq))
                       if i+width <= len(seq))
        bpp.write_text(''.join(f'{i+1} {base}\n' for i,base in enumerate(seq)))
        constraints.write_text(''.join(f'{i+1} x\n' for i in range(len(seq))))
        text = run('-x','-c',constraints,'--nmr-soft','--nmr-bulge-mode','fallback',
                   '--dd-nmr-pair-beam','1','--dd-nmr-witnesses','1',
                   '--stack-constraint',' '.join(kinds),bpp)
        assert ('uses 0 direct instance' in text) == expected, (seq,kinds,expected,text)
    # Native refinement preserves the original partial constraints and an
    # automatic threshold search skips infeasible threshold combinations.
    constraints.write_text('1 20\n3 x\n')
    run('-c', constraints, '-r', '1', '-B', output, '--dd-trace', trace, data/'nmr_constraints.fa')
    assert sum(json.loads(line)['event'] == 'summary' for line in trace.read_text().splitlines()) == 2
    values = read_pairs(output)
    assert values[0] == 19 and values[2] == -1
    run('-t','auto,auto','--base-pairs','AU=1',data/'nmr_constraints.fa')
    run('-t','auto,auto','--base-pairs','GC=99',data/'nmr_constraints.fa', success=False,
        message='No feasible structure for any DD threshold combination')
print('DD fixed structure, NMR counts, stacks, bulges, coaxial capacity and refinement passed')
