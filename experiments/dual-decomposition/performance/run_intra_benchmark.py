"""Compare single-DP engines against 1b12107, without competing build jobs."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--baseline', type=Path, required=True)
parser.add_argument('--candidate', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
cases = [
    ('exact_sparse', 'sparse', 'sparse', 1500, 12, 0, 0),
    ('exact_dense', 'sparse', 'dense', 1000, 12, 0, 0),
    ('exact_stacked', 'sparse', 'dense', 700, 12, 0, 1),
    ('beam_default', 'beam', 'sparse', 4000, 12, 100, 0),
    ('beam_wide', 'beam', 'sparse', 2000, 128, 512, 0),
    ('beam_dense', 'beam', 'dense', 1000, 12, 512, 0),
    ('improved_dense', 'improved', 'dense', 1000, 12, 512, 0),
    ('beam_stacked', 'beam', 'dense', 700, 12, 512, 1),
]
results = []
for name, engine, shape, n, degree, beam, stacked in cases:
    configurations = [('baseline', args.baseline, engine, 1), ('candidate', args.candidate, engine, 1)]
    if engine == 'sparse':
        configurations += [('candidate', args.candidate, 'wavefront', t) for t in (1, 2, 4)]
    else:
        configurations += [('candidate', args.candidate, engine, t) for t in (2, 4)]
    reference = None
    traces = {}
    for revision, binary, variant, threads in configurations:
        command = [str(binary), variant, shape, str(n), str(degree), str(beam), '5', str(threads), str(stacked)]
        result = json.loads(subprocess.check_output(command, text=True, timeout=120))
        result.update(case=name, revision=revision)
        if hasattr(os, 'sched_getaffinity'):
            result['cpu_affinity'] = sorted(os.sched_getaffinity(0))
        if reference is None:
            reference = result
        assert math.isclose(reference['objective'], result['objective'], rel_tol=1e-12, abs_tol=1e-9)
        if variant in traces:
            assert traces[variant] == (result['objective'], result['trace_hash'], result['pruned'])
        else:
            traces[variant] = (result['objective'], result['trace_hash'], result['pruned'])
        result['speedup_vs_baseline'] = reference['median_seconds'] / result['median_seconds']
        results.append(result)
        print(f"{name:16s} {revision:9s} {variant:9s} t={threads} "
              f"{result['median_seconds']*1000:8.2f} ms {result['speedup_vs_baseline']:.2f}x "
              f"workers={result.get('workers', 1)}", flush=True)
args.output.write_text(json.dumps(results, indent=2) + '\n')
