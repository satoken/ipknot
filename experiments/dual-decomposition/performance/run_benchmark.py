"""Compare identical deterministic graphs; run without competing build jobs."""
import argparse
import json
import math
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--baseline', type=Path, required=True)
parser.add_argument('--candidate', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
cases = [
    ('exact_sparse', 'dp', 'sparse', 1500, 12, 0, 1),
    ('exact_dense', 'dp', 'dense', 500, 12, 0, 1),
    ('beam_sparse', 'dp', 'sparse', 4000, 12, 100, 1),
    ('dd_exact', 'dd', 'sparse', 1500, 12, 0, 2),
    ('dd_beam', 'dd', 'sparse', 4000, 12, 100, 2),
]
results = []
for name, *case in cases:
    measurements = []
    for binary, threads in [(args.baseline, 1), (args.candidate, 1)]:
        command = [str(binary), *map(str, case), '5', str(threads)]
        result = json.loads(subprocess.check_output(command, text=True, timeout=120))
        result['revision'] = 'baseline' if binary == args.baseline else 'candidate'
        measurements.append(result)
    if case[0] == 'dd':
        command = [str(args.candidate), *map(str, case), '5', '2']
        result = json.loads(subprocess.check_output(command, text=True, timeout=120))
        result['revision'] = 'candidate'
        measurements.append(result)
        serial, parallel = measurements[1:]
        for field in ['objective', 'upper_bound', 'iterations']:
            assert serial[field] == parallel[field], (name, field, serial, parallel)
    if case[0] == 'dp':
        assert math.isclose(measurements[0]['objective'], measurements[1]['objective'], abs_tol=1e-9)
    baseline_seconds = measurements[0]['median_seconds']
    for result in measurements:
        result['case'] = name
        result['speedup'] = baseline_seconds / result['median_seconds']
        print(f"{name:14s} {result['revision']:9s} threads={result['threads']} "
              f"{result['median_seconds']:.6f}s {result['speedup']:.2f}x", flush=True)
        results.append(result)
args.output.write_text(json.dumps(results, indent=2) + '\n')
