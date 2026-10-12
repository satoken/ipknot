#!/usr/bin/env python3
"""Bound the column-work heuristic from logged graph sizes; not a timing run."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from benchmark_dataset import load_helper, options


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--cpu', type=int, default=0)
    args = parser.parse_args()
    root = args.root.resolve()
    dataset = load_helper(args.helper)
    entries = dataset.select(dataset.load_index(root), 3)
    os.sched_setaffinity(0, {args.cpu})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    result = dict(binary_sha256=dataset.digest(root/'bin/candidate'), cpu=args.cpu,
                  scope='Diagnostic with info logging, excluded from timings.',
                  proof='Each level/mask has at most total logged M pairs. sum min(100,left+1)/n <= 100*M/n. An upper bound below 8192 guarantees serial fallback for every decode of that graph, including recovery masks.',
                  threshold=8192, rows=[])
    for entry in sorted(entries, key=lambda e: (-e['length'], e['id'])):
        for mode in ('lpc-serial', 'improved-serial'):
            prediction = root/'column-diagnostic.bpseq'
            command = [str(root/'bin/candidate'), *options(mode), '--loglevel', 'info',
                       '-B', str(prediction), str(root/'inputs'/(entry['id']+'.fa'))]
            run = subprocess.run(command, capture_output=True, text=True, env=env)
            if run.returncode:
                raise RuntimeError(run.stderr)
            sizes = [int(s) for s in re.findall(r'DD: iterations=\d+, pairs=(\d+),', run.stdout+run.stderr)]
            if not sizes:
                raise ValueError('missing graph sizes')
            seq, pairs = dataset.read_bpseq(prediction.read_text())
            if seq != entry['sequence']:
                raise ValueError('changed sequence')
            upper = 100*max(sizes)/entry['length']
            result['rows'].append(dict(id=entry['id'], mode=mode, length=entry['length'],
                                       graph_pair_counts=sizes, average_offers_upper_bound=upper,
                                       proves_serial_fallback=upper < 8192,
                                       input_sha256=entry['input_sha256'], command=command,
                                       prediction_sha256=hashlib.sha256(json.dumps(pairs,separators=(',', ':')).encode()).hexdigest()))
            prediction.unlink()
            dataset.write_json(root/'column-diagnostics.json', result)
            print(entry['id'], mode, upper, flush=True)


if __name__ == '__main__':
    main()
