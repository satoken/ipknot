"""Finish short-RNA linear accuracy separately while the timed exact runs continue.

This does not modify the original final phase or its executable. The 190-nt
cutoff is fixed before these additional predictions are computed. Timings from
this phase are kept separate from the original homogeneous-X925 comparison.
"""
import concurrent.futures
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

from dataset import load_index, write_json, digest
from evaluate import COMMON, init_worker, invoke

ROOT = Path('/tmp/ipknot-nupack-evaluation-20261010')
BINARY = '/tmp/ipknot-nupack-build/ipknot'
LAUNCHER = str(ROOT / 'measure-r0')
CPUS = [0, 1, 2, 3, 4, 10, 11, 12, 13, 14]
PHASE = 'accuracy-short'
CUTOFF = 190
ENGINES = ['lnupack', 'lpc', 'lpv']


def main():
    entries = [r for r in load_index(ROOT) if r['length'] <= CUTOFF]
    entries.sort(key=lambda r: (r['length'], r['id']))
    manifest = dict(phase=PHASE, engines=ENGINES, common_options=COMMON,
        refinement=dict(lnupack=0, lpc=1, lpv=1), beam=100,
        address_space_limit_gib=8, timeout_seconds=1800, n_sequences=len(entries),
        cutoff_nt=CUTOFF, sample_ids=[r['id'] for r in entries], cpus=CPUS,
        binary=BINARY, binary_sha256=digest(BINARY), launcher_sha256=digest(LAUNCHER),
        wrapper_sha256=digest(__file__), evaluate_sha256=digest(Path(__file__).with_name('evaluate.py')),
        dataset_helper_sha256=digest(Path(__file__).with_name('dataset.py')),
        reference_index_sha256=digest(ROOT / 'references.json.gz'),
        started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        selection='All official RNAs <=190 nt; input length cutoff fixed before this phase; no accuracy selection',
        measurement='Separate Cortex-A725 phase for full accuracy availability. Do not merge these timings with the original Cortex-X925 final-phase timings.',
        scheduling='One RNA and its three engines per pinned A725 worker; shortest first, deterministic shuffled engine order')
    mp = ROOT / (PHASE + '-manifest.json')
    if mp.exists():
        old = json.loads(mp.read_text())
        for key in manifest.keys() - {'started_utc'}:
            if old[key] != manifest[key]:
                raise ValueError('incompatible resume: ' + key)
        manifest = old
    else:
        write_json(mp, manifest)
    (ROOT / PHASE).mkdir(exist_ok=True)
    tasks = [(r, str(ROOT), BINARY, LAUNCHER, 100, 1800, PHASE, ENGINES, 8, 60) for r in entries]
    start = last = time.monotonic()
    completed = 0
    failures = {e: 0 for e in ENGINES}
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(CPUS), initializer=init_worker,
        initargs=(CPUS,), mp_context=multiprocessing.get_context('fork')) as pool:
        pending = [pool.submit(invoke, t) for t in tasks]
        for future in concurrent.futures.as_completed(pending):
            for engine, code in future.result():
                completed += 1
                failures[engine] += int(code != 0)
            now = time.monotonic()
            if now - last > 30:
                print(json.dumps(dict(phase=PHASE, completed=completed, total=3*len(entries),
                    elapsed_seconds=round(now-start, 1), failures=failures)), flush=True)
                last = now
    manifest.update(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        completed_runs=completed, elapsed_seconds=time.monotonic()-start, failures=failures)
    write_json(mp, manifest)
    print('COMPLETE', completed, failures, flush=True)


if __name__ == '__main__':
    main()
