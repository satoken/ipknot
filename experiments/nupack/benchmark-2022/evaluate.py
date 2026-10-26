"""Compare LinearNUPACK, LPC and LPV on the complete official 2022 dataset.

All engines use the same binary and DD options; failures are retained. Each
sequence stays on one CPU and engine order is shuffled independently of scores.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import random
import subprocess
import time

from dataset import digest, load_index, metrics, crossing_pairs, read_bpseq, write_json

ENGINES = ['lnupack', 'lpc', 'lpv', 'nupack']
COMMON = ['--decoder', 'dd', '--dd-dp', 'beam', '--dd-beam', '100',
          '--dd-max-iter', '50', '--dd-patience', '0', '--dd-crossing-beam', '100',
          '--dd-witnesses', '16', '-t', 'auto,auto', '-n', '1']


def init_worker(cpus):
    identity = multiprocessing.current_process()._identity[-1]
    os.sched_setaffinity(0, {cpus[(identity - 1) % len(cpus)]})


def invoke(payload):
    entry, root, binary, launcher, beam, timeout, phase, engines, memory_gib, exact_timeout = payload
    root = Path(root)
    order = list(engines)
    random.Random(int(entry['id'], 16)).shuffle(order)
    results = []
    fasta = root / 'inputs' / (entry['id'] + '.fa')
    if digest(fasta) != entry['input_sha256']:
        raise ValueError('changed FASTA: ' + entry['id'])
    for engine in order:
        target = root / phase / f"{entry['id']}-{engine}.json"
        command = [binary, *COMMON, '--beam-size', str(beam), '-e', engine,
                   '-r', '0' if engine in ('lnupack', 'nupack') else '1',
                   '-B', str(target.with_suffix('.bpseq')), str(fasta)]
        if target.exists():
            old = json.loads(target.read_text())
            if old['command'] != command or old['input_sha256'] != entry['input_sha256']:
                raise ValueError('incompatible saved run: ' + str(target))
            results.append((engine, old['exit_code']))
            continue
        resource = target.with_suffix('.resources')
        env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
               'MKL_NUM_THREADS': '1'}
        limit = exact_timeout if engine == 'nupack' else timeout
        run = subprocess.run([launcher, str(resource), str(limit), str(memory_gib), *command],
                             capture_output=True, text=True, env=env)
        if not resource.exists():
            raise RuntimeError('measurement launcher failed: ' + run.stderr)
        row = json.loads(resource.read_text())
        row.update(id=entry['id'], mode=engine, repetition=0, phase=phase,
                   input_sha256=entry['input_sha256'], command=command,
                   cpu=next(iter(os.sched_getaffinity(0))), stderr=run.stderr,
                   stdout_sha256=hashlib.sha256(run.stdout.encode()).hexdigest())
        if row['exit_code'] == 0:
            seq, pairs = read_bpseq(target.with_suffix('.bpseq').read_text())
            if seq != entry['sequence']:
                raise ValueError('prediction sequence mismatch: ' + entry['id'])
            row.update(pairs=pairs, accuracy=metrics(entry['pairs'], pairs),
                       crossing_accuracy=metrics(crossing_pairs(entry['pairs']),
                                                 crossing_pairs(pairs)))
        write_json(target, row)
        target.with_suffix('.bpseq').unlink(missing_ok=True)
        resource.unlink()
        results.append((engine, row['exit_code']))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--launcher', type=Path, required=True)
    parser.add_argument('--beam', type=int, default=100)
    parser.add_argument('--engines', nargs='+', choices=ENGINES, default=ENGINES[:3])
    parser.add_argument('--memory-gib', type=int, default=8)
    parser.add_argument('--cpus', default='5,6,7,8,9,15,16,17,18,19')
    parser.add_argument('--timeout', type=int, default=1800)
    parser.add_argument('--exact-timeout', type=int, default=60)
    parser.add_argument('--phase', default='full100')
    parser.add_argument('--sample-per-stratum', type=int, default=0)
    args = parser.parse_args()
    from dataset import select
    entries = select(load_index(args.root), args.sample_per_stratum)
    cpus = list(map(int, args.cpus.split(',')))
    manifest = dict(phase=args.phase, engines=args.engines, common_options=COMMON,
                    refinement={e: 0 if e in ('lnupack', 'nupack') else 1 for e in args.engines},
                    address_space_limit_gib=args.memory_gib,
                    beam=args.beam, n_sequences=len(entries), workers=len(cpus), cpus=cpus,
                    timeout_seconds={e: args.exact_timeout if e == 'nupack' else args.timeout for e in args.engines},
                    binary=str(args.binary),
                    binary_sha256=digest(args.binary), launcher_sha256=digest(args.launcher),
                    script_sha256=digest(__file__), dataset_helper_sha256=digest(Path(__file__).with_name('dataset.py')),
                    reference_index_sha256=digest(args.root / 'references.json.gz'),
                    sample_ids=[r['id'] for r in entries], sample_per_stratum=args.sample_per_stratum,
                    uname=list(platform.uname()), started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    measurement='CLOCK_MONOTONIC wall, wait4 child CPU and peak RSS (KiB); includes posterior, DD and refinement; single run per condition',
                    scheduling='One sequence and its engines per pinned worker; deterministic shuffled engine order; longest first',
                    selection='All official entries, including ambiguous bases; no accuracy selection, parameter fitting or BPP cache')
    manifest_path = args.root / (args.phase + '-manifest.json')
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        for key in manifest.keys() - {'started_utc'}:
            if old[key] != manifest[key]:
                raise ValueError('resume manifest changed: ' + key)
        manifest = old
    else:
        write_json(manifest_path, manifest)
    (args.root / args.phase).mkdir(exist_ok=True)
    entries.sort(key=lambda r: (-r['length'], r['id']))
    tasks = [(r, str(args.root), str(args.binary), str(args.launcher), args.beam,
              args.timeout, args.phase, args.engines, args.memory_gib, args.exact_timeout) for r in entries]
    completed = 0
    failures = {e: 0 for e in args.engines}
    start = last = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(cpus), initializer=init_worker,
            initargs=(cpus,), mp_context=multiprocessing.get_context('fork')) as pool:
        futures = [pool.submit(invoke, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            for engine, code in future.result():
                completed += 1
                failures[engine] += int(code != 0)
            now = time.monotonic()
            if now - last > 30:
                print(json.dumps(dict(completed=completed, total=len(tasks) * len(args.engines),
                                      failures=failures, elapsed_seconds=round(now - start, 1))), flush=True)
                last = now
    manifest.update(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    elapsed_seconds=time.monotonic() - start, completed_runs=completed, failures=failures)
    write_json(manifest_path, manifest)
    print('COMPLETE', completed, failures, flush=True)


if __name__ == '__main__':
    main()
