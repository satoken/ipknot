#!/usr/bin/env python3
"""Paired decoder performance study on the official 2022 single sequences."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import platform
import random
import subprocess
import time

COMMON = ['--decoder', 'dd', '--dd-max-iter', '50', '--dd-patience', '0',
          '--dd-crossing-beam', '100', '--dd-witnesses', '16',
          '-t', 'auto,auto', '-r', '1', '-n', '1']


def options(mode):
    engine, variant = mode.split('-', 1)
    dp = 'nussinov' if engine == 'contra' else 'beam'
    if variant.startswith('wave'):
        dp = 'wavefront'
    if engine == 'improved':
        dp = 'improved-beam'
    result = COMMON + ['-e', 'CONTRAfold' if engine == 'contra' else 'lpc',
                       '--beam-size', '100', '--dd-dp', dp, '--dd-beam', '100']
    if variant != 'before':
        result += ['--dd-threads', '2' if variant == 'level2' else '1',
                   '--dd-dp-threads', '4' if variant in ('wave4', 'inner4') else '1']
    return result


def load_helper(path):
    spec = importlib.util.spec_from_file_location('dataset', path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def worker_init(groups, helper_path):
    global GROUP, DATASET
    number = multiprocessing.current_process()._identity[-1] - 1
    GROUP = groups[number % len(groups)]
    DATASET = load_helper(helper_path)


def measure(payload):
    entry, config, signature = payload
    root = Path(config['root'])
    fasta = root / 'inputs' / (entry['id'] + '.fa')
    if DATASET.digest(fasta) != entry['input_sha256']:
        raise ValueError('changed input: ' + str(fasta))
    outcomes = []
    env = {k: v for k, v in os.environ.items()
           if not k.startswith('IPKNOT_STUDY_') and k != 'IPKNOT_HIGHS_LOG'}
    env.update(OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
               OMP_PROC_BIND='FALSE', OMP_WAIT_POLICY=config['wait_policy'])
    if config['wait_policy'] == 'DEFAULT':
        env.pop('OMP_WAIT_POLICY', None)
        env.pop('GOMP_SPINCOUNT', None)
    # Bind the process explicitly; do not inherit site-specific runtime placement.
    env.pop('OMP_PLACES', None)
    env.pop('GOMP_CPU_AFFINITY', None)
    for rep in range(config['repetitions']):
        order = list(config['modes'])
        random.Random(int(entry['id'], 16) + rep * 7919).shuffle(order)
        for mode in order:
            threads = 4 if mode.endswith(('wave4', 'inner4')) else 2 if mode.endswith('level2') else 1
            os.sched_setaffinity(0, set(GROUP[:threads]))
            env['OMP_NUM_THREADS'] = str(threads)
            target = root / config['phase'] / f"{entry['id']}-{mode}-{rep}.json"
            prediction = target.with_suffix('.bpseq')
            resources = target.with_suffix('.resources')
            binary = config['baseline'] if mode.endswith('before') else config['candidate']
            cmd = [binary, *options(mode), '-B', str(prediction), str(fasta)]
            if target.exists():
                row = json.loads(target.read_text())
                if (row['configuration_sha256'] != signature or row['command'] != cmd or
                        row['input_sha256'] != entry['input_sha256']):
                    raise ValueError('incompatible saved run: ' + str(target))
                outcomes.append((mode, row['exit_code']))
                continue
            run = subprocess.run([config['launcher'], str(resources), str(config['timeout']),
                                  str(config['memory_gib']), *cmd],
                                 capture_output=True, text=True, env=env)
            if not resources.exists():
                raise RuntimeError('resource launcher failed: ' + run.stderr)
            row = json.loads(resources.read_text())
            if run.returncode != row['exit_code']:
                raise RuntimeError('launcher/child exit-code mismatch')
            row.update(id=entry['id'], mode=mode, phase=config['phase'], repetition=rep,
                       cpus=sorted(os.sched_getaffinity(0)), command=cmd,
                       input_sha256=entry['input_sha256'], configuration_sha256=signature,
                       stderr=run.stderr, stdout_sha256=hashlib.sha256(run.stdout.encode()).hexdigest())
            if row['exit_code'] == 0:
                seq, pairs = DATASET.read_bpseq(prediction.read_text())
                if seq != entry['sequence']:
                    raise ValueError('prediction sequence mismatch: ' + entry['id'])
                pp = set(pairs)
                if any((i-1, j+1) not in pp and (i+1, j-1) not in pp for i, j in pairs):
                    raise ValueError('isolated pair in output: ' + entry['id'])
                if any(seq[i-1] + seq[j-1] not in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG') for i, j in pairs):
                    raise ValueError('noncanonical pair in output: ' + entry['id'])
                row.update(pairs=pairs, accuracy=DATASET.metrics(entry['pairs'], pairs),
                           crossing_accuracy=DATASET.metrics(DATASET.crossing_pairs(entry['pairs']),
                                                             DATASET.crossing_pairs(pairs)))
            DATASET.write_json(target, row)
            prediction.unlink(missing_ok=True)
            resources.unlink()
            outcomes.append((mode, row['exit_code']))
    return outcomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--launcher', type=Path, required=True)
    parser.add_argument('--phase', required=True)
    parser.add_argument('--modes', nargs='+', required=True)
    parser.add_argument('--cpu-groups', default='5;6;7;8;9;15;16;17;18;19')
    parser.add_argument('--sample-per-stratum', type=int, default=0)
    parser.add_argument('--repetitions', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=1800)
    parser.add_argument('--memory-gib', type=int, default=8)
    parser.add_argument('--wait-policy', choices=['ACTIVE', 'PASSIVE', 'DEFAULT'], default='DEFAULT')
    args = parser.parse_args()
    for name in ('root', 'helper', 'baseline', 'candidate', 'launcher'):
        setattr(args, name, str(getattr(args, name).resolve()))
    dataset = load_helper(args.helper)
    root = Path(args.root)
    entries = dataset.select(dataset.load_index(root), args.sample_per_stratum)
    groups = [[int(c) for c in g.split(',')] for g in args.cpu_groups.split(';')]
    flat = [c for group in groups for c in group]
    if len(flat) != len(set(flat)) or not set(flat) <= os.sched_getaffinity(0):
        raise ValueError('unavailable/overlapping CPU groups')
    needed = max(4 if m.endswith(('wave4', 'inner4')) else 2 if m.endswith('level2') else 1 for m in args.modes)
    if min(map(len, groups)) < needed:
        raise ValueError('not enough cores per worker')
    cpu_info = {str(c): {
        'midr': Path(f'/sys/devices/system/cpu/cpu{c}/regs/identification/midr_el1').read_text().strip(),
        'max_khz': int(Path(f'/sys/devices/system/cpu/cpu{c}/cpufreq/cpuinfo_max_freq').read_text()),
    } for c in flat}
    if len({v['midr'] for v in cpu_info.values()}) != 1:
        raise ValueError('CPUs are not homogeneous')
    config = vars(args)
    definition = dict(config)
    definition.update(options={m: options(m) for m in args.modes},
                      n_sequences=len(entries), sample_ids=sorted(e['id'] for e in entries),
                      cpu_groups=groups, cpu_info=cpu_info,
                      hashes={k: dataset.digest(v) for k, v in {
                          'baseline': args.baseline, 'candidate': args.candidate,
                          'launcher': args.launcher, 'helper': args.helper,
                          'script': __file__, 'references': root/'references.json.gz'}.items()})
    signature = hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest()
    manifest = dict(**definition, configuration_sha256=signature, uname=list(platform.uname()),
                    started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    measurement='Fresh process: BPP, auto thresholds, refinement; no BPP cache/trace. CLOCK_MONOTONIC wall, wait4 CPU and peak RSS.',
                    scheduling='Longest first; same RNA uses the same CPU group. Serial uses group[0], parallel the first N cores. Hash-shuffled condition order.',
                    selection='All official entries; timing subset uses the established timing2022/member hash, three per dataset/size/PK stratum, independent of accuracy.')
    path = root / (args.phase + '-manifest.json')
    if path.exists():
        old = json.loads(path.read_text())
        if old['configuration_sha256'] != signature:
            raise ValueError('resume configuration changed')
        manifest = old
    else:
        dataset.write_json(path, manifest)
    (root / args.phase).mkdir(exist_ok=True)
    entries.sort(key=lambda e: (-e['length'], e['id']))
    tasks = [(e, config, signature) for e in entries]
    total = len(entries) * len(args.modes) * args.repetitions
    completed = 0
    failures = {m: 0 for m in args.modes}
    start = last = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=len(groups), initializer=worker_init, initargs=(groups, args.helper),
            mp_context=multiprocessing.get_context('fork')) as pool:
        pending = {pool.submit(measure, t) for t in tasks}
        while pending:
            done, pending = concurrent.futures.wait(pending, timeout=30,
                                                    return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                for mode, code in future.result():
                    completed += 1
                    failures[mode] += int(code != 0)
            now = time.monotonic()
            if now-last >= 30 or not pending:
                progress = dict(phase=args.phase, completed=completed, total=total,
                                failures=failures, elapsed_seconds=round(now-start, 1))
                dataset.write_json(root / (args.phase + '-progress.json'), progress)
                print(json.dumps(progress), flush=True)
                last = now
    for name in ('baseline', 'candidate', 'launcher', 'helper'):
        if dataset.digest(config[name]) != manifest['hashes'][name]:
            raise ValueError(name + ' changed during measurement')
    manifest.update(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    invocation_elapsed_seconds=time.monotonic()-start, completed_runs=completed,
                    failures=failures)
    dataset.write_json(path, manifest)
    print('COMPLETE', completed, failures, flush=True)
    if any(failures.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
