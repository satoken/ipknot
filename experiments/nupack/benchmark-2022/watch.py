"""Persist audited progress and the complete benchmark package on completion."""
import collections
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

root = Path('/tmp/ipknot-nupack-evaluation-20261010')
output = Path('/tmp/ipknot-nupack-evaluation-report')
destination = Path('/home/sato-lab.org/satoken/.codex/worktrees/nupack-linear/ipknot/experiments/nupack/benchmark-2022')
scripts = Path(__file__).resolve().parent
refs = {r['id']: r for r in json.loads(gzip.open(root / 'references.json.gz', 'rt').read())}
records = {}
last_report = 0
while True:
    for path in (root / 'final').glob('*.json'):
        if path.name not in records:
            records[path.name] = json.loads(path.read_text())
    manifest = json.loads((root / 'final-manifest.json').read_text())
    done = 'completed_utc' in manifest
    if time.monotonic() - last_report > 300 or done:
        command = [sys.executable, str(scripts / 'report.py'), '--root', str(root), '--output', str(output)]
        if not done:
            command += ['--allow-incomplete', '--no-bootstrap']
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        for path in output.iterdir():
            if path.is_file():
                shutil.copy2(path, destination / path.name)
        shutil.copy2(root / 'final-manifest.json', destination / 'final-manifest.json')
        last_report = time.monotonic()
    values = list(records.values())
    outcome = collections.Counter((r['mode'], 'success' if r['exit_code'] == 0 else
        'timeout' if r['timeout'] else 'memory_limit' if 'bad_alloc' in r['stderr'] else 'error') for r in values)
    progress = dict(runs=len(values), total=47372,
        complete_sequences=sum(n == 4 for n in collections.Counter(r['id'] for r in values).values()),
        minimum_length=min((refs[r['id']]['length'] for r in values), default=None),
        outcomes={e: {status: n for (engine, status), n in outcome.items() if engine == e}
                  for e in ['nupack', 'lnupack', 'lpc', 'lpv']},
        ambiguous_successes={e: sum(r['exit_code'] == 0 and r['mode'] == e and
            bool(set(refs[r['id']]['sequence']) - set('ACGU')) for r in values)
            for e in ['nupack', 'lnupack', 'lpc', 'lpv']})
    errors = sum(v.get('error', 0) for v in progress['outcomes'].values())
    print(f"評価中: {progress['complete_sequences']:,}/11,843配列; "
          f"{len(values):,}/47,372 runs; 最短{progress['minimum_length']} nt; "
          f"NUPACK {progress['outcomes']['nupack']}; その他エラー {errors}", flush=True)
    if done:
        env = dict(__import__('os').environ, MPLCONFIGDIR='/tmp/ipknot-nupack-plot-cache',
                   PYTHONPATH='/tmp/ipknot-paper2022/python')
        subprocess.run(['/tmp/ipknot-nmr-python/cpython-3.12.14-linux-aarch64-gnu/bin/python3.12',
            str(scripts / 'plot.py'), '--input', str(output / 'per-rna.csv'),
            '--output', str(output / 'figures')], check=True, env=env)
        shutil.copytree(output / 'figures', destination / 'figures', dirs_exist_ok=True)
        hashes = {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sorted(destination.rglob('*')) if p.is_file() and p.name != 'artifacts.json'}
        (destination / 'artifacts.json').write_text(json.dumps(hashes, indent=2) + '\n')
        print('COMPLETE: audited report, predictions, figures and hashes saved', flush=True)
        break
    time.sleep(55)
