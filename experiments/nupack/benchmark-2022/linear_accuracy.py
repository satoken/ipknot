"""Audit full linear-engine accuracy without merging measurements across CPUs.

RNAs <=190 nt come from accuracy-short (A725). Longer RNAs come from the
original final phase (X925). This length rule precedes accuracy evaluation.
Both phases use the identical frozen binary, decoder and per-engine options.
The original final phase remains the source of all performance comparisons.
"""
import collections
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from dataset import load_index, write_json, digest

ROOT = Path('/tmp/ipknot-nupack-evaluation-20261010')
OUTPUT = Path('/tmp/ipknot-nupack-evaluation-report')
DEST = Path('/home/sato-lab.org/satoken/.codex/worktrees/nupack-linear/ipknot/experiments/nupack/benchmark-2022')
AUDIT = Path('/tmp/ipknot-nupack-accuracy-audit')
CUTOFF = 190
ENGINES = ['lnupack', 'lpc', 'lpv']
LABELS = {'lnupack': 'LinearNUPACK (r=0)', 'lpc': 'LPC (r=1)', 'lpv': 'LPV (r=1)'}


def main():
    refs = load_index(ROOT)
    sources = [(r, engine, ROOT / ('accuracy-short' if r['length'] <= CUTOFF else 'final') /
                f"{r['id']}-{engine}.json") for r in refs for engine in ENGINES]
    while True:
        missing = sum(not p.exists() for _, _, p in sources)
        print(f'線形3方式・全件精度: {len(sources)-missing:,}/{len(sources):,} runs', flush=True)
        if not missing:
            break
        time.sleep(55)
    short_manifest = json.loads((ROOT / 'accuracy-short-manifest.json').read_text())
    final_manifest = json.loads((ROOT / 'final-manifest.json').read_text())
    for key in ['binary_sha256', 'reference_index_sha256', 'launcher_sha256', 'beam', 'common_options']:
        if short_manifest[key] != final_manifest[key]:
            raise ValueError('phase settings disagree: ' + key)
    for e in ENGINES:
        if short_manifest['refinement'][e] != final_manifest['refinement'][e]:
            raise ValueError('different refinement')
    AUDIT.mkdir(exist_ok=True)
    (AUDIT / 'final').mkdir(exist_ok=True)
    (AUDIT / 'references.json.gz').unlink(missing_ok=True)
    (AUDIT / 'references.json.gz').symlink_to(ROOT / 'references.json.gz')
    manifests = {'accuracy-short': short_manifest, 'final': final_manifest}
    same_predictions = 0
    for ref, engine, path in sources:
        row = json.loads(path.read_text())
        if row['exit_code'] != 0:
            raise ValueError('linear accuracy run failed: ' + str(path))
        if ref['length'] <= CUTOFF:
            other = ROOT / 'final' / path.name
            if other.exists():
                baseline = json.loads(other.read_text())
                if baseline['exit_code'] == 0:
                    if baseline['pairs'] != row['pairs']:
                        raise ValueError('predictions differ between CPUs: ' + path.name)
                    same_predictions += 1
        target = AUDIT / 'final' / path.name
        target.unlink(missing_ok=True)
        os.link(path, target)
    write_json(AUDIT / 'final-manifest.json', dict(phase='final', engines=ENGINES,
        sample_ids=[r['id'] for r in refs], source_manifests=manifests,
        selection='All <=190 nt: accuracy-short; all >190 nt: final; cutoff chosen before this auxiliary phase',
        scope='Accuracy only; mixed CPU classes. Use the original final phase for runtime and memory.'))
    audited = AUDIT / 'report'
    subprocess.run([sys.executable, str(Path(__file__).with_name('report.py')), '--root', str(AUDIT),
        '--output', str(audited)], check=True, stdout=subprocess.DEVNULL)
    result = json.loads((audited / 'summary.json').read_text())
    if not result['complete'] or result['common_linear_n'] != 11843:
        raise ValueError('incomplete full linear comparison')
    groups = {}
    for group, engines in result['groups'].items():
        if not group.startswith('common_linear'):
            continue
        groups[group] = {e: {k: engines[e][k] for k in ['attempts', 'successes', 'accuracy',
            'crossing_accuracy_pk_only', 'pk_sequences']} for e in ENGINES}
    comparisons = {}
    for name, item in result['comparisons'].items():
        if item is not None:
            comparisons[name] = {k: item[k] for k in ['n', 'macro_f1_delta_b_minus_a',
                'bootstrap_95_percentile', 'better', 'equal', 'worse']}
    provenance = dict(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        scope='Full-dataset accuracy only. Timings from different CPU classes are excluded.',
        cutoff_nt=CUTOFF, selection='All <=190 nt: accuracy-short/A725; all >190 nt: final/X925',
        source_manifests=manifests, wrapper_sha256=digest(Path(__file__).with_name('complete_linear.py')),
        report_script_sha256=digest(__file__), audit_script_sha256=digest(Path(__file__).with_name('report.py')),
        overlapping_final_short_predictions_identical=same_predictions,
        ambiguous_sequences=sum(bool(set(r['sequence'])-set('ACGU')) for r in refs),
        frozen_binary_sha256=final_manifest['binary_sha256'])
    summary = dict(complete=True, n_sequences=11843, n_runs=len(sources), groups=groups,
        comparisons=comparisons, provenance=provenance, metric_convention=result['metric_convention'])
    write_json(OUTPUT / 'linear-accuracy.json', summary)
    with gzip.open(OUTPUT / 'linear-accuracy-raw.jsonl.gz', 'wt') as f:
        for _, _, path in sources:
            f.write(json.dumps(json.loads(path.read_text()), separators=(',', ':')) + '\n')
    lines = ['# 線形3方式・全11,843配列の精度比較', '',
        '状態: 全件完了。公式2022年IPknot単一配列セット11,843配列、各方式とも全件正常終了。',
        'LinearNUPACKはrefinementなし、LPC・LPVはrefinement 1。全方式beam 100、同一DD設定。曖昧塩基を含む41配列もすべて含む。', '',
        'これは精度を先行集計した結果です。190 nt以下は別のCortex-A725コア、それより長い配列は元のCortex-X925測定から採用し、実行バイナリと設定は同一です。長さによる採用規則は追加評価の前に固定しました。CPUの異なる実行時間・メモリは混合しません。',
        '通常NUPACKの60秒・8 GiB評価と、全方式をCortex-X925に固定した速度比較は [RESULTS.md](RESULTS.md) に順次反映されます。', '']
    for key, title in [('common_linear', '全配列'), ('common_linear/short', '短鎖（150 nt以下）'),
        ('common_linear/medium', '中鎖（151–500 nt）'), ('common_linear/long', '長鎖（500 nt超）'),
        ('common_linear/Rfam14.5', 'Rfam 14.5'), ('common_linear/bpRNA-1m', 'bpRNA-1m')]:
        lines += ['## ' + title, '', '| 方式 | 配列数 | Macro F1 | PPV | Sensitivity | PK F1 |',
                  '|---|---:|---:|---:|---:|---:|']
        for e in ENGINES:
            r = groups[key][e]; a = r['accuracy']; pk = r['crossing_accuracy_pk_only']
            lines.append(f"| {LABELS[e]} | {r['successes']} | {a['macro_f1']:.6f} | {a['macro_ppv']:.6f} | {a['macro_sen']:.6f} | {pk['macro_f1']:.6f} |")
        lines.append('')
    lines += ['## LinearNUPACKとの差（配列ごとのpaired bootstrap）', '',
        '| 比較元 | Macro F1差 (LinearNUPACK − 比較元) | 95% CI | LinearNUPACKが上 / 同じ / 下 |',
        '|---|---:|---|---:|']
    for name in ['lpc_vs_lnupack', 'lpv_vs_lnupack']:
        c = comparisons[name]; ci = c['bootstrap_95_percentile']
        lines.append(f"| {name.split('_vs_')[0].upper()} | {c['macro_f1_delta_b_minus_a']:+.6f} | [{ci[0]:+.6f}, {ci[1]:+.6f}] | {c['better']} / {c['equal']} / {c['worse']} |")
    lines += ['', result['metric_convention'], '',
        '保存予測は全件、参照ハッシュ・塩基対座標・塩基種・DDのstack条件・独立再計算したスコアを検証済み。rawは `linear-accuracy-raw.jsonl.gz`、条件と集計は `linear-accuracy.json`。速度比較には元の `final` phaseのみを使用します。', '']
    (OUTPUT / 'LINEAR_ACCURACY.md').write_text('\n'.join(lines))
    for name in ['LINEAR_ACCURACY.md', 'linear-accuracy.json', 'linear-accuracy-raw.jsonl.gz']:
        shutil.copy2(OUTPUT / name, DEST / name)
    print('COMPLETE: full linear accuracy audited and saved', flush=True)
    print(json.dumps({e: groups['common_linear'][e] for e in ENGINES}), flush=True)


if __name__ == '__main__':
    main()
