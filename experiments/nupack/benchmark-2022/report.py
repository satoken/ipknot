"""Audit raw runs and compare both NUPACK engines with LPC and LPV.

Separate common sets prevent resource-limited exact runs from making the
linear comparison depend on the exact engine's successful sequence lengths.
"""
import argparse
import collections
import csv
import gzip
import json
import math
from pathlib import Path
import random
import statistics

from dataset import load_index, metrics, crossing_pairs, write_json

ENGINES = ['nupack', 'lnupack', 'lpc', 'lpv']
LABELS = dict(nupack='NUPACK (r=0)', lnupack='LinearNUPACK (b=100, r=0)',
              lpc='LPC (b=100, r=1)', lpv='LPV (b=100, r=1)')


def failure(row):
    if row['exit_code'] == 0:
        return 'success'
    if row['timeout']:
        return 'timeout'
    if 'bad_alloc' in row['stderr'] or 'length_error' in row['stderr']:
        return 'memory_limit'
    return 'other_error'


def accuracy(rows, key='accuracy'):
    if not rows:
        return None
    scores = [r[key] for r in rows]
    counts = {k: sum(r[k] for r in scores) for k in ['tp', 'fp', 'fn']}
    tp, fp, fn = (counts[k] for k in ['tp', 'fp', 'fn'])
    return dict(**counts, macro_f1=statistics.mean(r['f1'] for r in scores),
                macro_ppv=statistics.mean(r['ppv'] for r in scores),
                macro_sen=statistics.mean(r['sen'] for r in scores),
                micro_f1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 1)


def summarize(rows, refs):
    good = [r for r in rows if r['exit_code'] == 0]
    pk = [r for r in good if refs[r['id']]['pk']]
    return dict(attempts=len(rows), successes=len(good),
                outcomes=dict(collections.Counter(failure(r) for r in rows)),
                accuracy=accuracy(good), crossing_accuracy_pk_only=accuracy(pk, 'crossing_accuracy'),
                pk_sequences=len(pk),
                macro_f1_failures_zero=sum(r['accuracy']['f1'] for r in good) / len(rows) if rows else None,
                wall_sum_all_seconds=sum(r['wall_seconds'] for r in rows),
                cpu_sum_all_seconds=sum(r['user_seconds'] + r['system_seconds'] for r in rows),
                wall_sum_success_seconds=sum(r['wall_seconds'] for r in good),
                wall_median_success_seconds=statistics.median(r['wall_seconds'] for r in good) if good else None,
                wall_max_success_seconds=max((r['wall_seconds'] for r in good), default=None),
                rss_median_success_mib=statistics.median(r['rss_kib'] for r in good) / 1024 if good else None,
                rss_max_success_mib=max((r['rss_kib'] for r in good), default=0) / 1024)


def paired(a, b):
    ids = sorted(a.keys() & b.keys())
    if not ids:
        return None
    delta = [b[i]['accuracy']['f1'] - a[i]['accuracy']['f1'] for i in ids]
    rng = random.Random(20261010)
    # A paired bootstrap reflects sequence-level variability, not run timing.
    # Avoid quadratic overhead on a full 11,843-entry comparison.
    interval = []
    for _ in range(1000):
        interval.append(sum(delta[rng.randrange(len(delta))] for _ in delta) / len(delta))
    interval.sort()
    return dict(n=len(ids), macro_f1_delta_b_minus_a=statistics.mean(delta),
                bootstrap_95_percentile=[interval[24], interval[974]],
                better=sum(d > 1e-12 for d in delta), equal=sum(abs(d) <= 1e-12 for d in delta),
                worse=sum(d < -1e-12 for d in delta),
                wall_sum_ratio_b_over_a=sum(b[i]['wall_seconds'] for i in ids) / sum(a[i]['wall_seconds'] for i in ids),
                wall_ratio_geomean_b_over_a=math.exp(statistics.mean(math.log(b[i]['wall_seconds'] / a[i]['wall_seconds']) for i in ids)),
                rss_ratio_median_b_over_a=statistics.median(b[i]['rss_kib'] / a[i]['rss_kib'] for i in ids))


def table(rows, refs, ids):
    return {engine: summarize([r for r in rows if r['mode'] == engine and r['id'] in ids], refs)
            for engine in ENGINES}


def markdown_table(t, engines=ENGINES):
    lines = ['| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for engine in engines:
        r = t[engine]
        if not r['successes']:
            continue
        f1 = r['accuracy']['macro_f1']
        pk = r['crossing_accuracy_pk_only']
        pk_text = f"{pk['macro_f1']:.6f}" if pk else '—'
        lines.append(f"| {LABELS[engine]} | {r['successes']} | {f1:.6f} | {pk_text} | "
                     f"{r['wall_median_success_seconds']:.4f} | {r['wall_sum_success_seconds']:.1f} | "
                     f"{r['rss_median_success_mib']:.1f} / {r['rss_max_success_mib']:.1f} |")
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--phase', default='final')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--allow-incomplete', action='store_true')
    ap.add_argument('--no-bootstrap', action='store_true')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    refs = {r['id']: r for r in load_index(args.root)}
    manifest = json.loads((args.root / (args.phase + '-manifest.json')).read_text())
    expected = {(i, e) for i in manifest['sample_ids'] for e in manifest['engines']}
    rows = []
    seen = set()
    for p in sorted((args.root / args.phase).glob('*.json')):
        r = json.loads(p.read_text())
        key = (r['id'], r['mode'])
        if key in seen or key not in expected:
            raise ValueError('unexpected/duplicate run: ' + str(key))
        seen.add(key)
        ref = refs[r['id']]
        if r['input_sha256'] != ref['input_sha256']:
            raise ValueError('input hash mismatch')
        if r['exit_code'] == 0:
            pairs = list(map(tuple, r['pairs']))
            endpoints = [i for pair in pairs for i in pair]
            if len(pairs) != len(set(pairs)) or len(endpoints) != len(set(endpoints)) or any(not 1 <= i < j <= ref['length'] for i, j in pairs):
                raise ValueError('invalid prediction')
            if any(ref['sequence'][i-1] + ref['sequence'][j-1] not in ['AU', 'UA', 'CG', 'GC', 'GU', 'UG'] for i, j in pairs):
                raise ValueError('noncanonical or ambiguous predicted pair: ' + str(key))
            if metrics(ref['pairs'], pairs) != r['accuracy'] or metrics(crossing_pairs(ref['pairs']), crossing_pairs(pairs)) != r['crossing_accuracy']:
                raise ValueError('incorrect metrics')
            pair_set = set(pairs)
            if any((i-1, j+1) not in pair_set and (i+1, j-1) not in pair_set for i, j in pairs):
                raise ValueError('DD prediction has an isolated pair')
        rows.append(r)
    missing = expected - seen
    if missing and not args.allow_incomplete:
        raise ValueError(f'{len(missing)} missing measurements')
    by_engine = {e: {r['id']: r for r in rows if r['mode'] == e and r['exit_code'] == 0} for e in ENGINES}
    common3 = set.intersection(*(set(by_engine[e]) for e in ['lnupack', 'lpc', 'lpv']))
    common4 = set.intersection(*(set(by_engine[e]) for e in ENGINES))
    groups = {}
    for name, ids in [('all_attempts', set(manifest['sample_ids'])), ('common_linear', common3), ('common_all_four', common4)]:
        groups[name] = table(rows, refs, ids)
        for size in ['short', 'medium', 'long']:
            groups[name + '/' + size] = table(rows, refs, {i for i in ids if refs[i]['size'] == size})
        for dataset in ['Rfam14.5', 'bpRNA-1m']:
            groups[name + '/' + dataset] = table(rows, refs, {i for i in ids if refs[i]['dataset'] == dataset})
    comparisons = {}
    if not args.no_bootstrap:
        for a, b in [('lpc', 'lnupack'), ('lpv', 'lnupack'), ('nupack', 'lnupack'), ('lpc', 'nupack'), ('lpv', 'nupack')]:
            comparisons[a + '_vs_' + b] = paired(by_engine[a], by_engine[b])
    result = dict(complete=not missing, expected_runs=len(expected), actual_runs=len(rows),
                  missing_runs=len(missing), common_linear_n=len(common3), common_all_four_n=len(common4),
                  common_linear_length_range=[min((refs[i]['length'] for i in common3), default=None), max((refs[i]['length'] for i in common3), default=None)],
                  common_all_four_length_range=[min((refs[i]['length'] for i in common4), default=None), max((refs[i]['length'] for i in common4), default=None)],
                  groups=groups, comparisons=comparisons,
                  failures=[dict(id=r['id'], member=refs[r['id']]['member'], engine=r['mode'],
                                 category=failure(r), exit_code=r['exit_code'], stderr=r['stderr'])
                            for r in rows if r['exit_code'] != 0],
                  metric_convention='Exact pair equality; per-sequence macro F1; empty/empty F1=1; undefined PPV/SEN=0. PK F1: crossing pairs on reference-PK sequences. Failed runs are excluded from common sets; explicit all-attempt failure-zero macro also recorded.',
                  timing_convention='One fresh process per input/condition, single thread; includes posterior, DD and specified refinement. Wall sums are sums of child times, not elapsed benchmark duration. Shared-host parallel measurements; the auxiliary accuracy-only phase uses disjoint A725 cores during part of the run, and its timings are not merged. Censored exact failures are not successful latencies.')
    write_json(args.output / 'summary.json', result)
    with (args.output / 'per-rna.csv').open('w', newline='') as f:
        fields = ['id', 'member', 'dataset', 'length', 'size', 'pk', 'engine', 'status',
                  'wall_seconds', 'cpu_seconds', 'rss_mib', 'f1', 'ppv', 'sensitivity',
                  'crossing_f1', 'common_linear', 'common_all_four']
        writer = csv.DictWriter(f, fields, lineterminator='\n')
        writer.writeheader()
        for r in rows:
            ref = refs[r['id']]
            row = {k: ref[k] for k in ['id', 'member', 'dataset', 'length', 'size', 'pk']}
            row.update(engine=r['mode'], status=failure(r), wall_seconds=r['wall_seconds'],
                       cpu_seconds=r['user_seconds'] + r['system_seconds'], rss_mib=r['rss_kib'] / 1024,
                       common_linear=r['id'] in common3, common_all_four=r['id'] in common4)
            if r['exit_code'] == 0:
                row.update(f1=r['accuracy']['f1'], ppv=r['accuracy']['ppv'],
                           sensitivity=r['accuracy']['sen'], crossing_f1=r['crossing_accuracy']['f1'])
            writer.writerow(row)
    with gzip.open(args.output / 'raw-results.jsonl.gz', 'wt') as f:
        for r in rows:
            f.write(json.dumps(r, separators=(',', ':')) + '\n')
    lines = ['# NUPACK・LinearNUPACK・LPC・LPV ベンチマーク', '',
             f"状態: {'完了' if result['complete'] else '計測中'}。{len(rows):,}/{len(expected):,} runs。", '',
             '公式2022年IPknot単一配列セット全11,843配列（bpRNA-1m 3,437 / Rfam 14.5 8,406、12–4,381 nt）。',
             'NUPACK・LinearNUPACKはrefinement 0、LPC・LPVはrefinement 1。全方式でDD beam 100 / 50 iterations、crossing beam 100 / witnesses 16、自動閾値。線形エンジンのbeamは100。',
             '通常NUPACKは1配列60秒、他方式は1,800秒。全方式の仮想アドレス空間上限は8 GiB。曖昧塩基は位置を保持して対を作らない。', '',
             '## 完了率と全件の状況', '',
             '| 方式 | 成功 | 時間切れ | メモリ上限 | その他エラー |', '|---|---:|---:|---:|---:|']
    for e in ENGINES:
        outcomes = groups['all_attempts'][e]['outcomes']
        lines.append(f"| {LABELS[e]} | {outcomes.get('success',0)} | {outcomes.get('timeout',0)} | {outcomes.get('memory_limit',0)} | {outcomes.get('other_error',0)} |")
    if (args.output / 'LINEAR_ACCURACY.md').exists():
        lines += ['', '線形3方式の全11,843配列の精度は [LINEAR_ACCURACY.md](LINEAR_ACCURACY.md) に先行集計済みです。以下の実行時間は元のX925 phaseのみです。', '']
    lines += ['', f"## 線形3方式の共通成功集合（{len(common3):,}配列）", '',
              markdown_table(groups['common_linear'], ['lnupack', 'lpc', 'lpv']), '',
              f"## 4方式の共通成功集合（{len(common4):,}配列）", '', markdown_table(groups['common_all_four']), '',
              f"4方式の共通集合の長さ範囲: {result['common_all_four_length_range']} nt。通常版が完了した短い配列に偏るため、この集合を線形版の全体精度と混同しない。各集合内では全方式が同じ配列で比較される。", '']
    for size, label in [('short', '短鎖（150 nt以下）'), ('medium', '中鎖（151–500 nt）'), ('long', '長鎖（500 nt超）')]:
        lines += [f'## 線形共通集合: {label}', '', markdown_table(groups['common_linear/' + size], ['lnupack', 'lpc', 'lpv']), '']
    for dataset in ['Rfam14.5', 'bpRNA-1m']:
        lines += [f'## 線形共通集合: {dataset}', '', markdown_table(groups['common_linear/' + dataset], ['lnupack', 'lpc', 'lpv']), '']
    lines += ['## 指標・測定条件', '', result['metric_convention'], '', result['timing_convention'], '',
              'CPU: Cortex-X925の10コアに固定（5–9,15–19）。追加の精度評価は別のA725コアを使用し、その測定時間はこの表に含めない。同じ配列の各方式は同じCPUで実行し、方式順を入力ハッシュに基づきシャッフル。GCC 13.3、Release -O3 -DNDEBUG、ENABLE_ILP=OFF。', '',
              '精度で配列を選択せず、パラメータ調整やBPPキャッシュも行わない。通常版の打ち切りはその構造が誤りという意味ではなく、指定資源内で結果が得られなかったことを示す。', '',
              '通常NUPACKには既存のパラメータ読み込み・疑似結び目周辺確率の修正が入る。LinearNUPACKは制限付き文法とbeam pruningを使う近似であり、幅100でも通常版と一致する保証はない。', '',
              '生記録: `raw-results.jsonl.gz`。配列別結果: `per-rna.csv`。集計・paired bootstrap: `summary.json`。条件・ハッシュ: `final-manifest.json`, `provenance.json`。再実行方法: `README.md`。', '', '## グラフ', '']
    if result['complete']:
        lines += ['![Accuracy by length](figures/accuracy.png)', '', '![Runtime and memory by length](figures/scaling.png)', '']
    else:
        lines += ['グラフは全件の測定が完了してから生成します。', '']
    (args.output / 'RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps({k: result[k] for k in ['complete', 'actual_runs', 'missing_runs', 'common_linear_n', 'common_all_four_n']}))
    for e, r in groups['all_attempts'].items():
        print(e, r['outcomes'], r['accuracy'])


if __name__ == '__main__':
    main()
