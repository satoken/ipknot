#!/usr/bin/env python3
"""Audit complete measurements and export compact, reproducible dataset results."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics as stats

from benchmark_dataset import load_helper, options


def audit(root, phase, dataset, refs):
    manifest = json.loads((root / (phase + '-manifest.json')).read_text())
    expected = {(i, m, rep) for i in manifest['sample_ids'] for m in manifest['modes']
                for rep in range(manifest['repetitions'])}
    if dataset.digest(root/'references.json.gz') != manifest['hashes']['references']:
        raise ValueError('changed reference index')
    for name in ('baseline', 'candidate', 'launcher', 'helper'):
        if dataset.digest(manifest[name]) != manifest['hashes'][name]:
            raise ValueError('changed ' + name)
    if manifest['options'] != {m: options(m) for m in manifest['modes']}:
        raise ValueError('changed options')
    rows = {}
    for path in sorted((root / phase).glob('*.json')):
        row = json.loads(path.read_text())
        key = row['id'], row['mode'], row['repetition']
        if key not in expected or key in rows:
            raise ValueError('unexpected/duplicate result: ' + str(key))
        ref = refs[row['id']]
        binary = manifest['baseline'] if row['mode'].endswith('before') else manifest['candidate']
        command = [binary, *options(row['mode']), '-B', str(path.with_suffix('.bpseq')),
                   str(root/'inputs'/(row['id']+'.fa'))]
        threads = 4 if row['mode'].endswith(('wave4', 'inner4')) else 2 if row['mode'].endswith('level2') else 1
        if (row['phase'] != phase or row['command'] != command or
                row['input_sha256'] != ref['input_sha256'] or
                row['configuration_sha256'] != manifest['configuration_sha256'] or
                row['cpus'] not in [sorted(g[:threads]) for g in manifest['cpu_groups']]):
            raise ValueError('incorrect provenance: ' + str(key))
        for field in ('wall_seconds', 'user_seconds', 'system_seconds', 'rss_kib'):
            if not math.isfinite(row[field]) or row[field] < 0:
                raise ValueError('invalid resource measure: ' + str(key))
        if row['exit_code'] == 0:
            pairs = set(map(tuple, row['pairs']))
            endpoints = [v for p in pairs for v in p]
            if (len(pairs) != len(row['pairs']) or len(set(endpoints)) != len(endpoints) or
                    any(not 1 <= i < j <= ref['length'] for i, j in pairs)):
                raise ValueError('invalid pair coordinates: ' + str(key))
            if any(ref['sequence'][i-1] + ref['sequence'][j-1] not in
                   ('AU', 'UA', 'CG', 'GC', 'GU', 'UG') for i, j in pairs):
                raise ValueError('noncanonical/ambiguous pair: ' + str(key))
            if any((i-1, j+1) not in pairs and (i+1, j-1) not in pairs for i, j in pairs):
                raise ValueError('isolated pair: ' + str(key))
            if (dataset.metrics(ref['pairs'], pairs) != row['accuracy'] or
                    dataset.metrics(dataset.crossing_pairs(ref['pairs']),
                                    dataset.crossing_pairs(pairs)) != row['crossing_accuracy']):
                raise ValueError('incorrect accuracy: ' + str(key))
        rows[key] = row
    if expected != rows.keys():
        raise ValueError(f'{phase}: missing {len(expected-rows.keys())} runs')
    if manifest.get('completed_runs') != len(rows):
        raise ValueError('runner did not finish: ' + phase)
    for identifier in manifest['sample_ids']:
        group = [rows[identifier, m, rep] for m in manifest['modes']
                 for rep in range(manifest['repetitions'])]
        if len({min(r['cpus']) for r in group}) != 1:
            raise ValueError('different CPU groups for paired RNA: ' + identifier)
        for mode in manifest['modes']:
            successes = [r['pairs'] for r in group if r['mode'] == mode and r['exit_code'] == 0]
            if any(p != successes[0] for p in successes):
                raise ValueError('nondeterministic repetitions: ' + identifier + ' ' + mode)
    return manifest, rows


def accuracy(rows, field='accuracy'):
    if not rows:
        return None
    counts = {k: sum(r[field][k] for r in rows) for k in ('tp', 'fp', 'fn')}
    denom = 2*counts['tp'] + counts['fp'] + counts['fn']
    return dict(**counts, macro_f1=stats.mean(r[field]['f1'] for r in rows),
                macro_ppv=stats.mean(r[field]['ppv'] for r in rows),
                macro_sen=stats.mean(r[field]['sen'] for r in rows),
                micro_f1=2*counts['tp']/denom if denom else 1.)


def aggregate(rows, refs):
    good = [r for r in rows if r['exit_code'] == 0]
    return dict(n=len(rows), failures=len(rows)-len(good), accuracy=accuracy(good),
                crossing_accuracy_pk_only=accuracy([r for r in good if refs[r['id']]['pk']], 'crossing_accuracy'),
                wall_sum_seconds=sum(r['wall_seconds'] for r in rows),
                cpu_sum_seconds=sum(r.get('cpu_seconds', r['user_seconds']+r['system_seconds']) for r in rows),
                wall_median_seconds=stats.median(r['wall_seconds'] for r in rows),
                rss_median_mib=stats.median(r['rss_kib'] for r in rows)/1024,
                rss_max_mib=max(r['rss_kib'] for r in rows)/1024)


def compare(a, b):
    identifiers = sorted(a.keys() & b.keys())
    good = [i for i in identifiers if a[i]['exit_code'] == b[i]['exit_code'] == 0]
    delta = [b[i]['accuracy']['f1']-a[i]['accuracy']['f1'] for i in good]
    return dict(n=len(identifiers), common_success=len(good),
                speedup_wall_sum=sum(a[i]['wall_seconds'] for i in good)/sum(b[i]['wall_seconds'] for i in good),
                speedup_wall_geomean=math.exp(stats.mean(math.log(a[i]['wall_seconds']/b[i]['wall_seconds']) for i in good)),
                cpu_ratio=sum(b[i].get('cpu_seconds', b[i]['user_seconds']+b[i]['system_seconds']) for i in good)/
                          sum(a[i].get('cpu_seconds', a[i]['user_seconds']+a[i]['system_seconds']) for i in good),
                same_structure=sum(a[i]['pairs'] == b[i]['pairs'] for i in good),
                macro_f1_delta=stats.mean(delta), better=sum(d>1e-12 for d in delta),
                equal=sum(abs(d)<=1e-12 for d in delta), worse=sum(d < -1e-12 for d in delta))


def groups(refs, ids):
    result = {'all': sorted(ids)}
    for size in ('short', 'medium', 'long'):
        result[size] = sorted(i for i in ids if refs[i]['size'] == size)
    for source in ('Rfam14.5', 'bpRNA-1m'):
        result[source] = sorted(i for i in ids if refs[i]['dataset'] == source)
    for pk in (False, True):
        result['PK' if pk else 'PKF'] = sorted(i for i in ids if refs[i]['pk'] == pk)
    for source in ('Rfam14.5', 'bpRNA-1m'):
        for size in ('short', 'medium', 'long'):
            for pk in (False, True):
                result[f'{source}/{size}/{"PK" if pk else "PKF"}'] = sorted(
                    i for i in ids if refs[i]['dataset'] == source and refs[i]['size'] == size and refs[i]['pk'] == pk)
    return result


def repeated_summary(manifest, rows, refs):
    medians = {}
    for mode in manifest['modes']:
        medians[mode] = {}
        for i in manifest['sample_ids']:
            reps = [rows[i, mode, rep] for rep in range(manifest['repetitions'])]
            row = dict(reps[0])
            for field in ('wall_seconds', 'user_seconds', 'system_seconds', 'rss_kib'):
                row[field] = stats.median(r[field] for r in reps)
            row['cpu_seconds'] = stats.median(r['user_seconds']+r['system_seconds'] for r in reps)
            medians[mode][i] = row
    result = {}
    for name, ids in groups(refs, set(manifest['sample_ids'])).items():
        by = {m: {i: medians[m][i] for i in ids} for m in medians}
        comparisons = {}
        for engine in ('contra', 'lpc', 'improved'):
            before, serial = engine+'-before', engine+'-serial'
            for mode in by:
                if not mode.startswith(engine+'-'):
                    continue
                if before in by and mode != before:
                    comparisons[f'{before} -> {mode}'] = compare(by[before], by[mode])
                if serial in by and mode not in (before, serial):
                    comparisons[f'{serial} -> {mode}'] = compare(by[serial], by[mode])
        if 'contra-wave1' in by and 'contra-wave4' in by:
            comparisons['contra-wave1 -> contra-wave4'] = compare(by['contra-wave1'], by['contra-wave4'])
        result[name] = dict(modes={m: aggregate(list(v.values()), refs) for m, v in by.items()}, comparisons=comparisons)
    return medians, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--historical-contra', type=Path)
    parser.add_argument('--historical-lpc', type=Path)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    dataset = load_helper(args.helper)
    refs = {r['id']: r for r in dataset.load_index(root)}
    full_manifest, full = audit(root, 'full', dataset, refs)
    timing_manifest, timing = audit(root, 'timing', dataset, refs)
    default_manifest, default = audit(root, 'default-policy', dataset, refs)
    for name in ('baseline', 'candidate', 'launcher', 'helper', 'references'):
        if not full_manifest['hashes'][name] == timing_manifest['hashes'][name] == default_manifest['hashes'][name]:
            raise ValueError('full/timing executables or inputs differ')
    if timing_manifest['sample_ids'] != default_manifest['sample_ids']:
        raise ValueError('wait-policy samples differ')
    for (i, mode, rep), r in [*timing.items(), *default.items()]:
        if (i, mode, 0) in full and r['exit_code'] == full[i, mode, 0]['exit_code'] == 0:
            if r['pairs'] != full[i, mode, 0]['pairs']:
                raise ValueError('full/timing prediction differs: ' + i + ' ' + mode)
    by = {m: {i: full[i, m, 0] for i in full_manifest['sample_ids']} for m in full_manifest['modes']}
    full_groups = {name: {m: aggregate([by[m][i] for i in ids], refs) for m in by}
                   for name, ids in groups(refs, set(full_manifest['sample_ids'])).items()}
    medians, timing_groups = repeated_summary(timing_manifest, timing, refs)
    default_medians, default_groups = repeated_summary(default_manifest, default, refs)
    invariant_pairs = [('lpc-serial', 'lpc-level2'), ('lpc-serial', 'lpc-inner4'),
                       ('contra-serial', 'contra-level2'), ('contra-wave1', 'contra-wave4'),
                       ('improved-serial', 'improved-level2'), ('improved-serial', 'improved-inner4')]
    for repeated_modes in (medians, default_medians):
        for a, b in invariant_pairs:
            if a in repeated_modes and b in repeated_modes:
                c = compare(repeated_modes[a], repeated_modes[b])
                if c['same_structure'] != c['common_success']:
                    raise ValueError('thread-dependent output: ' + a + ' ' + b)
    policy_comparisons = {}
    for mode in default_medians:
        c = compare(medians[mode], default_medians[mode])
        if c['same_structure'] != c['common_success']:
            raise ValueError('wait-policy-dependent output: ' + mode)
        c['same_cpu_group'] = sum(medians[mode][i]['cpus'] == default_medians[mode][i]['cpus']
                                  for i in default_manifest['sample_ids'])
        policy_comparisons[mode] = c
    failures = [r for r in [*full.values(), *timing.values(), *default.values()] if r['exit_code'] != 0]
    historical = {}
    for model, source, old_mode, new_mode in [
            ('contra', args.historical_contra, 'ddcontranuss50', 'contra-serial'),
            ('lpc', args.historical_lpc, 'default', 'lpc-serial')]:
        if source is None:
            continue
        old = {}
        with gzip.open(source, 'rt') as stream:
            for line in stream:
                r = json.loads(line)
                if r['mode'] == old_mode and r['phase'] == 'full' and r['repetition'] == 0:
                    if r['id'] in old or r['input_sha256'] != refs[r['id']]['input_sha256']:
                        raise ValueError('incorrect historical input')
                    if dataset.metrics(refs[r['id']]['pairs'], r['pairs']) != r['accuracy']:
                        raise ValueError('incorrect historical accuracy')
                    old[r['id']] = r
        if old.keys() != by[new_mode].keys():
            raise ValueError('incomplete historical comparison')
        c = compare(old, by[new_mode])
        # Historical binaries predate the designated performance baseline. Only
        # use their predictions as accuracy references, never as matched timings.
        for k in ('speedup_wall_sum', 'speedup_wall_geomean', 'cpu_ratio'):
            c.pop(k)
        historical[model] = dict(comparison=c, source_sha256=dataset.digest(source),
                                 source=str(source.resolve()), old_accuracy=accuracy(list(old.values())),
                                 new_accuracy=accuracy(list(by[new_mode].values())),
                                 scope='Historical full-dataset accuracy reference; not the 454db3b matched performance baseline.')
    column_diagnostics = json.loads((root/'column-diagnostics.json').read_text())
    if column_diagnostics['binary_sha256'] != full_manifest['hashes']['candidate']:
        raise ValueError('diagnostic binary differs')
    expected_diagnostics = {(i, m) for i in timing_manifest['sample_ids']
                            for m in ('lpc-serial', 'improved-serial')}
    diagnostic_keys = [(r['id'], r['mode']) for r in column_diagnostics['rows']]
    if (len(diagnostic_keys) != len(expected_diagnostics) or
            set(diagnostic_keys) != expected_diagnostics or column_diagnostics['threshold'] != 8192):
        raise ValueError('incomplete column diagnostic')
    for r in column_diagnostics['rows']:
        ref = refs[r['id']]
        upper = 100 * max(r['graph_pair_counts']) / ref['length']
        if (r['input_sha256'] != ref['input_sha256'] or r['length'] != ref['length'] or
                r['average_offers_upper_bound'] != upper or r['proves_serial_fallback'] != (upper < 8192)):
            raise ValueError('incorrect column diagnostic bound')
        f = full[r['id'], r['mode'], 0]
        if hashlib.sha256(json.dumps(f['pairs'],separators=(',', ':')).encode()).hexdigest() != r['prediction_sha256']:
            raise ValueError('diagnostic prediction differs')
    summary = dict(full=full_groups, timing=timing_groups, default_wait_policy=default_groups,
                   wait_policy_comparisons=policy_comparisons, historical_accuracy=historical,
                   column_diagnostics=dict(n=len(column_diagnostics['rows']),
                                           guaranteed_serial_fallback=sum(r['proves_serial_fallback'] for r in column_diagnostics['rows']),
                                           maximum_average_offers_upper_bound=max(r['average_offers_upper_bound'] for r in column_diagnostics['rows'])),
                   audit=dict(full_runs=len(full), timing_runs=len(timing), default_wait_policy_runs=len(default), failures=len(failures),
                              nondeterministic_repetitions=0, invalid_structures=0,
                              full_timing_serial_mismatches=0, thread_dependent_structures=0),
                   notes=['Full dataset timing sums are aggregate child times, not study elapsed time.',
                          'Timing uses medians of three repetitions per RNA; all probability computations and refinement are included.',
                          'Short/medium/long timing strata have equal sample counts; the unweighted timing sum is not the full-dataset speedup.',
                          'Exact backends may choose different equal-score structures; the global DD procedure remains approximate.'])
    out.mkdir(parents=True, exist_ok=True)
    dataset.write_json(out/'summary.json', summary)
    with gzip.open(out/'raw-results.jsonl.gz', 'wt') as stream:
        for row in [*full.values(), *timing.values(), *default.values()]:
            stream.write(json.dumps(row, separators=(',', ':'))+'\n')
    fields = ['phase', 'id', 'mode', 'repetition', 'dataset', 'length', 'size', 'pk',
              'exit_code', 'wall_seconds', 'cpu_seconds', 'rss_mib', 'macro_f1', 'prediction_sha256']
    with (out/'per-rna.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for row in [*full.values(), *timing.values(), *default.values()]:
            ref = refs[row['id']]
            writer.writerow({**{k: row[k] for k in ('phase', 'id', 'mode', 'repetition', 'exit_code', 'wall_seconds')},
                             **{k: ref[k] for k in ('dataset', 'length', 'size', 'pk')},
                             'cpu_seconds': row['user_seconds']+row['system_seconds'],
                             'rss_mib': row['rss_kib']/1024,
                             'macro_f1': row.get('accuracy', {}).get('f1'),
                             'prediction_sha256': hashlib.sha256(json.dumps(row.get('pairs'),separators=(',', ':')).encode()).hexdigest()})
    for name in ('references.json.gz', 'inventory.json', 'provenance.json', 'dataset-audit.json', 'full-manifest.json', 'timing-manifest.json', 'default-policy-manifest.json', 'column-diagnostics.json',
                 'baseline-flags.make', 'candidate-flags.make', 'baseline-link.txt', 'candidate-link.txt'):
        if (root/name).exists():
            if name.endswith('-flags.make'):
                (out/name).write_text((root/name).read_text().rstrip()+'\n')
            else:
                shutil.copy2(root/name, out/name)
    for path in (root/'bin').iterdir():
        with path.open('rb') as src, gzip.open(out/(path.name+'.gz'), 'wb') as dest:
            shutil.copyfileobj(src, dest)
    for phase, source, manifest in [('full', root/'full-runner.py', full_manifest),
                                    ('timing', root/'timing-runner.py', timing_manifest),
                                    ('default-policy', Path(__file__).with_name('benchmark_dataset.py'), default_manifest)]:
        if dataset.digest(source) != manifest['hashes']['script']:
            raise ValueError('runner source hash differs: ' + phase)
        with source.open('rb') as src, gzip.open(out/(phase+'-runner.py.gz'), 'wb') as dest:
            shutil.copyfileobj(src, dest)
    t = timing_groups['all']
    d = default_groups['all']
    contra_gain = t['comparisons']['contra-before -> contra-serial']
    lpc_parallel = d['comparisons']['improved-serial -> improved-level2']
    standard = full_groups['all']['lpc-serial']
    improved = full_groups['all']['improved-serial']
    lpc_algorithm_speedup = standard['wall_sum_seconds'] / improved['wall_sum_seconds']
    lpc_parallel_algorithm_speedup = (d['modes']['lpc-level2']['wall_sum_seconds'] /
                                      d['modes']['improved-level2']['wall_sum_seconds'])
    lines = [
        '# 2022公開データでのNussinov高速化検証', '',
        '改修前 `454db3b` と高速化後 `eb7cc1c` を比較。全11,843配列を高速化後の厳密DP・標準/改良beamで各1スレッド実行し、'
        '従来と同じ層別36配列では改修前後・並列方式を各3回測定した。', '',
        f"全件 {len(full):,} 実行、ACTIVE反復 {len(timing):,} 実行、既定待機反復 {len(default):,} 実行。失敗 {len(failures)}、"
        '不正構造・反復間の構造差・対応する逐次版とのスレッド数による構造差はいずれも0。', '',
        f"CONTRAfold＋厳密Nussinovのsparse逐次版は、対応する改修前比で {contra_gain['speedup_wall_sum']:.3f} 倍。"
        f"通常設定でのLPC＋改良LinearNussinovのレベル2並列は、改良版逐次比で {lpc_parallel['speedup_wall_sum']:.3f} 倍、"
        f"CPU時間比は {lpc_parallel['cpu_ratio']:.3f}。速度は36配列の対応比較であり、全件速度比ではない。", '',
        '## 設定と測定範囲', '',
        '- 2022公式ZIPの全単一配列：bpRNA-1m 3,437、Rfam14.5 8,406、12–4,381塩基。'
        '配列を除外せず、曖昧塩基の位置と参照の非標準対を保持。ZIPから全参照・FASTA hashを再照合。',
        '- 自動閾値、refinement 1、DD最大50反復、patience 0、交差幅100、witnesses 16。'
        '孤立対は禁止。追加PK/NMR/NOE・上界・交換・回収は既定の無効設定。',
        '- CONTRAfoldは厳密Nussinov、LPCはposterior beam 100、標準/改良LinearNussinovはbeam 100。'
        'LPC-DDの推奨は改良版`improved-beam`。標準版も比較対象として測定した。',
        '- Release、GCC 13.3、aarch64、ILPなし。全件はX925の10コアに独立した1スレッド予測を分散。'
        '時間総和は子プロセスの合計であり、全件評価の経過時間とは異なる。',
        '- 反復測定はX925の独立した4コア群2組。各RNAの全条件・反復は同じ群、'
        '逐次版は先頭1コア、レベル並列は先頭2コア、DP内並列は4コアに固定。'
        '条件順は入力hashと反復番号で入れ替えた。動作周波数は固定していない。'
        'ACTIVEと既定待機のパスは分けて実施した。',
        '- 全計測はFASTAから新しいプロセスを起動し、BPP計算・閾値探索・refinementを含む。'
        'BPPキャッシュ・traceなし。1実行1,800秒、仮想空間8 GiB。CLOCK_MONOTONIC実時間、wait4 CPU/RSS。', '',
        '## 全件の精度・資源', '',
        '| 方式（1スレッド） | n | macro F1 | micro F1 | CPU合計 (s) | 実時間合計 (s) | 最大RSS (MiB) |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for mode, s in full_groups['all'].items():
        lines.append(f"| {mode} | {s['n']} | {s['accuracy']['macro_f1']:.6f} | {s['accuracy']['micro_f1']:.6f} | "
                     f"{s['cpu_sum_seconds']:.3f} | {s['wall_sum_seconds']:.3f} | {s['rss_max_mib']:.2f} |")
    lines += ['', '| 対象 | n | CONTRAfold Nussinov F1 | LPC beam F1 | LPC improved F1 |',
              '| --- | ---: | ---: | ---: | ---: |']
    for group in ('short', 'medium', 'long', 'Rfam14.5', 'bpRNA-1m', 'PKF', 'PK'):
        s = full_groups[group]
        lines.append(f"| {group} | {s['contra-serial']['n']} | {s['contra-serial']['accuracy']['macro_f1']:.6f} | "
                     f"{s['lpc-serial']['accuracy']['macro_f1']:.6f} | {s['improved-serial']['accuracy']['macro_f1']:.6f} |")
    lines += ['', 'short ≤150、medium 151–500、long >500。PK分類は参照構造による。'
              'F1は位置が完全一致する塩基対で計算し、macroは配列ごとの平均、microはTP/FP/FNを合算。'
              '空参照・空予測のF1は1、分母0のPPV/SENは0。', '',
              '## 通常の待機方式での並列効果', '',
              '`OMP_WAIT_POLICY`と`GOMP_SPINCOUNT`を未設定にし、通常のGNU OpenMP既定値で新版を再測定した。'
              '各RNAの全条件・3反復は同じCPU群。数値は配列ごとの3回の中央値から集計した。'
              '推奨判断には、このパス内での逐次・並列比較を用いる。', '',
              '| 方式 | 時間中央値の総和 (s) | 新版逐次比の速度 (倍) | CPU比（並列/逐次） | 最大RSS (MiB) |',
              '| --- | ---: | ---: | ---: | ---: |']
    for mode, s in d['modes'].items():
        serial = mode.split('-')[0]+'-serial'
        c = d['comparisons'].get(serial+' -> '+mode)
        lines.append(f"| {mode} | {s['wall_sum_seconds']:.3f} | " +
                     (f"{c['speedup_wall_sum']:.3f}" if c else '1.000') + ' | ' +
                     (f"{c['cpu_ratio']:.3f}" if c else '1.000') + f" | {s['rss_max_mib']:.2f} |")
    lines += ['', '`wave4`と`sparse serial`の比較は厳密DP方式の違いも含む。'
              '同じwavefront方式内での1/4スレッド対照は、下のACTIVEパスに含めた。'
              'どの測定もBPP計算とrefinementを含む。', '',
              '| 長さ | n | CONTRA serial (s) | CONTRA level2 (s) | wave4 (s) | LPC serial (s) | LPC level2 (s) | improved serial (s) | improved level2 (s) |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for group in ('short', 'medium', 'long'):
        s = default_groups[group]['modes']
        lines.append(f"| {group} | {s['contra-serial']['n']} | " + ' | '.join(
            f"{s[m]['wall_sum_seconds']:.3f}" for m in ('contra-serial', 'contra-level2', 'contra-wave4', 'lpc-serial', 'lpc-level2', 'improved-serial', 'improved-level2')) + ' |')
    lines += ['', '## 改修前との対応反復測定（ACTIVE）', '',
              '最初の比較はDP内の短いバリアで待機方式を揃えるため、`OMP_WAIT_POLICY=ACTIVE`で行った。'
              'GNU libgompでは未設定時とACTIVE時で待機のspin回数が異なるため、'
              'CPU時間の解釈と通常設定の推奨用に上の既定待機パスを追加した。'
              '仕様は[GNU libgomp GOMP_SPINCOUNT](https://gcc.gnu.org/onlinedocs/libgomp/GOMP_005fSPINCOUNT.html)を参照。', '',
              '各配列3回の実時間中央値を合計した。RSSも各配列3回の中央値を求め、その最大値を示す。'
              '36配列はデータ源×長さ×参照PK有無の12層から各3配列。'
              '各層のサンプル数を揃えているため、この合計の速度比を全11,843配列の速度比と同一視できない。', '',
              '| 方式 | 時間中央値の総和 (s) | 改修前比の速度 (倍) | 同じ新版逐次方式比 (倍) | CPU比（新版/改修前） | 最大RSS (MiB) |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    for mode, s in t['modes'].items():
        engine = mode.split('-')[0]
        key = engine+'-before -> '+mode
        before = t['comparisons'].get(key)
        serial = t['comparisons'].get(engine+'-serial -> '+mode)
        lines.append(f"| {mode} | {s['wall_sum_seconds']:.3f} | " +
                     (f"{before['speedup_wall_sum']:.3f}" if before else '1.000') + ' | ' +
                     (f"{serial['speedup_wall_sum']:.3f}" if serial else '—') + ' | ' +
                     (f"{before['cpu_ratio']:.3f}" if before else '1.000') + f" | {s['rss_max_mib']:.2f} |")
    lines += ['', '`before`＝改修前1スレッド、`serial`＝新版1スレッド、`level2`＝独立レベル2スレッド、'
              '`wave1/4`＝厳密wavefrontの1/4スレッド、`inner4`＝列内4スレッド要求。'
              '同じ新版逐次方式比はNussinovのsparse→wavefrontの方式差を含む。'
              'wavefront自体の並列効果は次の比で切り分けられる。', '',
              f"wavefront 1→4スレッドの実時間総和比は {t['comparisons']['contra-wave1 -> contra-wave4']['speedup_wall_sum']:.3f} 倍。", '',
              '| 長さ | n | CONTRA 改修前 (s) | sparse新版 (s) | level2 (s) | wave4 (s) | LPC改修前 (s) | LPC新版 (s) | LPC level2 (s) |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for group in ('short', 'medium', 'long'):
        s = timing_groups[group]['modes']
        lines.append(f"| {group} | {s['contra-before']['n']} | " + ' | '.join(
            f"{s[m]['wall_sum_seconds']:.3f}" for m in ('contra-before', 'contra-serial', 'contra-level2', 'contra-wave4', 'lpc-before', 'lpc-serial', 'lpc-level2')) + ' |')
    lines += ['', '## 待機方式の比較', '',
              'ACTIVEと既定待機は別時刻のパスであり、方式間の比は参考値。'
              '同じCPU群を使ったRNA数も示す。通常の並列効果は、上の既定待機パス内部の対応比較で判断する。', '',
              '| 方式 | 既定/ACTIVEの速度 (倍) | CPU比（既定/ACTIVE） | 同じCPU群 | 同一構造 |',
              '| --- | ---: | ---: | ---: | ---: |']
    for mode, c in policy_comparisons.items():
        lines.append(f"| {mode} | {c['speedup_wall_sum']:.3f} | {c['cpu_ratio']:.3f} | {c['same_cpu_group']}/{c['n']} | {c['same_structure']}/{c['n']} |")
    lines += ['', '## このデータでの推奨', '',
              '厳密Nussinovはsparse逐次版を基本設定とする。通常待機でのレベル2並列は新版逐次比 '
              f"{d['comparisons']['contra-serial -> contra-level2']['speedup_wall_sum']:.3f} 倍、"
              f"wavefront 4スレッドは {d['comparisons']['contra-serial -> contra-wave4']['speedup_wall_sum']:.3f} 倍であり、"
              f"後者のCPU時間は {d['comparisons']['contra-serial -> contra-wave4']['cpu_ratio']:.3f} 倍だった。"
              '単一DP内のwavefront並列化は実装済みだが、この推奨パラメータでの全体処理には追加の利点が小さい。', '',
              'LPC-DDは`--dd-dp improved-beam --dd-beam 100`を推奨する。'
              f"全11,843配列で標準版→改良版のmacro F1は {standard['accuracy']['macro_f1']:.6f}→{improved['accuracy']['macro_f1']:.6f}、"
              f"micro F1は {standard['accuracy']['micro_f1']:.6f}→{improved['accuracy']['micro_f1']:.6f}、"
              f"1スレッドでの実時間総和比は {lpc_algorithm_speedup:.3f} 倍だった。"
              '中長配列のmacro F1と交差対のみのF1は微減しており、精度向上は全体評価での結果である。', '',
              '単一配列の待ち時間を短縮する場合は、改良beam 100のレベル2並列を推奨する。'
              f"通常設定では改良版逐次比 {lpc_parallel['speedup_wall_sum']:.3f} 倍、CPU時間比 {lpc_parallel['cpu_ratio']:.3f}。"
              f"両方式をレベル2並列にした場合も改良版は標準版比 {lpc_parallel_algorithm_speedup:.3f} 倍高速だった。"
              '多数のRNAを処理する場合は、各プロセス1スレッドで独立したRNAをコアへ割り当てる運用を基本とする。'
              '列内並列は推奨beam幅で逐次フォールバックとなるため、その要求スレッド数を増やす効果は確認できなかった。', '',
              '対応するCLI設定（並列版はOpenMP対応ビルドが必要）：', '', '```sh',
              'ipknot '+' '.join(options('contra-serial'))+' input.fa',
              'ipknot '+' '.join(options('improved-level2'))+' input.fa',
              '```', '',
              '`--dd-dp-threads`は1とし、LPCのレベル並列には`--dd-threads 2`を使う。'
              '測定結果の`lpc-*`は標準beamの比較対象、`improved-*`は推奨する改良beamを表す。']
    lines += ['', '## 正しさと列内並列の診断', '',
              '全件で対の座標、端点の排他性、正準対、孤立対禁止、参照から再計算した精度を検証。'
              '反復間と全件/反復の逐次版が一致。レベル並列と対応する逐次版、beam列内要求と逐次版、'
              'wavefrontの1/4スレッドはそれぞれ同じ構造だった。'
              '厳密DPの方式間や改修前後では等得点の選択とDDの軌跡が変わり得る。'
              'oracleの設計と独立DP/全探索テストは[性能検証資料](../../INTRA_DP.md)を参照。', '',
              '| 対応比較（36配列） | 同一構造 | F1差 macro | 改善 / 同点 / 悪化 |',
              '| --- | ---: | ---: | ---: |']
    for key, c in t['comparisons'].items():
        if key.split(' -> ')[0].endswith('before'):
            lines.append(f"| {key} | {c['same_structure']}/{c['n']} | {c['macro_f1_delta']:+.6f} | {c['better']} / {c['equal']} / {c['worse']} |")
    diag = summary['column_diagnostics']
    lines += ['', f"標準/改良beamの36配列、計{diag['n']}ケースはすべて列内の逐次フォールバック条件を満たした。"
              f"全レベルを合わせた候補数から求めた平均prefix offer数の上限でも最大 {diag['maximum_average_offers_upper_bound']:.3f}、"
              '切替基準は8192。回収のmaskは候補を減らすので、この上限は回収にも適用できる。'
              'この推奨設定で`inner4`を指定した結果を、実働4スレッドの並列速度として解釈できない。'
              'infoログを使う診断はA725の別CPUで実施し、時間集計には含めていない。', '',
              '過去の全件出力は精度のみの参照として照合した。過去の実行ファイルは指定の改修前revisionと異なるため、'
              '速度の比較には使用していない。']
    lines.append('')
    for model, h in historical.items():
        c = h['comparison']
        lines.append(f"- {model}: 同一構造 {c['same_structure']}/{c['n']}、macro F1差 {c['macro_f1_delta']:+.6f}。")
    lines += ['', '## 保存物', '',
              '- `summary.json`：全体・長さ・データ源・PK・全12層、対応速度/精度。',
              '- `raw-results.jsonl.gz`：全予測、コマンド、終了状態、資源、入力hash。`per-rna.csv`：各実行の表。',
              '- `references.json.gz`、`inventory.json`、`dataset-audit.json`、各manifest、`provenance.json`、`SHA256SUMS`：入力/環境/ハッシュ。',
              '- `baseline.gz` / `candidate.gz` / `measure.gz`：凍結した実行ファイル。各`*-runner.py.gz`：使用した測定コード。',
              '- `column-diagnostics.json`：列内仕事量の上限と診断予測hash。診断の予測は全件パスと一致。', '',
              '全件開始/終了 UTC: '+full_manifest['started_utc']+' / '+full_manifest['completed_utc']+'。',
              'ACTIVE反復開始/終了 UTC: '+timing_manifest['started_utc']+' / '+timing_manifest['completed_utc']+'。',
              '既定待機反復開始/終了 UTC: '+default_manifest['started_utc']+' / '+default_manifest['completed_utc']+'。', '']
    (out/'RESULTS.md').write_text('\n'.join(lines))
    (out/'SHA256SUMS').write_text(''.join(f'{dataset.digest(p)}  {p.name}\n' for p in sorted(out.iterdir())
                                        if p.is_file() and p.name != 'SHA256SUMS'))
    print(json.dumps(summary['audit']))
    print(json.dumps(summary['full']['all'], indent=2))
    print(json.dumps(summary['timing']['all'], indent=2))
    print(json.dumps(summary['default_wait_policy']['all'], indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
