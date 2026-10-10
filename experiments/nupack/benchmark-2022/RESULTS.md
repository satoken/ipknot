# NUPACK・LinearNUPACK・LPC・LPV ベンチマーク

状態: 計測中。10,772/47,372 runs。

公式2022年IPknot単一配列セット全11,843配列（bpRNA-1m 3,437 / Rfam 14.5 8,406、12–4,381 nt）。
NUPACK・LinearNUPACKはrefinement 0、LPC・LPVはrefinement 1。全方式でDD beam 100 / 50 iterations、crossing beam 100 / witnesses 16、自動閾値。線形エンジンのbeamは100。
通常NUPACKは1配列60秒、他方式は1,800秒。全方式の仮想アドレス空間上限は8 GiB。曖昧塩基は位置を保持して対を作らない。

## 完了率と全件の状況

| 方式 | 成功 | 時間切れ | メモリ上限 | その他エラー |
|---|---:|---:|---:|---:|
| NUPACK (r=0) | 0 | 991 | 1698 | 0 |
| LinearNUPACK (b=100, r=0) | 2693 | 0 | 0 | 0 |
| LPC (b=100, r=1) | 2695 | 0 | 0 | 0 |
| LPV (b=100, r=1) | 2695 | 0 | 0 | 0 |

線形3方式の全11,843配列の精度は [LINEAR_ACCURACY.md](LINEAR_ACCURACY.md) に先行集計済みです。以下の実行時間は元のX925 phaseのみです。


## 線形3方式の共通成功集合（2,692配列）

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 2692 | 0.408517 | 0.126519 | 2.5019 | 32151.0 | 173.8 / 4591.0 |
| LPC (b=100, r=1) | 2692 | 0.466341 | 0.103928 | 0.4987 | 6119.1 | 10.9 / 236.4 |
| LPV (b=100, r=1) | 2692 | 0.443365 | 0.069524 | 0.5534 | 6768.9 | 8.1 / 220.8 |

## 4方式の共通成功集合（0配列）

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|

4方式の共通集合の長さ範囲: [None, None] nt。通常版が完了した短い配列に偏るため、この集合を線形版の全体精度と混同しない。各集合内では全方式が同じ配列で比較される。

## 線形共通集合: 短鎖（150 nt以下）

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 143 | 0.494632 | 0.341389 | 0.8520 | 121.4 | 86.0 / 103.9 |
| LPC (b=100, r=1) | 143 | 0.541431 | 0.149144 | 0.1835 | 24.3 | 7.1 / 10.3 |
| LPV (b=100, r=1) | 143 | 0.524529 | 0.182381 | 0.1863 | 25.3 | 6.1 / 7.1 |

## 線形共通集合: 中鎖（151–500 nt）

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 1733 | 0.423181 | 0.151759 | 1.7034 | 3753.1 | 136.2 / 471.8 |
| LPC (b=100, r=1) | 1733 | 0.488294 | 0.132409 | 0.3647 | 766.5 | 8.9 / 22.9 |
| LPV (b=100, r=1) | 1733 | 0.472591 | 0.089077 | 0.3909 | 827.6 | 7.1 / 17.3 |

## 線形共通集合: 長鎖（500 nt超）

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 816 | 0.362282 | 0.087026 | 29.9096 | 28276.5 | 1456.0 / 4591.0 |
| LPC (b=100, r=1) | 816 | 0.406557 | 0.061910 | 5.7999 | 5328.4 | 76.0 / 236.4 |
| LPV (b=100, r=1) | 816 | 0.367074 | 0.039621 | 6.3963 | 5916.0 | 64.5 / 220.8 |

## 線形共通集合: Rfam14.5

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 1303 | 0.441051 | 0.157368 | 1.7034 | 3694.7 | 134.5 / 866.6 |
| LPC (b=100, r=1) | 1303 | 0.516847 | 0.127141 | 0.3610 | 725.0 | 8.9 / 45.5 |
| LPV (b=100, r=1) | 1303 | 0.496747 | 0.076164 | 0.3913 | 796.5 | 7.2 / 36.5 |

## 線形共通集合: bpRNA-1m

| 方式 | 完了数 | Macro F1 | PK F1 | 実時間中央値 (s) | 合計実時間 (s) | RSS中央値 / 最大 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| LinearNUPACK (b=100, r=0) | 1389 | 0.377997 | 0.084984 | 5.3409 | 28456.3 | 363.1 / 4591.0 |
| LPC (b=100, r=1) | 1389 | 0.418961 | 0.072672 | 1.2013 | 5394.1 | 17.8 / 236.4 |
| LPV (b=100, r=1) | 1389 | 0.393288 | 0.060584 | 1.3007 | 5972.4 | 13.7 / 220.8 |

## 指標・測定条件

Exact pair equality; per-sequence macro F1; empty/empty F1=1; undefined PPV/SEN=0. PK F1: crossing pairs on reference-PK sequences. Failed runs are excluded from common sets; explicit all-attempt failure-zero macro also recorded.

One fresh process per input/condition, single thread; includes posterior, DD and specified refinement. Wall sums are sums of child times, not elapsed benchmark duration. Shared-host parallel measurements; the auxiliary accuracy-only phase uses disjoint A725 cores during part of the run, and its timings are not merged. Censored exact failures are not successful latencies.

CPU: Cortex-X925の10コアに固定（5–9,15–19）。追加の精度評価は別のA725コアを使用し、その測定時間はこの表に含めない。同じ配列の各方式は同じCPUで実行し、方式順を入力ハッシュに基づきシャッフル。GCC 13.3、Release -O3 -DNDEBUG、ENABLE_ILP=OFF。

精度で配列を選択せず、パラメータ調整やBPPキャッシュも行わない。通常版の打ち切りはその構造が誤りという意味ではなく、指定資源内で結果が得られなかったことを示す。

通常NUPACKには既存のパラメータ読み込み・疑似結び目周辺確率の修正が入る。LinearNUPACKは制限付き文法とbeam pruningを使う近似であり、幅100でも通常版と一致する保証はない。

生記録: `raw-results.jsonl.gz`。配列別結果: `per-rna.csv`。集計・paired bootstrap: `summary.json`。条件・ハッシュ: `final-manifest.json`, `provenance.json`。再実行方法: `README.md`。

## グラフ

グラフは全件の測定が完了してから生成します。
