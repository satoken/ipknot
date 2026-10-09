# 既存の形状スコアの重みを調整し、独立データで評価する

projectedで定義した形状スコアは、crossingでもそのまま使える。今回の係数調整は、BPPモデル、候補閾値、交差制約、補助変数の定式化を変更しない。

\[
s(m)=\beta_0+\beta_s\min(A,B)
 -\beta_L\sum_{r=1}^3\log(1+L_r)
 +\beta_c I[L_2=0].
\]

全体倍率lambdaと各betaを同時に調整すると同じスコアを重複して探索するため、lambda=1に固定する。crossingでは、各交差辺に `max_m s(m)/(A*B)` を配分し、実際に選択された辺の寄与を既存の支持行に対応する連続補助変数で目的関数へ加える。配分と実現条件については [crossing.md](crossing.md) を参照。

betaはIPknot目的関数の単位で、kcal/molではない。とくにbeta_cは中央ループがゼロという形状の特徴であり、配列依存のcoaxial stackingエネルギーではない。

## 事前に固定した評価方法

原著者が公開した [IPknot++参照データ](https://zenodo.org/records/4923158) を使用する。出典は Sato and Kato, *Briefings in Bioinformatics* 23(1), bbab395 (2022), [doi:10.1093/bib/bbab395](https://doi.org/10.1093/bib/bbab395)。参照データのライセンスはCC-BY-4.0。ZIPのMD5は `01f491361b31b4520f5fd527d4ad5cbf`。

- 12〜200塩基の単一配列のみを対象とし、曖昧な塩基を含む配列と同一配列の重複を除く。BPSEQの添字・相手の範囲・相互性を検証する。
- 調整用はfrom_bpRNA-1mのPK陽性141配列とPKフリー141配列、合計282配列。独立評価用はfrom_Rfam14.5のPK陽性300配列とPKフリー300配列、合計600配列。
- PKフリー配列は、対応するPK陽性セットと50塩基幅の長さ区分ごとの配列数を一致させる。抽出順は `SHA256(pk-score-v1:ZIP内パス)` で決め、予測結果は使わない。
- 調整用の17条件は [calibrate.py](calibrate.py) のgridに事前定義した。ゼロ重み、以前のbeta0=0.2/betaL=0.03とその弱い倍率、ループの強い減点、負の開始スコア、短い方のステム長の加点、中央ループゼロの加点を含む。
- 通常のFASTA入力でLinearPartition-C、二レベル、自動閾値、refinement=1。BPP候補閾値は元のまま。自動閾値選択のPKスコア倍率kappaは通常の1を使い、弱い二条件では0も調整用で比較する。
- BPPの計算アルゴリズムは同じだが、refinementは最初の予測に依存するため、二回目のBPPの値まで同一とは限らない。固定BPPでのデコーダ単体比較ではなく、通常設定の最終予測を検証する。
- 調整用の全塩基対micro F1が基準版以上、かつPKフリー配列のPK誤検出数が基準版以下の条件から、PK塩基対micro F1が最大のものを選ぶ。同点では特徴係数の絶対値の和が小さい方を選ぶ。ゼロ重みも選択候補に含む。
- 勝者をJSONに固定してから独立評価へ進む。独立評価では基準版、ゼロ重み、以前のcrossing、以前のprojected、選択されたcrossingを比較する。基準版は試作前に保存した実行ファイル。
- 調整用は各一回、独立評価は逐次三回。時間は配列ごとの中央値の合計。上限10秒。失敗・タイムアウトを黙って除外しない。
- 塩基対は添字の完全一致で採点する。PK塩基対は、一つ以上の塩基対と交差する塩基対の集合。F1はmicro平均で、PKフリー配列で予測した交差対もFPに含める。
- 基準版との差の95%区間は、同じRNAを対にして5000回bootstrap、乱数seed=20261003。ZIPには家族IDがないため、家族単位のbootstrapではない。

調整用には非AU/CG/GUの参照対が764/6846対ある。独立評価用の参照対13844対にはそれらはない。参照は改変せず全対を採点する。この差、およびBPPに存在しない対をスコアだけでは復元できないことは、調整と評価の解釈に影響する。

## 調整用の結果と固定した係数

17条件×282配列=4794実行に失敗・タイムアウトはなかった。基準版では全塩基対F1=0.549493、PK塩基対F1=0.192986、PKフリー配列の誤PK予測は82/141配列だった。

| crossingの係数 | 全塩基対F1 | PK塩基対F1 | PKフリー配列の誤検出 |
|---|---:|---:|---:|
| beta0=0.2, betaL=0.03 | 0.534450 | 0.188100 | 91/141 |
| beta0=0.02, betaL=0.003 | 0.551197 | 0.199186 | 84/141 |
| beta0=-0.02, betaL=0, betaS=0 | 0.551714 | 0.190195 | 71/141 |
| beta0=-0.02, betaS=0.01, betaL=0.003 | 0.552419 | 0.197215 | 80/141 |
| beta0=-0.05, betaS=0.025, betaL=0.0075 | 0.550982 | 0.200178 | 79/141 |

最後の条件を選び、独立評価の前に [calibration-selected.json](measurements/calibration-selected.json) に固定した。betaC=0、lambda=1、kappa=1。式は

\[
s(m)=-0.05+0.025\min(A,B)
 -0.0075\sum_{r=1}^3\log(1+L_r).
\]

両ステムが2対なら開始項とステム項が相殺し、ループがある形状は減点される。両方が長いほど加点できるため、以前の正の開始項だけの例より短い候補を過剰に優遇しにくい。これは候補形状の重みの解釈であって、選択構造に完全なステムが必要になる変更ではない。

```sh
build/ipknot --pk-h-formulation crossing \
  --pk-h-intercept -0.05 --pk-h-stem-reward 0.025 \
  --pk-h-loop-penalty 0.0075 sequence.fa
```

調整用の改善は選択に使用したデータ上の結果であり、汎化の証拠にはしない。全17条件の結果は [calibration-tune-summary.json](measurements/calibration-tune-summary.json)、事前定義は [calibration-protocol.json](measurements/calibration-protocol.json)、配列と参照対は [calibration-references.jsonl](measurements/calibration-references.jsonl)、各実行は [calibration-tune.jsonl](measurements/calibration-tune.jsonl) に保存した。

## 独立600配列の結果

ARM64、HiGHS 1.15.1、1スレッド。各条件で600配列×3回、合計9000実行。失敗・タイムアウトはゼロで、各配列の三回の予測はすべて一致した。ゼロ重みの予測は600配列すべてで基準版と一致した。

| 条件 | 全塩基対F1 | PK塩基対F1 | PK precision | PK recall | PKフリー配列の誤検出 | 600配列の合計時間 |
|---|---:|---:|---:|---:|---:|---:|
| 基準版 | 0.635633 | 0.244491 | 0.245388 | 0.243601 | 143/300 | 69.096 s |
| ゼロ重み | 0.635633 | 0.244491 | 0.245388 | 0.243601 | 143/300 | 69.129 s |
| projected、beta0=0.2/betaL=0.03 | 0.619065 | 0.237970 | 0.201907 | 0.289718 | 175/300 | 75.875 s |
| crossing、beta0=0.2/betaL=0.03 | 0.617128 | 0.216071 | 0.186124 | 0.257502 | 175/300 | 154.017 s |
| crossing、調整済みの固定係数 | 0.634084 | 0.240930 | 0.239381 | 0.242498 | 136/300 | 118.754 s |

調整済み条件は146/600配列で予測が変わった。PKフリー配列での誤検出は減ったが、PK誤対の総数は3395→3492、正解PK対は1104→1099となり、PK塩基対F1は改善しなかった。全塩基対F1も微減した。

基準版との差のpaired RNA bootstrap 95%区間は、全塩基対F1で `[-0.004949,+0.002009]`、PK塩基対F1で `[-0.014060,+0.007012]`。点推定は低下し、改善の証拠はない。家族ごとの依存を扱った区間ではないため、この区間を家族間の汎化保証に使わない。

同じ特徴式の係数を変えるだけで実装できたが、今回の有限の探索では独立データの精度向上を確認できなかった。推奨する非ゼロのデフォルト値は得られていないため、デフォルトはゼロのままにする。上のCLIは調整用で選ばれた条件の再現用であり、精度向上を保証する推奨設定ではない。

### 計算量と実測時間

調整済み条件で最大111個の連続変数、375行の追加制約、1784個の重み付き交差辺。新しい整数変数はゼロ。補助変数と行数の線形性は保たれているが、通常設定の合計時間は基準版の約1.72倍となった。以前の強いcrossing重みでは約2.23倍だった。projectedは追加変数・行がなくても約1.10倍となった。線形個の補助変数と同じ計算次数だけでは、求解時間が不変であることを保証できない。

全条件の指標は [calibration-test-summary.json](measurements/calibration-test-summary.json)、差の区間と予測変更数は [calibration-test-comparisons.json](measurements/calibration-test-comparisons.json)、全実行は [calibration-test.jsonl](measurements/calibration-test.jsonl)。実行ファイル・参照データ・事前定義のSHA256と環境情報は [calibration-fingerprint.json](measurements/calibration-fingerprint.json) に保存した。

## 補足確認：既存NMRと実在する長鎖RNA

重みを再選択せず、以前から使用しているNMR 16配列と、公開Rfam参照データの実在する長鎖RNA二配列も逐次三回ずつ比較した。合計108実行で失敗・タイムアウトはゼロ。長鎖RNAは501〜800塩基の候補から `SHA256(pk-score-long-v1:ZIP内パス)` 順でPK陽性とPKフリーを一つずつ選んだ。NMRセットは過去に結果を見ているので、独立のブラインド評価として扱わない。

| NMR 16配列 | 全塩基対F1 | PK塩基対F1 | 正解PK対 | 誤PK対 | 合計時間 |
|---|---:|---:|---:|---:|---:|
| 基準版 | 0.826021 | 0.170213 | 8/41 | 45 | 0.883 s |
| 調整済みcrossing | 0.834094 | 0.177778 | 8/41 | 41 | 1.554 s |

この小さなセットでは誤対を減らしてF1が上がったが、正解PK対数は変わらない。独立600配列での改善を示す結果にはならない。

| 実在するRfam配列 | 塩基数 | 基準版の時間 | 調整済みcrossingの時間 | 連続補助変数 | 追加行 |
|---|---:|---:|---:|---:|---:|
| URS0000D6935B_12908_1-615、PK陽性 | 615 | 1.021 s | 1.223 s | 48 | 155 |
| AEYC01000066.1_4818-4311、PKフリー | 508 | 3.639 s | 4.835 s | 120 | 388 |

615塩基の配列では全塩基対F1が0.475862→0.479167、PK塩基対F1が0.123711→0.116505。508塩基の配列の予測は同一で全塩基対F1=0.045977、PK誤対86対。参照精度自体が低い例を含む二例だけで、長鎖への汎化を主張しない。実行時間はそれぞれ約1.20倍、1.33倍で、線形個の補助変数でも同じ予測を得るまでの時間は増え得る。

結果と参照は [calibration-external-separated.json](measurements/calibration-external-separated.json)、[calibration-external-references.jsonl](measurements/calibration-external-references.jsonl)。再実行には [check_external.py](check_external.py) を使う。

HiGHSで既存テストを含む54 CTestが通った。実データのFASTA検証で、RF00005は10本の72〜82塩基の配列であり、739塩基の一配列ではないことを確認した。以前の合成入力の測定と説明は [README.md](README.md) に訂正を記録した。

## 再実行

```sh
curl -L --fail https://zenodo.org/api/records/4923158/files/ipknot_dataset.zip/content \
  -o /tmp/ipknot-dataset.zip
python3 experiments/pk-score/calibrate.py prepare
python3 experiments/pk-score/calibrate.py tune
python3 experiments/pk-score/calibrate.py test
python3 experiments/pk-score/check_external.py
```

生の結果、予測、ログは `results/calibration/` に保存する。同じ実行ファイル・データ・事前定義を使う場合は中断後に再開できる。変更を検出した場合は、結果を混ぜず新しい出力先を要求する。
基準版は試作前の `7d560bb` のソースに対応する保存済みバイナリを使用した。別環境ではその基準版をビルドし、`--baseline /absolute/path/to/ipknot` で指定する。スコア実装は `a55d476` 以降のソースにある。

## 適用範囲

これはPK陽性を増やした短いRNAの部分集合での検証であり、自然なPK出現率、全公開データ、長鎖RNAを代表する結果ではない。調整と評価はデータソースを分け、完全一致配列の重複を除いたが、今回新たに80%配列同一性で家族をクラスタリングしたわけではない。有限の17条件の探索であるため、この形状モデル全体の最適な係数を保証しない。

形状候補の上限は既存の12塩基のステム・30塩基の各ループのまま。crossingは選択された交差辺のスコアを厳密に加えるが、その重みの元になったステム全体・非対合ループが成立することは要求しない。完全なHモチーフの物理エネルギーを評価した実験ではない。
