# DDでのprojected PKスコア

**DDでも `--pk-h-formulation projected --pk-h-allocation blocks` を使える。**
既存の塩基対目的係数へ補正を加え、交差支持行は元のまま保つ。crossing型の積因子は
作らない。形状、DP、CC、学習済みBPP補正、形状/BPP統合を利用でき、固定対・NMR制約、
refinementの経路にも対応する。

同じ形状/BPP統合重みを用いた600配列の比較では、crossingより約9%速く、PK対F1は
両標本でわずかに高かった。全対F1はわずかに低く、精度差の95%区間はいずれも0を
含むため、精度の優位は確認できていない。

## 係数と交差条件

[crossingで使ったモデルと重み](dd_integration.md)を変更せず、固定最大候補ステム
ブロックA、Bの接点係数を計算する。

```text
s = -0.05 + 0.025 min(A,B) - 0.0075 Σ log(1+Lr)
w_AB = s/(A B) + 0.05 (θ・f)/anchor_length
```

ブロックAの上位レベルkの各候補対uについて、下位レベルlの相手ブロックBからの提案は
`w_AB × N_B,l`。`N_B,l` はBの当該レベルの全候補対数であり、最終的に選択した対数や、
支持行に保持した対数ではない。各lで採点された相手の最大提案を取り、lごとの値を
合計してuの目的係数へ加える。BからAへの提案も同じ規則で処理する。

全提案が負の場合は負の最大値を保つ。提案がない場合、または接点係数が0の場合は
補正しない。学習済みモデルでは、固定対がBPPに存在しないブロックの学習済み成分を
使わず、hybridの形状成分は利用する。crossing側もILPと同じ扱いに揃えた。

ブロック対は既存の交差beamの走査中に採点し、全ブロック対を総当たりしない。
`--dd-witnesses` で絞る前に採点するため、支持相手数を変えても係数は変わらない。
ただし `--dd-crossing-beam` は採点するブロック対プールに影響する。同じプールなら
ILPの固定ブロックprojectedと同じ係数になるが、通常の有界DDがILPの全プールを
必ず使うわけではない。`--dd-crossing-beam 0` は全プールの診断用で線形性を保証しない。

支持行は補正前の塩基対係数で構築する。したがってprojectedを有効にしても既存の
支持行の相手は変わらず、選択された上位対は各下位レベルの選択対と実際に交差する。
一方、採点の根拠になった最良の相手ブロックが選ばれる条件までは課さない。相手が
部分的に成立した場合も、上位対に配分された候補全体の補正が入る。

追加の変数・支持行・積因子を作らず、補正後の係数をDP、回復、交換、上下界に共通して
使う。係数生成にはブロック対の採点とキャッシュが必要だが、固定された最大ステム長、
beam、支持相手数、レベル数、反復数のもとで期待O(n+M)を保つ。

## 精度・時間

以前に評価した二つの300配列を再使用した。各標本はPK陽性150・陰性150、12–200nt。
新しい独立検証ではなく、DDへの対応確認と探索的な比較である。モデルの再学習やDD
向けの調整はしていない。主比較は以前の重み `gamma=1, lambda=.05`。以前ILPで試した
半分の重み `gamma=.5, lambda=.025` も副比較として測定した。

初期LPC100のBPPは全条件で同じキャッシュ、refinementはオフ。改良DD beam100、
最大50反復、patience0、交差beam100、各支持行16相手、自動閾値auto/autoを使う。
時間は条件順を回した逐次3反復のRNAごとの中央値の合計で、BPP計算を含まない。
CLI起動、BPP読み込み、出力は含む。aarch64、Release、HiGHS、ViennaRNA 2.7.2、
CPU affinityは固定せず許可された0–19番を使用した。

先の300配列。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK | 時間 |
|---|---:|---:|---:|---:|
| DD、スコアなし | 0.624238 | 0.228545 | 65 | 3.613秒 |
| crossing、統合 | 0.627324 | 0.236665 | 65 | 3.929秒 |
| projected、統合 | 0.626984 | 0.237342 | 64 | 3.574秒 |
| projected、半分の重み | 0.625350 | 0.230734 | 64 | 3.617秒 |

後の300配列。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK | 時間 |
|---|---:|---:|---:|---:|
| DD、スコアなし | 0.645019 | 0.232143 | 81 | 3.762秒 |
| crossing、統合 | 0.644453 | 0.234266 | 79 | 4.091秒 |
| projected、統合 | 0.644127 | 0.237619 | 80 | 3.723秒 |
| projected、半分の重み | 0.644022 | 0.233188 | 80 | 3.757秒 |

projected−crossingのPK対F1差は先で `+0.000677`、95%区間
`[-0.006307,+0.007521]`、後で `+0.003354`、区間 `[-0.004551,+0.014382]`。
全対F1差は先で `-0.000341`、後で `-0.000327`。予測が変わったのは先で27本、後で22本。
各標本の対応付き80%配列同一性クラスタbootstrapを5000回行った探索的な区間である。

時間比projected/crossingは先で0.910、95%区間 `[0.898,0.922]`、後で0.910、
区間 `[0.895,0.925]`。600本の合計はcrossing 8.020秒、projected 7.298秒。
スコアなしDDの7.375秒と近く、この標本ではprojectedによる時間増は観測しなかった。
目的係数が変わると停止までの反復も変わるため、採点処理自体が無料という意味ではない。
FASTAからのBPP計算込みで同じ時間比になるとは限らない。

半分の重みでは両標本のPK対F1がcrossingより低く、このDD比較では元の重みを上回らなかった。
独立データで再検証する前に精度を保証する設定とは扱わない。
[集計・係数監査](measurements/dd-projected-summary.json)、
[配列ごとの測定](measurements/dd-projected-per-rna.json)。

## 検証

- 正負の形状・DP・CC・BPP・hybridを2/3レベルで検査。全プールでILPの既存係数と一致。
- 交差beam0/1/100、支持相手数1/無制限で、projectedの支持行がスコアなしと一致。
  支持相手数を変えてもprojected係数は一致。部分成立と負の補正も検査。
- 全600配列で、重み0のprojectedはスコアなしDDの予測・グラフ・目的値・反復・停止理由と一致。
- 基準DDとcrossingの600配列×3反復は、前回の予測と各グラフの数値に完全一致。
- 24本のR0と4本のR1、5,325反復状態で、排他性、レベル内非交差、スタック、各下位
  レベルの交差支持、符号付き補正込みの目的値を独立に再計算。最大目的値誤差は `1.42e-14`。
  R0の16,575個の係数はこの標本ではILP全プール計算と一致し、最大差は `1.39e-17`。
  これは別の配列でも有界なプールが全プールと一致する保証ではない。
- 固定対・NMR制約のlinear/full経路、追加回復・交換・上下界、R1、CLIを含む79件のCTestが成功。

追加診断では交差beamを無制限にした4本のR0について、5,470個の係数がILPの全プール
計算と一致した。562・1000・2883ntの3本でもLPC100、R1、改良DD100、最大50反復で
実行した。562・2883ntで非ゼロの選択PK補正を確認し、1000ntではこのH型スコアの
採点対象がなく補正0だった。長鎖は各1回の動作確認であり、精度・時間の比較標本ではない。
[追加診断記録](measurements/dd-projected-extra.json)。

## 使用方法

```sh
build/ipknot --decoder dd --dd-dp improved-beam --dd-beam 100 \
  --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  -e lpc --beam-size 100 -r 0 -t auto,auto \
  --pk-h-formulation projected --pk-h-allocation blocks \
  --pk-h-max-stem 200 --pk-h-max-loop 200 --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/hybrid-r0-v1-posterior.txt \
  --pk-learned-scale 0.05 --pk-hybrid-shape \
  --pk-h-intercept -0.05 --pk-h-stem-reward 0.025 \
  --pk-h-loop-penalty 0.0075 --pk-h-weight 1 input.fa
```

`-r 1` ならrefinementあり。既定のDD formulationは引き続きcrossing/blocks。
ログの `DD PK projection` に採点ブロック対数と非ゼロの係数数を出力する。
`Selected PK score` は選択された対の補正を合計した値である。

以前のBPP・中立特徴キャッシュを用いた再現手順。

```sh
python3 experiments/pk-score/dd_projected.py run
python3 experiments/pk-score/dd_projected.py audit
python3 experiments/pk-score/dd_projected.py summarize
python3 experiments/pk-score/dd_projected_extra.py
```

バイナリ、ソース、モデル、入力のSHA256と固定設定は集計記録に保存した。
