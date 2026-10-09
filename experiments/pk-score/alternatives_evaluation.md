# 短いステム区間・排他性補正・最良相手・構造候補の順位学習

四つの案を実装し、crossing/projectedで使える条件を比較した。新規300 RNAでは、
最良相手crossing＋順位学習がPK対F1最大の0.24242（baseline 0.22064）だった。
ただしPKなし150 RNAの偽PKは66→79件、全対F1は0.65414→0.65008。
PK対F1差の95%区間も0を含むため、安定した改善や既定値への採用は支持しない。
DP付き排他性補正crossingはPK対F1 0.23541、偽PK69件、時間1.097倍。

## 条件と分割

- refinement=0、初期BPP=LPC100。DDはimproved-beam100、最大50反復、patience0、
  crossing-beam100、witnesses16。全条件でBPP、候補cut、基本制約を共有する。
- スコア・学習への入力はBPP、競合BPP、対の座標、長さ、構造だけ。
  A/C/G/U、GC、k-mer、配列embedding、RNA名・家族名を特徴にしない。
  配列文字は従来LPCへの入力と完全一致除外のハッシュに限って使う。
- 調整用は従来の282 RNA（PK141、陰性141）、281 identity群、既存5 fold。
  DPのeta/tauは以前の同じfoldのモデルを使用し、held-out RNAの係数をそのRNAでfitしない。
  幅・倍率・順位学習の正則化は調整用OOF結果で選ぶため、nested CVの独立推定ではない。
- 119条件×282 RNA＝33,558回のnative調整。
  最大PK対F1を、全対F1≧baseline・陰性偽PK数≦baselineの条件で選ぶ。補正0も残す。
- 順位学習は12系列の候補を出力する3,384回のnative実行と、数値特徴だけの学習。
  推論時のDD/BPP呼出し数は増えない。
- 全設定と主比較を凍結後、新規300 RNA（PK150、長さ帯を合わせた陰性150）を選択。
  過去3,000 IDと完全一致ハッシュを除外。候補母集団にはPK249 RNAが残っていた。
  近縁・家族の独立性は保証せず、bootstrapはaccessionメタデータ群による。
- 以前の最後の300 RNAは既知の工学的回帰確認。新規300は条件順を回転する逐次3反復。
  16条件、合計19,200回。失敗・除外・再試行は0。
  時間はRNAごとの中央値の和で、BPP計算を除き、起動・モデル読込・閾値探索を含む。
  5,000回のpaired bootstrap。副比較には多重比較補正をしていない。

## 実装したスコア

### 1. 短い区間による採点

`--pk-core-width 2/3` は、候補の最大連続ステムを重ならない短い区間に分ける。
末尾に1対だけ残る場合は直前へ結合するので、最大幅は指定値＋1。
塩基対候補・変数は増やさず、各対の所属区間は入力から一意に決まる。
全区間の列挙や、選択後に形状を取り直す処理は行わない。
ILPとDDのcrossing/projectedで共通に使用できる。

各区間のBPP平均、長さ、ループ形状で既存の局所DPスコアを計算し、
区間長A×Bで配分する。区間分割は適用範囲と長さ正規化の両方を変える。
局所のループ項を複数区間対へ付けるため、構造全体の厳密なDPエネルギーにはならない。

### 2. 排他性を除く局所モデル

`IPKNOT_PK_EXCLUSION_V1`。候補区間の平均BPPをa,bとし、観測分布の00/10/01を
`1-a-b, a, b` と置く。本来独立な潜在事前から11を除いたという仮定の下で、
11の重みを復元し、Jで傾ける。

```text
J = eta - tau * G_loop/(RT)
Q = [1-a-b, a, b, a*b*exp(J)/(1-a-b)]
P = Q / sum(Q)
Phi = kappa * [A*(P10+P11-a) + B*(P01+P11-b)]
w_AB = Phi/(A*B)
```

`a*b=0` または `1-a-b<.01` なら適用を見送る。
log-weightで計算し、極端なJで指数がoverflowしないようにする。
core平均はステム成立確率のproxyであり、この式はPKフリーBPPから真のPK確率を
同定するものではない。期待対数の局所増分を目的係数に使う近似である。

flatはeta=tau=0。DPはオリジナルDirks–Pierce 2003のループ項だけを使い、
以前の調整結果eta=1.8523222707481422、tau=.043312241784810926、37℃を共有する。
RNAごとの調整実行では対応するfold外fitを使う。配列依存NNエネルギーは加えない。
新しい潜在事前についてeta/tauを再fitした比較ではなく、kappaを別途調整した比較。
crossingは実際の選択接点、projectedは以前と同じ潜在最良相手へのunary配分に使う。

### 3. 実際の最良相手

`--pk-best-partner`。上位対uと各下位レベルの支持行について、

```text
x_u * max_{B: B内に選択相手あり} sum_{v in B} w_uv*x_v
```

を加える。未選択の相手区間は最大値に参加しない。選ばれた相手に負のスコアしか
なければ負の最大値を使う。選ばれた未採点区間の重み0は参加する。
接点ごとの積因子を支持行の一つの因子へ置き換え、その因子内で交差支持条件も課す。
下位対コピーのunary重みの正部分を集計し、各相手群の最良非空選択を比較するので
局所求解はO(d＋群数)。2^d列挙は不要で、作業バッファも反復間で再利用する。

今回はDD crossingだけで実装。通常のILP、NMR制約付きDD、additive-contact用の
global/joint bound・exchange・追加recoveryとは混用できないことを明示的に検査する。
通常DDの基本制約・beam・反復上限は同じ。zero scoreでは従来の支持行へ戻る。

### 4. 実現可能な構造候補の順位学習

`--pk-rank-model` と `--pk-rank-scale`。
既存auto閾値探索の最大10構造を、元の選択値＋lambda×8特徴の内積で選ぶ。
教師は局所ステムの占有率ではなく、その構造の
`全対F1 + PK対F1 - .1*陰性RNAの偽PK`。
RNA内の候補比較を均等化したpairwise squared-hingeとL2正則化で重みをfitする。

特徴はpseudoF、pseudoPKF、塩基対の個数/配列長、PK対の個数/配列長、PK対平均BPP、PK対平均競合差、
DP構成要素ループ項/(10×長さ)、PK補正値/長さの8個。
fold外で学習した重みによる調整用結果から、正則化.001/.01/.1/1、
lambda=0/.25/.5/1/2/4を選ぶ。lambda=0では元の選択と一致する。
共通の候補を学習して選ぶ方式で、候補集合を確率混合した以前のMEA方式とは異なる。

## 新規300 RNA

短い区間の二条件は、調整で補正0が選ばれたため、最良の非ゼロ診断条件を載せる。
他の新スコアは制約を満たした調整用勝者。順位学習は非ゼロが選ばれた4系列だけを
追加実行し、残り8系列のlambda=0は親と同じ条件へまとめた。

| 条件 | 全対F1 | PK対F1 | 陰性150の偽PK | 時間/baseline |
|---|---:|---:|---:|---:|
| baseline | .654143 | .220635 | 66 | 1.000 |
| 従来局所DP crossing k=1 | .653267 | .231975 | 69 | 1.104 |
| 従来局所DP projected k=2 | .648292 | .218757 | 69 | 1.008 |
| core3 crossing k=1（診断） | .650344 | .229157 | 72 | 1.168 |
| core3 projected k=1（診断） | .646405 | .215208 | 75 | 1.041 |
| 排他flat crossing k=.4 | .653783 | .227516 | 69 | 1.106 |
| 排他flat projected k=1 | .647630 | .218330 | 70 | 1.011 |
| 排他DP crossing k=.4 | .653772 | .235408 | 69 | 1.097 |
| 排他DP projected k=.4 | .649097 | .220086 | 71 | 1.011 |
| 最良相手・局所DP crossing k=2 | .651229 | .234780 | 69 | 1.094 |
| 最良相手・排他flat crossing k=1 | .651315 | .228392 | 69 | 1.097 |
| 最良相手・排他DP crossing k=.4 | .650327 | .225542 | 69 | 1.091 |
| 従来局所DP projected＋rank lambda=1 | .648023 | .227979 | 74 | 1.014 |
| 排他flat projected＋rank lambda=.5 | .647612 | .222222 | 73 | 1.010 |
| 排他DP projected＋rank lambda=.5 | .649077 | .223955 | 74 | 1.020 |
| 最良相手・局所DP crossing＋rank lambda=.25 | .650080 | .242424 | 79 | 1.096 |

調整用で凍結した主比較は最後の行。PK対F1差+.02179の95%区間は
[-.000380, +.045534]、全対F1差-.00406の区間は[-.012389, +.004234]。
PK対recallは.21142→.27371、precisionは.23069→.21755。
陰性で新たに偽PKになる15 RNAと偽PKが消える2 RNAがあり、偽PKが13件増える。
順位学習単独の親に対するPK対F1差の区間も[-.012168, +.030170]で0を含む。

排他DP crossingは、全対F1をほぼ保ちながらPK対precision/recallが両方上がる。
ただしPK対F1差+.01477の95%区間は[-.002007, +.032910]で0を含む。
既知300では.23578→.25372だが、これを新しい独立確認とは扱わない。

## 同じ倍率での最良相手の対照

最良相手の勝者k=2とflat k=1は新規集合選択前に調整用で決まった値。
同じ倍率の通常crossingを追加し、スコア型の差を倍率の差から切り分けた。
この対照を見て主比較・設定を変更していない。
初回の2対照は2,400実行、その後の通常/最良相手の交互実行は4,800実行。
初回の集計もresults内に保存し、精度の分母から失敗を除く処理は行っていない。

| 新規300、同じk | 通常crossingのPK対F1 | 最良相手のPK対F1 | 最良相手−通常の95%区間 |
|---|---:|---:|---|
| 局所DP k=2 | .235838 | .234780 | [-.006934, +.004941] |
| 排他flat k=1 | .231758 | .228392 | [-.009016, +.002053] |
| 排他DP k=.4 | .235408 | .225542 | 主表の同じ倍率の比較 |

最良相手だけの独立した精度改善は見られない。
速度も親と対照を同じバッチで交互に実行した比較では、大きな差は確認できない。
主表とは別バッチの時間を割った最初の集計は比較用に使わず、
[交互実行記録](measurements/alternatives-matched-interleaved.json)を使用する。

## 採点範囲と実行時間

調整用の正解交差接点8,204個のうち、両対がいずれかの元の候補cutを超えるものは
1,405個（17.13%）。形状の潜在的被覆は最大ステム1,340個（16.33%）、
core2 1,387個（16.91%）、core3 1,381個（16.83%）。
前の設計説明で重視した「長い最大ステムのため採点から落ちる」問題の寄与は小さく、
この設定では候補化されない正解対が主な制限だった。
これらは物理的な潜在被覆で、レベル・保持witnessの同時実現可能性や予測recallの上限ではない。

新規300の合計DD反復はbaseline 39,162、従来局所crossing 42,105、
排他DP crossing 41,629、最良相手crossing 40,644、core3 crossing 45,922。
core3は採点接点を85,161→137,723へ増やすので、定数幅でも時間が増える。
順位学習は親と反復・グラフが完全に同じで、追加の時間比は1.00163
（95%区間[.98674,1.01414]）。crossing側の約10%の増分を順位学習のコストと混同しない。

## 検証・成果物

- CTest 83/83成功。局所因子5,000例の全列挙、潜在事前逆変換5,000例、
  小さいDDモデル40例の全列挙による実現可能性・目的値・上界を検査。
  core0/2/3、2/3レベル、通常・局所DP・排他モデルのILP/DD係数整合性を検査。
- 実RNA2,544グラフ、27,702復元状態、10,563係数を独立計算。
  最大係数誤差7.22e-16。容量、同レベル平面性、stack、全下位レベルの交差支持を検査。
- 数値rank候補2,560個を検査。新規300×3反復×4 rank系列＝3,600実行で
  親とDD統計が完全一致。既知300のbaseline・従来crossing/projectedの900実行が
  以前の予測・DD統計と一致。
- [全結果と凍結設定](measurements/alternatives-summary.json)、
  [実行由来](measurements/alternatives-evidence.json)、
  [独立検算](measurements/alternatives-audit.json)、
  [回帰・時間区間](measurements/alternatives-regression.json)、
  [RNA別の圧縮記録](measurements/alternatives-per-rna.json.gz)。
- [凍結したモデル](models/alternatives-v1/README.md)。すべて任意指定で、既定値は変更しない。

## 再実行

既存old282キャッシュと `/tmp/ipknot-dataset.zip` が必要。NumPyのあるpython3.12を使用。
過去のmanifestが増えればfreshの完全一致除外集合も変わるため、保存済みfresh manifestを
そのまま使用することが今回の評価の再現条件になる。

```sh
python3.12 experiments/pk-score/evaluate_alternatives.py prepare
python3.12 experiments/pk-score/evaluate_alternatives.py tune --jobs 4
python3.12 experiments/pk-score/evaluate_alternatives.py rank
python3.12 experiments/pk-score/evaluate_alternatives.py fresh
python3.12 experiments/pk-score/evaluate_alternatives.py cache
python3.12 experiments/pk-score/evaluate_alternatives.py evaluate
python3.12 experiments/pk-score/evaluate_alternatives.py summary
python3.12 experiments/pk-score/audit_alternatives.py
python3.12 experiments/pk-score/matched_alternatives.py
```
