# 配列を特徴に使わないPKスコアの実装と精度評価

提案した4係数のモデルを実装し、refinementなし、DD改良beam 100・LPC 100・最大50反復で評価した。
調整用282配列で選択したlog-projectedは、新規300配列のPK対F1を0.24583から0.24739へ
上げたが、差の95%区間は−0.00700〜+0.00949で、改善を確証できなかった。
既存の形状＋BPP統合projectedの0.25011にも届かない。新モデルを既定値にはしない。

crossingとCCの調整では補正0が選ばれた。効果を確認する非ゼロ比較も実施した。
log-crossingの最良非ゼロ設定は新規300配列で0.25153と有望な点推定だが、
差の区間は0を含む。これは主評価の後に定義した探索的な追加比較であり、
元の選択結果を置き換えるものではない。

## 実装したスコア

入力は初期BPP、塩基対座標、ステム長、ループ長だけ。
塩基の種類、GC率、k-mer、配列embedding、RNA名・家族名、絶対位置を特徴にしない。
学習時もBPPファイルの塩基列と既存exportの配列欄を読み飛ばす。
配列を使うのは従来のLPC計算と評価データの完全一致重複を除くハッシュだけである。
BPPに含まれる配列由来の情報は利用する。既存の12特徴モデルも配列を直接特徴にしていない。

ステムXの長さをX、塩基対確率をpとし、各対の両端に接する別の塩基対のうち
最大のBPPをqとする。自分自身を除外し、閾値で候補を切る前の疎BPP全体からqを求める。

```text
h(p) = p/(1+p)
e_X = mean_X h(p)
m_X = mean_X [h(p)-h(q)]
E = min(e_A,e_B)
M = min(m_A,m_B)
C_log = sum_r log(1+L_r)
F = -b0 + bE*E + bM*M - bC*C
S = clip(F,-T,T)
w_AB = S/(A*B)
```

4係数は非負。競合を固定した場合、自分のBPPを上げて補正が減ることを防ぐ。
これらはPKの同時成立確率やエネルギーではなく、既存目的関数への修正量である。
CC版はCだけを `max(0, native CC06/DP-fallback q)` に置き換える。
既存の幾何依存ループコストを使い、配列依存の塩基対・スタッキングエネルギーや
CC09の三次接触表は加えない。

crossingでは、保持した接点に `w_AB*x_u*x_v` を加える。全A×B接点が選ばれれば
ブロック対の補正はS、部分的な成立では接点の選択率に応じた補正となる。
DDが保持しなかった接点には課金しない。

projectedでは同じw_ABから、上位対の各下位レベルlについて
`max_B(w_AB*N_B,l)` を既存の単項係数へ加える。
Nは相手ブロックのそのレベルの全候補対数で、実際に選ばれた対数ではない。
負の最大値は負のまま保ち、スコア0の相手は候補から外す。
既存の交差支持制約は維持するが、採点に使った最良相手の成立を要求しない。
各下位レベルについて単項補正の絶対値はT/A以下、上位ブロック合計はT以下。
構造全体や全ブロック対の合計をT以下に制限するモデルではない。
`--pk-learned-scale`を1以外にすると、この上限もその絶対値で拡大する。

新しいモデル形式は `IPKNOT_PK_BOUNDED_V1`。
必須項目は `bias, support, competition, loop_cost, cap`、任意項目は
`loop_model log|cc`（省略時log）。不正値・重複・欠落はエラー、cap=0は補正を無効化する。
従来の `IPKNOT_PK_LINEAR_V1` の計算は維持した。
正規化だけを変える対照用に `IPKNOT_PK_LINEAR_AB_V1` も追加した。
bounded形式と従来のhybrid形状項・特徴exportの併用はエラーにする。
BPPにない強制対を含むブロックでは、bounded補正を棄権する。

## 学習・設定選択とデータ

学習は従来のbpRNA調整用282配列だけで行った。PK陽性141、陰性141、
提供済みの80%同一性クラスタ281群を5foldに分けた。
refinementなしの中立的な候補exportから3,160ブロック対を得た。
候補があるRNAは178、ないRNAは104。候補を持つRNA内の重みの和を1にする。

教師値は `(正解対を含む割合_A)*(正解対を含む割合_B)` である。
部分的な正解ステムを扱うためのsoft occupancyで、物理的な同時成立確率ではない。
訓練側RNAの教師値平均から固定した事前log-oddsをoffsetとし、
その上の補正Fを非負係数のロジスティック回帰で学習する。
事前offsetはnativeの補正Fに加えない。4係数すべてにL2正則化を置き、
16通りの制約面を列挙して解き、KKT誤差が1e−7以下であることを確認した。

正則化候補は0.01と0.1。5foldの保留ブロック対loglossで両モデルとも0.01を選択した。
次に、RNAごとにそのfoldを除外して学習したモデルでDD予測し、
Tを `[0,0.025,0.05,0.1,0.2]` から選んだ。
全対F1と陰性RNAの偽PK数を悪化させない条件の下で、PK対F1が最大のTを選ぶ。
補正0は必ず比較に含め、同点では小さいTを選ぶ。
正則化とTの選択に同じ5foldを使うため、これは調整の推定値であり、
完全なnested CVによる独立な検証値ではない。

| モデル | b0 | bE | bM | bC | crossingの選択T | projectedの選択T |
|---|---:|---:|---:|---:|---:|---:|
| log | 0 | 0.149595 | 0.534538 | 0.003262 | 0 | 0.025 |
| CC | 0 | 0.177604 | 0.503611 | 0.003016 | 0 | 0 |

選択前の基準PK対F1は0.18493、選択log-projectedは0.18841。
それ以外の非ゼロ設定は基準を上回らなかった。
全係数とTを固定してから評価した。評価900配列で再学習・再調整はしていない。

評価は次の3集合。それぞれPK陽性150、長さを50nt区間で合わせた陰性150である。

- previous300：以前から使っているRfamの300配列。
- later300：前回追加したRfamの300配列。
- new300：今回の設定固定後に、過去のPK実験で使っていないRfamから新しく抽出した300配列。

すべて12〜200nt。new300は過去の全PK実験manifestのIDと完全一致配列ハッシュを除いた。
塩基配列による類似度計算は追加していない。利用可能なデータに家族IDがないため、
**近縁RNAの排除や家族を分離した外挿性能は保証できない**。
既知600配列は回帰・技術比較、新規300配列が今回の保留評価である。
新規300のbootstrapは配列類似度ではなくaccession/URS識別子の接頭辞299群を使う。
その区間は家族内相関を十分扱えていない可能性がある。

## 主評価

同じRNAでは全条件で同一の初期LPC100 BPPキャッシュを使い、`-r 0`。
DDは `improved-beam 100, max-iter 50, patience 0, crossing-beam 100, witnesses 16`。
自動閾値 `auto,auto`、固定最大ステム・ループ200、選択用PK重み0、learned-scale 1。
基準、既存hybrid 2方式、AB正規化対照2方式、選択モデル4設定、
事前に指定した最大非ゼロT=0.2の診断4設定、合計13条件を比較した。

全対・PK対F1はmicro集計。PK対は、同じ構造内で別の対と交差するすべての対と定義する。
偽PK RNA数は、正解にPKのない150 RNAのうち、PKを予測したRNAの数である。
13条件×900配列×3回=35,100実行がすべて成功し、予測とDDグラフ統計は3回一致した。

新規300配列の主な結果：

| 条件 | 全対F1 | PK対F1 | 偽PK RNA /150 | DD等の合計時間(s) |
|---|---:|---:|---:|---:|
| 補正なし | 0.62621 | 0.24583 | 61 | 1.7983 |
| 既存hybrid crossing | 0.62674 | 0.24792 | 59 | 1.9749 |
| 既存hybrid projected | 0.62784 | 0.25011 | 59 | 1.7829 |
| 既存hybrid AB対照 crossing | 0.62621 | 0.24660 | 58 | 1.9863 |
| 既存hybrid AB対照 projected | 0.62620 | 0.24960 | 58 | 1.8008 |
| 新log projected、選択T=0.025 | 0.62636 | 0.24739 | 59 | 1.7857 |
| 新log crossing / CC両方式、選択T=0 | 0.62621 | 0.24583 | 61 | 1.816〜1.822 |

log-projected−補正なしのPK対F1差は+0.00156、
5,000回の対応付きbootstrap 95%区間は **[−0.00700,+0.00949]**。
全対F1差は+0.00014、区間 **[−0.00319,+0.00330]**。予測が変わったのは29 RNA。
既存hybrid projectedとの差はPK対F1で−0.00273、区間 **[−0.01271,+0.00539]**。

| 既知集合 | 補正なしPK対F1 | 既存hybrid projected | 新log projected |
|---|---:|---:|---:|
| previous300 | 0.22855 | 0.23734 | 0.23201 |
| later300 | 0.23214 | 0.23762 | 0.22595 |

既知2集合でも新projectedの変化は一様ではなく、両差の区間は0を含む。
新規300の長さ別PK対F1は次のとおり。小さい区分の記述統計であり、
この表を使って設定は選んでいない。

| 長さ(nt) | RNA数 | 補正なしPK対F1 | 既存hybrid projected | 新log projected | 追加log crossing |
|---|---:|---:|---:|---:|---:|
| 1–50 | 4 | 0.20690 | 0.20690 | 0.20690 | 0.20690 |
| 51–100 | 188 | 0.28818 | 0.29187 | 0.28571 | 0.28392 |
| 101–150 | 78 | 0.27030 | 0.26761 | 0.27070 | 0.27649 |
| 151–200 | 30 | 0.08673 | 0.10804 | 0.10230 | 0.11503 |

T=0.2を固定した事前指定の診断では、新規300のPK対F1は
log-crossing 0.24223、log-projected 0.24154、CC-crossing 0.22249、
CC-projected 0.23335。強く補正すれば良くなるという結果ではなかった。
これらはTの選択結果ではない。

## 最良非ゼロ設定の追加比較

主評価の後、最大Tだけの診断とは別に、元の調整結果で最良だった**適格な非ゼロT**を
選択T=0の3系列で評価した。log-crossing T=0.025、CC-crossing T=0.025、
CC-projected T=0.05。CC-projected T=0.025は元の調整で全対F1等の条件を満たさない。
評価集合の精度でTを選び直していない。ただし追加比較の発案は主評価の後であり、
**独立な確認実験・新しい採用設定とは扱わない**。元の選択T=0は保存した。

同じ900配列、同じBPPを用い、基準・既存hybrid両方式・主設定log-projectedも再実行した。
7条件×900×3回=18,900実行が成功。再実行対照の予測とDDグラフは主評価と一致した。

| 非ゼロ診断 | T | previous300 PK対F1 | later300 PK対F1 | new300 PK対F1 | new300 全対F1 | new300 偽PK /150 |
|---|---:|---:|---:|---:|---:|---:|
| log crossing | 0.025 | 0.23339 | 0.23537 | 0.25153 | 0.62765 | 59 |
| CC crossing | 0.025 | 0.23226 | 0.23725 | 0.24471 | 0.62618 | 59 |
| CC projected | 0.05 | 0.22483 | 0.22824 | 0.24226 | 0.62483 | 58 |

log-crossingの新規PK対F1差は+0.00570、探索的95%区間 **[−0.00425,+0.01585]**。
CC-crossingは−0.00112 **[−0.01123,+0.00857]**、
CC-projectedは−0.00357 **[−0.01403,+0.00641]**。
各区間に0を含み、複数比較による調整はしていない。
log-crossingの点推定は今後調べる価値があるが、今回の結果から精度向上を確定しない。

## 時間・検算・解釈

時間は条件の順をRNAと反復ごとに回転した逐次実行で測定し、RNAごとの3回の中央値を合計した。
初期BPP計算を含めず、プロセス起動・入力・モデル読込・DD・出力を含む。
Release、HiGHS、aarch64、ViennaRNA 2.7.2。CPU affinityは固定していない。
過去の別の測定期間との絶対時間比較には使わない。

主評価log-projected / 補正なしの時間比は0.993、95%区間[0.976,1.011]で、ほぼ同じ。
追加比較log-crossingの時間比は1.088 [1.064,1.111]、CC-crossingは1.095 [1.071,1.118]。
projectedは既存係数だけを更新する。crossingはDDの接点部分問題とその乗数更新が残る。
それぞれ同方式の既存hybridと比べた時間比は約1であり、新特徴の計算が大きな負担になっていない。

top-two競合索引、ステム要約、ブロック対スコアを再利用する。
候補対数M、長さnに対して、固定beam・相手数・最大ステム長・レベル数・反復数の下では
期待O(n+M)を維持する。ILPやDDの最適化が同じ時間で終わることまでは保証しない。
`--dd-crossing-beam 0` 等の全候補診断はこの計算量の前提から外れる。

検証結果：

- native CTest **79/79成功**。新形式の読込、正負・ゼロ、log/CC、DD/ILP、2/3レベルを検証。
- 12 RNA×5設定=60 traceを別の数値計算で検算。11,811回復状態、27,697接点、
  35,205 projected係数を確認。接点最大誤差6.94e−18、projection 1.39e−17、目的値5.86e−14。
- その標本ではprojectedの全相手プール計算との差は0。
  有界DDプールと全プールが一般に一致する保証ではない。
- 新規300すべてでFASTAからのLPC予測と保存BPPからの基準予測が一致。
- 全900でcap=0の予測・グラフが基準と一致。既知600の基準・旧hybrid両方式は前回結果と一致。
- 主評価35,100＋追加18,900実行の全成功・3回一致、traceの最終予測と時間測定時予測の一致を確認。

min支持・競合差は、PKフリーBPPで弱く見える真のPKステムにも罰則を与える可能性がある。
AB正規化は、従来のアンカー正規化で長い相手へ加わっていた補正を弱める。
ABだけを変更した対照も既知600では旧projectedより低い点推定であり、
単調性や上限を置けば精度が上がるとは言えない。これらは結果の解釈で、原因を確定する実験ではない。
今回の上限200nt、PKを半数含む構成と近縁排除の限界を超えた一般化は検証していない。

## 保存物と再実行

主評価の[固定設定・学習・集計・区間](measurements/sequence-free-summary.json)、
[配列別測定](measurements/sequence-free-per-rna.json)、
追加比較の[設定・集計](measurements/sequence-free-best-nonzero-summary.json)と
[配列別測定](measurements/sequence-free-best-nonzero-per-rna.json)、
[長さ別集計](measurements/sequence-free-length-breakdown.json)、
[数値監査](measurements/sequence-free-audit.json)を保存した。
各集計には実行ファイル、ソース、入力、モデルのSHA256を含めた。
生ログ・キャッシュ・traceは `results/sequence-free/` と `results/sequence-free-best-nonzero/` にある。

保存モデル：

- [log T=0.025](models/sequence-free-log-v1-T0.025.txt)：選択projected、追加crossing診断で共通。
- [log T=0](models/sequence-free-log-v1-T0.txt)：crossingの元の選択結果。
- [CC T=0](models/sequence-free-cc-v1-T0.txt)：CC両方式の元の選択結果。
- [CC T=0.025](models/sequence-free-cc-v1-T0.025.txt)、[T=0.05](models/sequence-free-cc-v1-T0.05.txt)：追加診断。
- [旧モデルのAB正規化対照](models/hybrid-r0-v1-posterior-ab.txt)：特徴・係数は旧モデルと同じ。

選択log-projectedを手元のFASTAで試す例（実験用。精度改善が確立した設定ではない）：

```sh
build/ipknot --decoder dd --dd-dp improved-beam --dd-beam 100 \
  --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  -e lpc --beam-size 100 -r 0 -t auto,auto \
  --pk-h-formulation projected --pk-h-allocation blocks \
  --pk-h-max-stem 200 --pk-h-max-loop 200 --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/sequence-free-log-v1-T0.025.txt \
  --pk-learned-scale 1 input.fa
```

評価再実行には従来実験のmanifest/export/BPPキャッシュと公開データarchiveが必要。
[sequence_free_evaluate.py](sequence_free_evaluate.py) は入力を数値特徴へ変換し、
学習・選択・新規抽出・キャッシュ・測定・集計を別段階で行う。
入力や設定が変わればfingerprintが拒否するため、各段階の`--output`に新しい結果ディレクトリを指定する。
下記は既定ディレクトリへ初めて実行する場合のコマンドで、測定済みディレクトリで学習をやり直さない。
学習時点のスクリプトそのものも生結果に `training-script.py` として保存した。

```sh
c++ -std=c++17 -O2 -Isrc experiments/pk-score/sequence_free_geometry.cpp \
  src/pk_energy.cpp -o /tmp/ipknot-sequence-free-geometry
python3 experiments/pk-score/sequence_free_evaluate.py prepare
python3 experiments/pk-score/sequence_free_evaluate.py train
python3 experiments/pk-score/sequence_free_evaluate.py tune
python3 experiments/pk-score/sequence_free_evaluate.py fresh
python3 experiments/pk-score/sequence_free_evaluate.py cache_fresh
python3 experiments/pk-score/sequence_free_evaluate.py test
python3 experiments/pk-score/sequence_free_evaluate.py summarize
python3 experiments/pk-score/sequence_free_audit.py
# 元の調整から最良非ゼロTを選ぶ探索的な追加比較
python3 experiments/pk-score/sequence_free_nonzero.py run
python3 experiments/pk-score/sequence_free_nonzero.py summarize
```

モデル・監査・長さ集計の公開用コピーは
[sequence_free_postprocess.py](sequence_free_postprocess.py) が作成する。
これは学習や設定選択をせず、測定結果・3回一致・SHA256を検査する。
