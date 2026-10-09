# オリジナルDPの線形補正・局所期待値・DD候補集合の比較

3方式をnative実装し、重みを従来282配列だけで調整してから、既知900配列と新規300配列で評価した。
新規300で局所Boltzmann/crossingのPK対F1は0.24494→0.24799だったが、差の95%区間は0を含む。
線形補正は改善せず、projectedへ切り離すと局所モデルでもPK対F1が下がった。
候補集合方式は追加DD/BPP求解なしで約2%の時間増に収まるが、PK再現率が落ちる。
いずれも任意指定の実験機能とし、既定値は変更していない。

## 共通の条件とエネルギー

- LPC100の初期BPPを全条件で共有、refinementは0。
- DD improved-beam100、最大50反復、patience0、crossing-beam100、witnesses16。
- 配列文字・塩基種・GC含量などはスコア、学習、モデル選択に渡さない。
  新規RNAの配列は既存LPCの入力と完全一致の除外用ハッシュに限って使用する。
- 原著DPのパラメータを固定し、変換の少数の係数だけを282配列で調整した。
  設定固定後に新規300を抽出した。既知900は工学的比較として扱う。
- 各集合はPK陽性150、陰性150、長さ12–200nt。新規300は以前の2,400 IDと完全一致配列を除外し、
  陰性を50nt幅の長さ階級で対応させた。近縁配列や同じファミリーの完全除外は保証しない。

[Dirks–Pierce (2003), Eq.16](https://www.its.caltech.edu/~niles/chebe163/papers/dirks03.pdf)の
PK固有ループ項は開始9.6 kcal/mol（外側）、15.0（多分岐/入れ子）、境界ステム0.1、未対合塩基0.1。
単純な外側H型の二本のステムでは `G=9.8+0.1*U` とした。
`RT=0.6163314008174533 kcal/mol`、`q=G/(RT)`。
DP09ではない。配列依存の通常のnearest-neighbor項を加えた完全自由エネルギーでもない。

1・2は候補ブロックのループ**長の和**をUのproxyにするため、ループ内の対合を正確には数えない。
3は実際に選択された対合からUを数える。交差する最大ステムの連結成分を一つのPKループとして扱い、
外側/包囲された文脈と非交差の子ステムを数える。単純なH型のPK固有項には一致するが、
複雑なPKや入れ子の交差成分の分解は近似である。未対応の複雑な形状をコスト0に戻す方式は採らない。

## 1・2の実装と重み

候補ステム長をA,B、各ステムの初期BPPの平均をa,bとする。
1を旧実装の固定モチーフ重みだけで比較すると長さ配分も変わるため、
新しい1・2では同じ `1/(A*B)` 配分を使い、全ステムの塩基対数に対応する尺度で比較した。
旧DPのモチーフ重み固定方式も別の対照として残した。

```text
J = eta - tau*q
1: Phi = kappa*(A+B)*J
2: Q = [(1-a)*(1-b), a*(1-b), (1-a)*b, a*b*exp(J)]
   P = Q/sum(Q)
   Phi = kappa*[A*(P10+P11-a) + B*(P01+P11-b)]
crossing contact coefficient = Phi/(A*B)
```

2は安定したlogit/sigmoidと代数的に等価な期待値式を使い、指数のoverflowを避ける。
ステムのBPP平均は行列ごとのキャッシュで再利用する。
1・2とも新しい整数変数・支持制約を加えない。DDのcrossingでは既存接点にスコアを付け、
projectedでは各上位対の目的係数へ潜在相手の最大スコアを移す。projectedの支持候補は元のまま。
DDのwitnesses16ではcrossingスコアが証人の順位にも影響するため、残る証人や実行可能集合が
baselineと完全に一致するとは限らない。projectedは元の証人順位を保つ。

a,bはステム全体の成立確率のproxyで、独立な `a*b` をPKフリー分布から得た真の交差同時確率とは呼ばない。
2の期待値は局所的な全か無かのステムモデルに対するもの。部分ステムや重なるモチーフ、
projectedの相手切り離しにより全体の厳密なMEAにはならず、グローバルなBPPも置き換えない。
温度37℃は固定、tauは統計的なtempering。etaによる開始項の相殺を物理的導出とは扱わない。

同じ3,160個の数値ブロック例を使い、各RNA内の重みを揃え、正解の対合占有率の積をsoft labelとした。
`logit(a*b)+eta-tau*q` に二係数のridge logistic fittingを行う。
既存281 identity群の5 foldでridge .01/.1を比較し、.01を選択。
held-out local loglossは0.25103（.1では0.26446）。1・2は同じ学習済みJを使う。
その後、fold外のnative予測でkappaを選ぶ。
全対F1を下げず、陰性RNAの偽PK数を増やさない候補の最大PK対F1を採用し、補正0も候補に含めた。
同じfoldでridgeとkappaを選ぶため、この値は調整用の推定であり独立な評価ではない。

```text
eta = 1.8523222707481422
tau = 0.043312241784810926
kappa linear = .05 (crossing/projected)
kappa local  = 1   (crossing/projected)
```

linearのkappa候補は `.0005,.002,.01,.05`、localは `.025,.1,.4,1`。
変換後のスケールが異なるため共通のkappaを強制しない。全選択設定は非ゼロだった。
原著の係数は固定だが、この選択自体は経験的な較正である。

## 3の実装と追加費用

自動閾値探索で既に得ているDD最良解だけを保存し、物理構造を重複除去する。
DDの回復手順、反復、beam、候補グラフを変更せず、追加求解もしない。
今回はRNAごとに10閾値解、重複除去後は1–10構造、中央値5構造だった。
候補が1構造の6 RNAでは選び直しても構造は変わらない。

```text
F_common(S) = sum_selected_pairs alpha[level]*(initial_BPP - .25)
log W(S) = F_common(S)/T_score + eta*N_PKcomponents - tau*G_loop(S)/(RT)
P(S) = softmax(log W(S))
p_tilde(i,j) = sum_S P(S)*I[(i,j) in S]
selected = argmax_existing_S sum_(i,j in S) (p_tilde(i,j)-t_MEA)
```

全候補で同じ基準閾値.25を使う。生成元の異なる閾値の目的値をそのまま比較しない。
同一物理構造の別レベル分解は最良の共通目的値を採用し、閾値探索の訪問回数を重みにしない。
構造の混合なので各塩基の周辺確率の和は1以下。
これは有限候補集合に条件付けたハイブリッド分布で、全RNA構造のDP分配関数ではない。
MEAも候補内に限定するため追加DD求解がなく、候補外の新しい対合は作れない。

old282で候補集合方式内の設定を比較し、`tau=.01, eta=.15900536605796906, T_score=4, t_MEA=.5` を固定した。
etaは外側H型のU=0を基準にした値。グリッドは測定[protocol](measurements/dp-conversion-protocol.json)に記録した。
**この方式内の最良設定でも、old282のPK対F1はbaselineの0.18493に対して0.14631だった。**
比較のために最良の候補集合設定を評価したのであり、baselineを上回る採用候補としては選んでいない。

追加作業は候補の保存、重複除去、ループ成分の採点、確率集計、候補内MEA。
周辺確率集計は保存した対合数に線形。今回の単純な交差成分検出はステム数Sに対してO(S²)で、
候補数Kを合わせてO(K*S²)の項を持つ。追加費用が完全に0、あるいは全追加処理が線形とは主張しない。
この長さ・候補数ではその費用は小さかった。

## 新規300の精度と時間

全対/PK対F1は対合数を合計したmicro F1。偽PKは陰性150 RNAのうちPKを一つ以上予測したRNA数。
同じ初期BPPから3回逐次実行し、RNAごとに条件順を回転させた。
時間は各RNAの3回の中央値の合計をbaseline=1とした値。
初期BPP計算は含まないので、end-to-end FASTA速度の測定とは区別する。

| 条件 | 全対F1 | PK対F1 | 偽PK /150 | 時間比 |
|---|---:|---:|---:|---:|
| baseline | .63587 | .24494 | 72 | 1.000 |
| 1 linear crossing | .63485 | .23818 | 75 | 1.106 |
| 1 linear projected | .63647 | .24007 | 75 | 1.011 |
| 2 local crossing | .63633 | .24799 | 75 | 1.103 |
| 2 local projected | .63706 | .23857 | 75 | 1.011 |
| 3 finite ensemble | .64065 | .20267 | 34 | 1.022 |
| 3 同じ温度/MEAでeta=tau=0 | .63892 | .20027 | 39 | 1.009 |
| 旧DP固定モチーフ重み crossing | .63464 | .24262 | 72 | 1.098 |
| 2 物理因子 eta=0,tau=1,kappa=.1 | .63714 | .24100 | 70 | 1.110 |

2 crossingのPK対F1差は+0.00305、対応するRNA/メタデータ群のbootstrap 95%区間は
`[-.01357,+.01844]`。精度改善を確認できたとは言えない。
2 projectedとcrossingの差は−0.00943、95%区間 `[-.02289,+.00068]`。
軽くはなるが同等精度だとも断言できない。
3のbaselineとの差は−0.04227、95%区間 `[-.08005,-.00774]`。
偽PKが減る一方、真のPKも取り落とす。
エネルギーを外した候補集合でも大きく下がっており、低下全体をDP項だけに帰せない。

新規300のbaselineは3,000グラフ・44,089反復。1 crossingは47,664、2 crossingは47,550反復で、
それぞれ102,995/102,055の非ゼロ接点を持った。
反復数は約8%増え、さらに各接点の積因子処理が加わる。
1・2 projectedは42,812/43,026反復で、接点積因子の採点は0。
3はbaselineと同じ3,000グラフ・44,089反復。
3の時間比95%区間は `[1.0009,1.0463]`、2 crossingは `[1.0722,1.1351]`。
baselineの合計時間は2.024秒、3は2.069秒、2 crossingは2.231秒だった。

## 既知900とエネルギーだけを外す対照

既知集合は各300を1回ずつ評価した工学的比較。速度の主評価は3回測定した新規300に置いた。

| PK対F1 | 旧300 | 後続300 | 前回の新規300（現在は既知） |
|---|---:|---:|---:|
| baseline | .22855 | .23214 | .24583 |
| 1 crossing | .24000 | .23343 | .24373 |
| 1 projected | .23327 | .22662 | .24422 |
| 2 crossing | .23254 | .23685 | .25061 |
| 2 projected | .23077 | .23134 | .24544 |
| 3 finite ensemble | .19609 | .20690 | .17751 |

既知900でも2 crossingは各集合で小幅増だが、各差の95%区間は0を含む。
3は一貫して低い。新規300の結果を見て設定を選び直していない。

重み・温度・interceptを固定したままenergy slopeだけを0にする追加対照も行った。
ここには再学習・再調整はない。新規300のPK対F1は：

| 条件 | interceptだけ | DP energy slopeあり | 差の95%区間 |
|---|---:|---:|---:|
| 1 crossing | .23644 | .23818 | [-.00798,+.01258] |
| 2 crossing | .24234 | .24799 | [-.01132,+.02393] |
| 3 finite ensemble | .21780 | .20267 | [-.03710,+.00235] |

2のDP固有項だけの精度寄与も有意とは確認できない。
3では固定した正のpriorにenergy penaltyを戻すとPK対F1が下がる傾向がある。
自由なinterceptを含む全スコアの変化とDP項そのものの効果を区別する。
bootstrapは5,000回、複数方式の探索的比較で多重比較補正はしていない。
新規300の群はaccessionメタデータであり、family単位の独立性を保証する区間ではない。

## 検算、ファイル、実行

- CTestは81/81成功。数値テストは独立な4状態列挙3,000例、ゼロ補正、極端なlog-weight、
  符号付きCLI、singleton、重複除去、塩基容量、MAPとMEAの差、DPループ項を確認。
- RNA上で1,908グラフ、1,982ブロック係数を独立計算し、最大誤差1.94e−16。
  24,100回復構造の容量・同レベル非交差・下位レベル支持・目的値を検算した。
- 48 RNAの候補集合を出力し、共通BPP目的値、実現ループ項、確率、期待利得、選択を独立に検算。
  全1,200 RNAで3のDDグラフ/反復列がbaselineと一致。既知900のbaselineは過去の予測と完全一致。
- 主測定は16,200 native実行。失敗0、再試行や除外なし。数値を含む全RNA記録は
  [gzip JSON](measurements/dp-conversion-per-rna.json.gz)へ損失なく圧縮した。

モデルは `models/dp-loop-linear-v1.txt`, `models/dp-loop-local-v1.txt`。
両方式のモデルはheaderだけが異なり、eta/tauは同じ。
使用例（初期BPP入力、refinementなし）：

```sh
build/ipknot -x -r 0 -t auto,auto --decoder dd --dd-dp improved-beam \
  --dd-beam 100 --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  --pk-h-allocation blocks --pk-h-max-stem 200 --pk-h-max-loop 200 \
  --pk-h-formulation crossing --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/dp-loop-local-v1.txt \
  --pk-learned-scale 1 input.bpp

build/ipknot -x -r 0 -t auto,auto --decoder dd --dd-dp improved-beam \
  --dd-beam 100 --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  --pk-ensemble --pk-ensemble-scale .01 --pk-ensemble-intercept .15900536605796906 \
  --pk-ensemble-temperature 4 --pk-ensemble-threshold .5 input.bpp
```

1はモデルをlinearへ変更しscaleを.05、projectedはformulationを変更する。
3はDD専用。固定閾値では候補が一つなので予測を変更しない。
`--pk-ensemble-output FILE` で配列を含まない数値JSONLを出力できる。時間測定では出力しない。

[評価スクリプト](evaluate_dp_conversion.py)は `prepare → tune → fresh → cache → evaluate → summarize`。
既存のold282キャッシュ/数値ブロックと `/tmp/ipknot-dataset.zip` が必要。
モデル・実装・キャッシュをhashで固定し、変更がある場合は新しい出力先を使う。
[独立監査](audit_dp_conversion.py)の `audit`、同じscriptの `ablate`、
[圧縮と時間統計](export_dp_conversion.py)の順に実行する。
詳細は[測定結果](measurements/dp-conversion-summary.json)、
[独立監査](measurements/dp-conversion-native-audit.json)、
[energyだけの対照](measurements/dp-conversion-ablation.json)、
[反復/候補数](measurements/dp-conversion-runtime-audit.json)に記録した。
