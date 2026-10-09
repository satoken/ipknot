# IPknotの双対分解デコーダ

DD beam 200で反復上限50/100を全11,843配列で比較した結果は[beam 200の性能評価](dual-decomposition-beam200-performance.md)を参照。

長鎖での精度低下の原因診断とbeam幅の改善実験は[長鎖精度の検討](dual-decomposition-long-accuracy.md)を参照。

2022年論文の全公開単一配列データによる精度・速度・メモリの比較は[性能評価](dual-decomposition-performance-2022.md)を参照。

整数ギャップを縮める端点クラスタ上界と局所交換、その安全性・速度比較は[整数ギャップの検討](dual-decomposition-integer-gap.md)を参照。

上界・下界の追加改善と、その証明・実測比較は[上界・下界改善の検討](dual-decomposition-bounds.md)を参照。前回の実装・収束実験は以下に保存する。

`codex/pk-score` の `b10fa2f` から分岐した実装。ILPバックエンドは残し、
既定の `--decoder dd` でソルバーを呼ばない経路を選ぶ。`--decoder auto` はリンク済みの
ILPソルバーがあればILP、なければDDを選ぶ。

## 参考にした実装

[DAFS](https://github.com/satoken/dafs) の手元の `wip/linearization`、コミット
`1624ed1` を参照した。特に `src/nussinov.cpp` の `LinearNussinov::decode_impl`、
`src/gradient_manager.cpp`、`src/dd_recovery.h`、`src/dd_certificate.h`、
`docs/dafs_dd_acceleration_2026-09-19_ja.md` と
`docs/dafs_linfold_optimization_round4_2026-09-19_ja.md` を読んだ。
公開リポジトリの版と手元の線形化ブランチは区別する。

採用した考え方は、疎な候補の索引を反復中に再利用すること、右端ごとの区間開始点を
ビームで制限するNussinov、DPの作業領域の再利用、反復中に実行可能な主問題解を
保存すること、ビーム値と有効な上界を区別することである。

DAFSの追加高速化実験は一律には成功していない。射影ノルムのpilotは全体速度改善の
根拠にならず、column boundも採用ゲートを通っていない。この実装では容量制約の
大部分がゼロ乗数・余裕ありになり得るため、外向き成分を除くノルムを既定値とした。
DAFSからの性能の引き継ぎを主張せず、この問題で測定する。
`--dd-projected-norm=false` で全成分のノルムにも切り替えられる。

## 主問題とラグランジュ緩和

候補塩基対を `p=(i,j,l)`、選択を `x_p` とする。元の単項スコアは
`w_p = alpha[l]*(BPP(i,j)-threshold[l]) - 1[l>0]*pk_level_penalty`。
各レベルの `x` は、塩基を重複使用しない非交差構造である。

レベルをまたぐ制約は以下の2種類である。

1. 各塩基 `i` の使用数 `d_i(x) <= 1`。
2. 上位候補 `u` は**すべての**下位レベル `b<level(u)` に交差相手を持つ：
   `x_u <= sum_{v in W(u,b)} x_v`。

非負の乗数 `lambda_i`、`mu_ub` により、最大化主問題のラグランジアンは

```text
L(x,lambda,mu) = sum_p w_p*x_p
  + sum_i lambda_i*(1-d_i(x))
  + sum_ub mu_ub*(sum_{v in W(u,b)} x_v-x_u)
```

となる。単項係数を更新すれば、各レベルのDPは独立に解ける。最小化する双対側の
劣勾配は `g_lambda=1-d`、`g_mu=sum(x_v)-x_u`。
更新は `q <- max(0,q-eta*g)` で、違反すると乗数が増える。

## PKスコア

crossing / blocksとprojected / blocksを使える。候補の和集合から最大の連続ステムを
一度作り、その形状、DP/CCループエネルギー、または学習済みモデルを評価する。
固定ブロックの幾何スコアは `Phi/(A*B)`、従来の `IPKNOT_PK_LINEAR_V1` 学習済み補正は
アンカーの長さで割る。hybridでは両者を加算する。
新しい `IPKNOT_PK_BOUNDED_V1` は、BPP支持・競合との差・ループコスト・切片の
非負4係数で補正を作り、`[-T,T]` に制限してからA×Bで割る。
`loop_model log|cc`で対数ループ長または配列に依存しないCC06/DP-fallbackコストを使う。
塩基配列自体を特徴にせず、crossing/projectedへ同じ接点係数を渡す。
`IPKNOT_PK_LINEAR_AB_V1` は旧12特徴の正規化だけをA×Bへ変更する対照用形式。
新モデルのrefinementなし900配列の精度・時間評価は
[検証報告](../experiments/pk-score/sequence_free_evaluation.md)を参照。改善を確証できず、既定値にはしない。

保持した接触 `(u,v)` のスコアは `c_uv*x_u*x_v`。コピー `z_u,z_v` を導入し、
`z_u=x_u`、`z_v=x_v` を符号自由の乗数 `rho_u,rho_v` で緩和する。
各接触部分問題は次の4値の最大を選ぶだけである。

```text
(0,0): 0
(1,0): rho_u
(0,1): rho_v
(1,1): c_uv + rho_u + rho_v
```

DPには `-rho_u`、`-rho_v` を加える。劣勾配は `z-x`、更新は符号自由。
負のPKスコアもこの4状態で処理する。存在しない接触には課金しない。
`--pk-learned-model`、`--pk-hybrid-shape`、`--pk-energy-model`、
`--pk-level-penalty`、単純なH形状の各係数を利用できる。
DDではHスコアの既定formulation/allocationをcrossing/blocksにする。

`--pk-h-formulation projected --pk-h-allocation blocks` では、上位対の既存目的係数へ
符号付き補正を加える。ブロック対の接点係数を `w_AB`、ブロックBの下位レベルlの
候補対数を `N_B,l` とすると、Aの上位対に対するlからの提案は `w_AB N_B,l`。
各lで採点された相手から最大の提案を取り、lごとの値を合計する。全提案が負なら
負の最大値を保ち、0へ切り上げない。提案がない場合は0。逆向きも同じ規則で計算する。

相手ブロック対は既存の `--dd-crossing-beam` の走査から取得し、支持相手数
`--dd-witnesses` で絞る前に係数を計算する。相手ブロックの候補対数は支持行に
残った対数ではなく、当該レベルの全候補数を使う。全ブロック対を総当たりしない。
同じブロック対プールならILPの固定ブロックprojectedと同じ係数になるが、有界な
DDのプールがILPの全候補プールと常に一致する保証はない。
`--dd-crossing-beam 0` は全プールの診断用で、線形計算量の保証を外す。

交差支持行は補正前の目的係数で構築し、projectedによって相手の順位や保持する行を
変えない。既存の各下位レベルとの実際の交差は必要だが、採点に使われた最良の相手
ブロックが選ばれるとは限らない。積因子・追加変数・追加支持行は作らず、補正後の
係数をDP、回復、交換、上下界のすべてで使う。固定・NMR制約のlinear/full経路と
refinementにも対応する。固定されたbeam、支持相手数、最大ステム長、レベル数、
反復数のもとで期待O(n+M)を保つ。[検証と測定](../experiments/pk-score/dd_projected.md)。

exact/substems、supported、rerank、PK feature exportはDDでは未対応で、
明示的なエラーになる。ILPのcrossing auxiliaryの正規化・簡略化フラグはDDの積因子には
影響しない。ブロック幾何と係数が同じでも、保持する交差接触が異なる場合は総スコアも
変わる。全交差接触を維持したという保証はしない。

`PKPosteriorContext` の確率とステム要約のキャッシュは、順序を使わない検索だけなので
ハッシュに変更した。既存ILP側の特徴、スコアと出力の回帰テストも通す。

## 通常版と孤立対禁止版のDP

`F(i,j)` は完成した構造、`P(i,j)` は外側の塩基対が親のスタックから支援を受けても
よい構造とする。孤立対を許す場合は通常のNussinovで、完成対は
`w(i,j)+F(i+1,j-1)`。

孤立対を禁止するときは

```text
P(i,j) = w(i,j) + max(F(i+1,j-1), P(i+1,j-1))
S(i,j) = w(i,j) + P(i+1,j-1)
F(i,j) = max(F(i,j-1), F(i,k-1)+S(k,j))
```

を使う。`P` は存在する候補対にだけ持ち、欠損した内側対は使えない。負の内側対も
有利な外側対を支える場合は残す。空区間と隣接対も、入力候補にあれば扱える。
`-i` で孤立対を許せる。

ビームの区間内部は、直前列の開始点が `i+1` 以上の最良構造を使う。
先頭の未使用塩基を省くDAFSのsuffix queryを用いて、区間長全体の走査を避ける。
トレースは整数索引の非再帰処理なので長いステムでも呼び出しスタックを消費しない。

従来ILPの「隣接位置にも同じ向きの対がある」制約より、この連続塩基対の定義は厳しい。
バルジなどの扱いを含め、完全一致を目標にしていない。

## 線形性を保つ候補制限

疎なBPPが線形サイズでも、交差関係は二次個存在し得る。すべての交差辺を明示的に
持つモデルと、一般入力での線形計算量を同時には主張できない。この実装は以下を
明示的に制限する。

- `--dd-dp beam`：従来のbeam探索。`--dd-dp nussinov`はbeam探索を使わない厳密区間DP。
- `--dd-beam`：beam方式の右端ごとのDP区間開始点。既定100。nussinov方式では参照しない。
- `--dd-crossing-beam`：左端のsweep中の未閉鎖候補、各レベルに既定100。
- `--dd-witnesses`：上位対・下位レベルの1行あたりの交差相手、既定16。

同じ左端の候補は全照会の後で追加し、共有端点を交差と誤認しない。
active beamは単項スコア、witnessは単項スコア＋PK補正を優先し、同点は入力索引で
決める。捨てた候補と辺の件数はログに出す。ビームから捨てた候補の対自体はDPに
残り得るが、交差支援の辺は復活しない。

配列長 `n`、レベル込み候補数 `M`、レベル数 `K`、各固定ビーム幅と反復予算 `T` を
固定すると、前処理、1反復と実行可能解復元は `O(n+M)`、全DDは
`O(T*(n+M))`、記憶量も `O(n+M)` である。ハッシュ検索とselectionの平均計算量を
用いる。一般の密なBPPは `M=Theta(n^2)` なので配列長に対する線形性はない。
LinearPartitionなど固定幅の疎な生成器と組み合わせる必要がある。

3つの幅に0を指定すると制限を外す。短い入力での比較・検証用であり、線形性は失う。
DPの0は厳密Nussinov、交差の2つの0は全交差相手である。ただし有限反復の双対分解は
整数最適性を一般には保証しない。

自動閾値選択にも従来の全交差走査があった。DD経路では、固定幅active beamから作る
物理的な交差接触集合でPKのpseudo F値を評価する。これは従来のPK pseudo F値の
近似であり、閾値選択を変え得る。通常の塩基対pseudo F値は従来通り。
交差幅0ならこの評価も全交差を使う。

LinearPartition、閾値候補数、refinement回数も固定すれば、この経路は疎な入力規模に
対する線形計算量を維持する。ViennaRNA/CONTRAfoldなど別の確率生成器、無制限の
幅、可変反復予算には、この計算量の主張を適用しない。

refinement、アラインメントの確率平均とモデル混合も、行内の線形検索をハッシュ索引に
置き換えた。既存の候補の順序とfloatの加算順序は保存する。行あたりの候補数に
依存する二次時間を回避し、独立した旧方式とのビット単位の一致を検証した。

## 更新、上界、実行可能解

Polyak型の `eta=relaxation*(beam_L-best_primal)/||g||^2` を使う。
ビームの値は双対上界ではないため、主問題解以上にならない場合は減衰ステップに
切り替える。改善しない反復ではrelaxationも減らす。

ビームとは独立に、各塩基の正の単項係数の最大値の和を2で割った値をDPの上界に
使う。この値にラグランジアンの定数と4状態因子の最大を加えれば、現在の**保持した
交差グラフ**に対する上界になる。厳密DPではDP値そのものを使う。
グラフ制限前のILP最適値の上界としては使わない。ログは `graph_upper_bound` と
してこの範囲を示す。丸め誤差を考慮した形式的な区間演算証明ではない。

反復前に、交差支援と隣接スタックの不可能性をキューで伝播する。支持数がゼロになった
候補と、その候補を必要とする上位の対を取り除く。各候補・辺は一度処理するので
`O(M+E)` であり、ゼロになることが分かっているPK因子も取り除ける。

毎回、レベル0のDP解から始め、使用済みの塩基と下位の実際の交差支援を条件に、
上位のDPを再計算する。復元した実行可能構造のうち最良のものを返す。有限反復で
合意に至らなくても、塩基の重複、レベル内交差、下位レベルの支援不足を残さない。

停止は上界gap、劣勾配の停止、`--dd-patience`（既定0で停滞停止を無効）、
反復上限（既定50）。反復上限は閾値候補・refinementごとのsolveに適用する。
その他の既定値はDD beam 100、crossing beam 100、witnesses 16、
LPC beam 100、threshold `auto,auto`、refinement 1。
ビーム値の停留を整数最適性の証明とは解釈しない。停止名は `beam_stationary` /
`exact_stationary` に分ける。

`--dd-trace FILE` と `--dd-trace-state` で各反復を保存できる。
`--dd-unpruned-bound` は実際にDP枝刈りがなかったレベルの厳密値を上界に使い、
`--dd-schedule diminishing` は指数0.75の減衰ステップを使う。いずれも既定off。
整数ギャップ、双対最適化、回収、上界の余裕を切り分けた検査と推奨設定は
[収束検討](dual-decomposition-convergence.md)を参照。

構造固定、NMR個数、NMRスタック・バルジ・coaxial制約にも対応する。
既定の `--dd-constraints linear` はNMR候補・witnessの幅、伝播回数、回復状態数も
固定し、同一レベルの交差禁止はDPと線形のplanarity検査で守る。
NMR観測数・pattern長、全幅・反復/回復予算を固定すれば、制約構築と回復を含め
期待 `O(n+M)` 時間・線形記憶量を維持する。候補制限で実行可能解を失うことはあるが、
返した解のhard制約は検査する。`--dd-constraints full` は全制約モデルの診断用で、
線形保証はない。詳細・計算量の内訳・長さ倍増測定は
[制約付きDD](dual-decomposition-constraints.md)を参照。

## ビルドと実行

```sh
cmake -S . -B build-dd -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build-dd -j
ctest --test-dir build-dd --output-on-failure
build-dd/ipknot --decoder dd examples/drz_Ppac_1_1.fa
build-dd/ipknot --decoder dd -r 0 -t .5,.25 --pk-learned-model \
  experiments/pk-score/models/learned-r0-v1-all.txt examples/drz_Ppac_1_1.fa
```

ソルバー付きビルドなら同じ実行ファイルで `--decoder ilp` と `--decoder dd` を比較できる。
厳密な小規模診断は `--dd-beam 0 --dd-crossing-beam 0 --dd-witnesses 0
--dd-patience 0` を追加する。

## 検証

- 独立した物理マッチング全探索と、符号付き係数・マスク付きのNussinovを比較。
- 1〜3レベルの全探索で、DDの実行可能性、返した構造の目的値、上界を検証。
- 正負のPK接触、DPエネルギー、学習済みモデル、hybrid、学習項ゼロを検証。
- 全交差評価との比較、固定幅の制限、20,005塩基の非再帰トレースを検証。
- CLI、ネイティブBPPのroundtrip、閾値探索、refinement、不正な設定を検証。
- ILPなしの5テスト、HiGHS付きの63テスト、ASan/UBSanを実行。

速度と構造の比較は `experiments/dual-decomposition/benchmark.py` を使用する。
元ブランチの実行ファイルとDDに同じ初期BPPを渡し、refinement=0、各3回を
交互の順序で測定する。比較対象は既存の3サンプルと16個の既存NMR参照RNA。
別途、固定予算と疎な長距離候補・局所Hモチーフで長さを倍増する測定を行う。
生物学的な同等性、未観測RNAでの一般化、常にILPより速いことは主張しない。
具体的な結果は同ディレクトリのreport.jsonとresults.mdに保存する。
