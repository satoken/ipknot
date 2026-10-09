# シュードノットスコアの試作

このディレクトリは2026-10-09に `/tmp/ipknot-pk-score` から統合した実装の説明とモデルです。
大型の入力データ・計測結果・実験スクリプトは元worktreeに残しています。
以下の実験記録は元チャットでの評価であり、今回の統合検証とは別です。
統合範囲と検証結果は [統合記録](../../docs/pk-score-merge-20261009.md) を参照してください。

DP付き相互排他補正crossingを、過去の評価RNAを除いたRfam残余全件で検証した結果は [large_exclusion_benchmark.md](large_exclusion_benchmark.md)。未使用4,445本（PK陽性99、陰性4,346）で、PK対F1の改善は確認できず、追加のbalanced比較ではbaselineより低かった。

新log crossingの非ゼロ診断を同じ残余4,445本で評価した結果は [large_log_crossing_benchmark.md](large_log_crossing_benchmark.md)。T=.025、.2ともPK対F1はbaselineより低く、偽PK判定の減少と引き換えに正解PK対も落ちたため採用しない。

局所Boltzmann期待値補正のcrossing/projectedを、過去4コホートと残余Rfam全件を合わせた5,645本で評価した結果は [large_local_boltzmann_benchmark.md](large_local_boltzmann_benchmark.md)。残余4,445本では両方式ともPK対F1が有意に低下し、偽PKも増えた。既定値にはしない。

短いステム区間、排他性補正、成立した最良相手、構造候補の順位学習をcrossing/projectedで比較した結果は [alternatives_evaluation.md](alternatives_evaluation.md)。新規300でPK対F1最大は.22064→.24242だが、偽PK66→79件で差の95%区間も0を含む。DP付き排他crossingは.23541・時間1.097倍。同じ倍率では最良相手方式単独の改善は見られず、既定値は変更しない。

DP補正倍率を線形κ=3.2・局所κ=32まで広げた探索は [dp_expanded_scale.md](dp_expanded_scale.md)。選択が変わったのは局所projectedのκ=1→2。さらに新しい300配列ではPK対F1が.22970→.23060だが、差の95%区間は0を含む。倍率を増やすだけで安定したPK精度向上は確認できなかった。

オリジナルDirks–Pierceの線形補正・局所Boltzmann期待値・DD候補集合の3方式を実装し、refinementなしで比較した結果は [dp_conversion_evaluation.md](dp_conversion_evaluation.md)。新規300では局所crossingのPK対F1が0.24494→0.24799だが差の95%区間は0を含む。候補集合方式はDD/BPPを追加せず約2%の時間増で動く一方、PK対F1が下がる。既定値は変更しない。数式の設計記録は [dp_score_conversion.md](dp_score_conversion.md)。

塩基配列を特徴に使わない4係数スコアの実装、refinementなしDDでの900配列の評価は [sequence_free_evaluation.md](sequence_free_evaluation.md)。新規300で選択projectedのPK対F1は0.24583→0.24739だが差の95%区間は0を含む。既存hybridの0.25011には届かず、既定値にはしない。元の設計案は [sequence_free_score_design.md](sequence_free_score_design.md)。

DDのprojected対応とcrossingとの600配列の比較は [dd_projected.md](dd_projected.md)。既存の交差支持行を保ち、符号付きの係数補正だけを加える。追加の積因子は作らない。

DDブランチの統合と、指定された改良beam100・LPC100・最大50反復でのPKスコア検証は [dd_integration.md](dd_integration.md)。crossing型の形状/BPP統合、DP、CCを既存の重みで比較する。

refinementを切り、同じ形状＋BPP統合スコアをcrossingとprojectedへ配分した精度・時間比較は [projected_hybrid.md](projected_hybrid.md)。固定最大ステムブロックのprojectedも、追加変数・追加制約なしで使用できる。

Dirks–PierceとCao–Chenの物理項を構造補正に利用する検討と、調整用データだけで行うCCループ表の追加実験は [energy_models.md](energy_models.md)。PK固有のエントロピー、単位変換、交差辺への配分と表の適用範囲を扱う。

DP・CC06・CC09を交差支持行の補助変数に組み込んだ実装、未使用300配列の精度評価、長鎖RNAでの共通BPPによる時間比較は [energy_implementation.md](energy_implementation.md)。候補ステム対の重みを固定して交差辺へ配分する。

その結果を受けた、対合相手との競合・交差相手の支持を使う少数特徴の学習スコア案と、調整用データの誤予測診断は [score_design.md](score_design.md)。一対あたりの補正を既存の交差辺へ配分する設計で、新スコアの精度向上は未検証。

既存の交差制約一行につき一つの連続補助変数を置き、選ばれた交差相手へスコアを付ける追加実験は [crossing.md](crossing.md)。相手ごとのAND変数を作らず、元の支持候補と実行可能集合を保つ。

既存のレベル間交差制約を利用して、projectedの加点に対応する相手を実際に選ばせる追加実験は [supported.md](supported.md)。変数・制約行数を増やさず同時成立条件を一部保つが、実行可能集合を絞る方式である。

計算時間を重視するなら、まず **既存の自動閾値探索の候補解を、成立したモチーフのスコアで選び直す `rerank`** が扱いやすい。整数計画法の目的関数へ直接入れる場合は、追加変数・制約を作らない `projected` が軽い比較対象になる。ただし後者は候補の幾何情報を使う近似で、実際の相手ステムの成立を要求しない。モチーフの同時成立を厳密に入れる `exact` は実装できたが、今回の測定では求解時間が大幅に増える例があった。

**実験の状態。** すべて任意指定で、デフォルトの特徴重みはゼロ。非ゼロの推奨パラメータはまだ決めていない。現在の作業中だったNMR変更を基準コミット `7d560bb` に複製し、その上で実装した。基準コミットの親は `61900da`。元の `dev` の作業ファイルは変更していない。

**スコアの定義。** H型の配列上の順序を `S1-left, L1, S2-left, L2, S1-right, L3, S2-right` とする。ステムは連続した塩基対で長さ2以上。長さを `a,b`、三つの区間長を `L1,L2,L3` とし、現在の特徴スコアは次式。

\[
\Phi(m)=\beta_0+\beta_s\min(a,b)
 -\beta_L\sum_{r=1}^3\log(1+L_r)
 +\beta_c\mathbf{1}[L_2=0],\qquad s(m)=\lambda\Phi(m).
\]

対応する引数は `--pk-h-intercept`, `--pk-h-stem-reward`, `--pk-h-loop-penalty`, `--pk-h-coax-bonus`, `--pk-h-weight`。`beta0=0.2, betaL=0.03` は計算量と動作を比較するために置いた未学習の実験値。熱力学パラメータではない。ステム項はBPPと安定性を二重評価する可能性があるため、初期の比較ではゼロにした。中央区間ゼロの項も、配列依存のcoaxial stackingエネルギーではなく配置の指標である。後続の [係数調整と独立評価](calibration.md) では、同じ特徴式のステム係数も調整対象に含めた。

`--pk-h-table FILE` の行 `a b L1 L2 L3 score` は、該当する形状の `Phi` を置き換える。最終係数には `--pk-h-weight` が掛かる。未掲載の形状は上の特徴式へ戻るため、特徴重みがすべてゼロなら未掲載形状へのスコアはゼロ。重複行、負の長さ、非有限の値、空の表はエラー。表の記載だけではスコア対象の長さ範囲は広がらない。

**スコアの選択肢。**

| 設計 | 得られる情報 | 検証で注意する点 |
|---|---|---|
| 開始コスト・区間長の少数特徴 | シュードノットの頻度とコンパクトさ | 短すぎるループの立体的な不適合はlog項だけでは表せない |
| ステム長×ループ長の表・交互作用特徴 | BPPにないステムとループの幾何的な相関 | 疎な形状では平滑化し、RNAファミリーを分けて検証する |
| Cao–Chen等から作るエントロピー・自由エネルギーの表 | 物理モデルに基づく形状の評価 | モデルの適用範囲と単位の校正、BPPとの二重評価に注意する |
| ループ–ステム三次接触の特徴 | PK固有の安定化 | 配列・配置の一致を別途定義し、候補を増やしすぎない |

Cao–Chenのモデルはループ長だけでなく、それがまたぐステム長への依存を扱う。[2006年の論文](https://doi.org/10.1093/nar/gkl346)と[interhelix loopを扱う2009年の論文](https://pubmed.ncbi.nlm.nih.gov/19237463/)が参考になる。この試作には、それらの数値表や三次相互作用モデルはまだ移植していない。表を使う場合は、例えば形状で決まるループ形成の自由エネルギー寄与に対して `Phi=-DeltaG_PK_shape` を適切な適用範囲で作り、混合係数 `lambda` を調整する。現在の表のキーは長さだけなので、末端塩基の種類や配列に依存する寄与は直接区別しない。BPPの目的関数はエネルギーではないので、その和を正しい事後分布やMFEと解釈しない。

学習するなら、まず形状の開始・区間長・ステムとの交互作用だけを小さなモデルにし、正解PKだけでなくPKなしの配列と誤った交差候補を負例にする。近縁RNAを学習と評価の両方へ入れない。候補保持の閾値、最適化の `lambda`、最終選択の重みは別々に調整する。形状表の参照コストは表の大きさに依存するが、同じ候補形状への係数だけを更新する場合は整数計画法のサイズを増やす必要がない。初期の16配列、PK陽性3配列では汎用の形状表を学習・検証するには足りず、重みの学習は実施しなかった。後続の [calibration.md](calibration.md) で、282配列による既存特徴の係数調整と、別ソースの600配列による評価を記録する。

**整数計画法への入れ方。**

| 方式 | 追加変数・制約 | 評価するもの |
|---|---|---|
| `--pk-level-penalty` | なし | 上位レベルの塩基対への一律罰則。分解に依存する比較用の事前分布 |
| `--pk-h-formulation exact` | 共有する連続補助変数とAND制約 | 実際に選択された最大ステム二本の組合せ |
| `--pk-h-formulation projected` | なし | 上位レベルの塩基対へ配分した、可能なH型のスコア |
| `--pk-h-formulation rerank` | なし。追加のILP求解もなし | 既存の閾値探索が出す候補構造の、成立したH型スコア |

`exact` ではレベルの和 `u_ij=sum_l x_ij^l` を利用する。ステムは全構成塩基対が選択され、両端の隣接塩基対が選択されない場合だけ出現する。したがって長いステムの部分候補を多重に数えない。共通のステムと区間は共有する。ANDは両方向の制約を入れるため、負のスコアでも出現変数をゼロにして課金を避けられない。

各AND変数は `[0,1]` の連続変数でよい。整数の塩基対割当が決まれば、入力リテラルがすべて0/1なので、AND制約から補助変数も0/1になる。整数変数は増えないが、LP緩和の形と塩基対選択間の結び付きは変わるため、求解時間は増え得る。

`--pk-h-loop-mode unpaired`（既定）は三つのループ区間がすべて非対合の場合だけ成立する完全な単純H型。`span` は区間内の別の対合を許し、最大ステム二本のH型コアを評価する。`span` の長さは区間の幅であって、実際の非対合塩基数ではない。複雑なPKでは複数のステムの組を数えるので、「PK成分一個へのコスト」やループエントロピーとは扱わない。

`projected` は、各塩基対に対し、それを含む候補モチーフの `s(m)/当該ステム長` の最大値を取り、上位レベルの既存変数の目的係数へ加える。全候補の和にして候補密度で報酬が膨らむことを避けるための最大値である。既存のレベル制約から実際の交差は必要になるが、スコアを与えた相手ステムや区間が選ばれるとは限らない。レベル分解に依存し、`--no-levelwise` との併用はエラー。

`rerank` は各ILP解の最大ステムを取り出し、実際の構造へ同じ形状スコアを適用する。分解に依存しない。自動閾値探索の選択値を `pF+pF_pk-NMR_penalty+selection_weight*PK_score` とする。`--pk-selection-weight` はこの重みで、最適化の `lambda` とは別に調整できる。固定閾値では候補が一つだけなので `rerank` は予測を変えない。スコア込みの全構造の最適解を求める方式ではなく、既存の候補集合の中の選択である。refinementの選択が変わると次回のBPPと求解時間が変わり得るが、求解回数やBPP計算回数を増やさない。

**候補数と計算量。** 既存の疎な塩基対候補から、長さ上限 `B` 以下のすべての部分ステムを生成する。塩基対候補数を `E` とするとステム候補数は高々 `E(B-1)`。二本目のステムは左端と右端の局所範囲から検索し、全ステム候補同士の総当たりは行わない。`B` とループ上限 `L` を固定すれば、各候補から調べる座標・長さの組合せは配列長によらず有界になる。疎な索引の構築にはmapの対数時間がある。

`exact` の追加連続変数数は、使用するステム数・共有区間数・モチーフ数の和。制約行数はステムの長さと区間の長さに比例する項に、モチーフあたり高々6行を加えた数になる。非ゼロ要素数には各区間に入る既存候補の数も効く。`projected` はこれらの変数・行を作らない。`rerank` は候補の全列挙をせず、選択済みステムから局所範囲にある相手だけを調べる。

既定の対象範囲は `--pk-h-max-stem 12`, `--pk-h-max-loop 30`。これはスコアの対象範囲であり、それを超える構造を禁止する制約ではない。`exact` と `projected` には `--pk-h-max-motifs 10000` の資源上限があり、超えればエラーにする。途中の候補を黙って捨てない。対象外のPKにはこのH型スコアは入らないため、特に負のスコアは範囲を外れる構造への偏りに注意が必要。

`--pk-candidate-threshold` は目的関数の閾値と独立の候補保持閾値で、既定の `-1` は従来の候補カットを保つ。低い値を指定すれば負の単独係数を持つ候補も残せるが、候補数、既存の交差制約、求解時間が増える。BPPに載っていない塩基対は、この指定でも追加されない。

**測定条件。** ARM64、GCC 13.3、HiGHS 1.15.1、ViennaRNA 2.7.2、1スレッド、LinearPartition-C。手元のNMR参照16配列（47–111塩基、うちPK陽性3配列）を使用。各条件を順番に3回実行し、配列ごとの中央値を合計した。PK塩基対は「少なくとも一つの選択塩基対と交差する塩基対」と定義する。全塩基対F1、PK塩基対F1はmicro平均。小さな探索的実験であり、独立な評価データでの精度向上を示す結果ではない。

BPPを一度だけ計算して共通にし、固定閾値 `0.125,0.0625` で整数計画法を比較した結果は [pilot-summary.json](measurements/pilot-summary.json)。

| 条件 | 全塩基対F1 | PK塩基対F1 | 16配列の合計時間 |
|---|---:|---:|---:|
| 基準版 | 0.833 | 0.155 | 0.133 s |
| exact・非対合ループ | 0.844 | 0.163 | 3.196 s |
| exact・区間幅 | 0.836 | 0.152 | 2.181 s |
| projected | 0.849 | 0.162 | 0.142 s |
| 候補保持を0.025へ拡張、スコアなし | 0.822 | 0.172 | 0.725 s |
| 候補保持を0.025へ拡張、projected | 0.819 | 0.133 | 1.244 s |

この固定閾値条件では、候補を拡張しないスコア方式の正解PK塩基対数は8/41で変わらない。2MIYでは正解PK塩基対14個のうち6個がBPP出力にないため、目的関数だけでは回復できない。候補拡張は回復数を11個に増やしたが、誤検出も増えた。

通常のFASTA入力、自動閾値、refinement=1を含む比較は [rerank-summary.json](measurements/rerank-summary.json)。

| 条件 | 全塩基対F1 | PK塩基対F1 | 16配列の合計時間 |
|---|---:|---:|---:|
| 基準版 | 0.826 | 0.170 | 0.874 s |
| rerank・非対合ループ | 0.829 | 0.170 | 0.887 s |
| rerank・区間幅 | 0.828 | 0.165 | 0.925 s |
| projected | 0.847 | 0.167 | 0.947 s |

通常設定でも正解PK塩基対数は8/41で変わらず、projectedではPK誤検出が45個から47個へ増えた。今回の特徴・重みをデフォルトにする根拠にはしない。

長さによる性能比較は正解構造を用いず時間だけ測定した。[weight-sweep.json](measurements/weight-sweep.json) の739塩基の入力では、基準版0.988秒、projectedの重みを0.1倍にした条件0.761秒、0.25倍0.637秒、1倍3.196秒。係数の変更だけでも整数計画法の難しさは変わる。[scaling.json](measurements/scaling.json) では `exact` の両モードが同入力で10秒の外部タイムアウトになった。projectedの候補係数生成は約1ミリ秒以下。固定閾値比較の2N3Qではexactの定式化約0.7ミリ秒に対し求解約1.74秒で、主な増加は最適化にある。

訂正：初期の性能比較でRF00005と呼んだ739塩基の入力は、`examples/RF00005.fa` の10本のtRNAを連結した合成配列だった。元のFASTAは各72〜82塩基の別々の配列である。数値は合成入力での測定としてのみ解釈する。`benchmark.py` の読込を修正し、複数配列のFASTAを指定した場合は各配列を個別に扱うようにした。元の測定を再現するための入力は、次のように明示的に作れる。

```sh
python3 - <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, 'experiments/pk-score')
from benchmark import fasta_entries
sequence = ''.join(s for _, s in fasta_entries('examples/RF00005.fa'))
Path('/tmp/RF00005-concatenated.fa').write_text('>RF00005-concatenated\n' + sequence + '\n')
PY
```

実在する長鎖RNAでの追加測定は [calibration.md](calibration.md) に記録する。

**検証。** 既存50テストに加え、831通りの固定割当を独立な最大ステムのオラクルと照合した。正負のスコア、部分ステムの多重計数、ループ内の追加塩基対、空のループ、対象範囲、レベルの変更を含む。同じオラクルでrerankの構造スコアも検証。CLIでは候補カット、負の課金、表参照、ゼロ重み、自動選択、projected、rerank、無効引数を確認した。HiGHSで52テストが通り、GLPKでも同じ831割当のモデル検証が通った。他のバックエンドにもAPIを実装したが、この環境では実行検証していない。

**再実行。** このworktreeではビルド済みの試作が `build/ipknot`、変更前の実行ファイルが `experiments/pk-score/results/ipknot-baseline`。測定中は他のビルドやベンチマークを同時に走らせない。

```sh
cmake --build build -j 4
ctest --test-dir build --output-on-failure

# 共通BPP、固定閾値での比較
python3 experiments/pk-score/benchmark.py \
  --output experiments/pk-score/results/pilot \
  --variants baseline,off,level_penalty,h_compact,span_compact,projected_compact,expanded_off,expanded_projected

# 通常の自動閾値・refinementを含む比較
python3 experiments/pk-score/benchmark.py \
  --workflow default --threshold auto,auto \
  --output experiments/pk-score/results/rerank \
  --variants baseline,off,rerank_compact,rerank_span,projected_compact

# 成立した単純H型による候補解選択。以下の重みは未学習の実験値
build/ipknot --pk-h-formulation rerank \
  --pk-h-intercept 0.2 --pk-h-loop-penalty 0.03 \
  --loglevel info examples/drz_Ppac_1_1.fa

# 表のスコアを成立したH型へ適用。lambda=0ならHスコアを無効化
build/ipknot --pk-h-formulation rerank \
  --pk-h-table path/to/scores.tsv --pk-h-weight 0.1 \
  --loglevel info sequence.fa
```

測定の個別結果と予測塩基対は `measurements/*.json` に保存した。`results/` の実行ファイル、BPP、ログはGit管理の対象外。再現用の基準実行ファイルを別に作る場合は、基準コミット `7d560bb` をビルドして `--baseline` に指定する。
