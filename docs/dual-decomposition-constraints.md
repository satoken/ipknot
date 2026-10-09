# 固定構造・NMR制約付き双対分解

`codex/dual-decomposition` の `1d88e79` から分岐した `codex/dd-constraints` で実装。
ベンチマーク中の元worktree、実行ファイル、実験結果は変更していない。

`--decoder dd` と固定構造 `-c`、`--base-pairs`、`--stack-constraint` を併用できる。
既定の `--dd-constraints linear --dd-dp beam --dd-noe-solver relaxed` は
候補生成・制約構築・DD・主問題回復の全段階に
幅と予算を設ける。固定予算と固定個数・長さのNMR観測の下で、疎な入力規模に対する
期待線形時間・線形記憶量を保つ。ILPソルバーなしでも動作する。

## 対応する制約

- 固定対、非対合 `x`、未指定 `.`、相手方向 `<` / `>`、必須対合 `|`。
- NMR塩基対数のexact / lower-bound、非標準塩基対、canonical-neighbor候補フィルタ。
- スタックの通常・compact表記、1塩基バルジ、none / fallback / all。
- 3種類のflush coaxial stacking、実際の第三ヘリックス、直接の子でない構造の排除。
- 観測間の塩基対共有禁止・許可、coaxial helix face容量1。
- softの個数不足・超過とスタック欠測のペナルティ。
- crossing / blocksの符号付きPKスコア。

閾値以下の明示的な固定対は必ず候補に入れる。採用したNMR witnessに必要な対も、
閾値以下なら負の単項係数で追加する。hard制約を緩めた構造は返さない。
全列範囲・整数性・全制約行・同一レベルのplanarityを検査してから返す。

スタック支援は従来ILPの上流・下流端点近傍の行で表す。採用したバルジwitnessが
ある場合は距離2、ない場合は距離1。制約なしDDの直接隣接対を要求するDP条件とは
異なるため、制約付き・制約なし経路の予測が同じになるとは限らない。

## 線形性を保つ構築

全塩基対・全交差・全NMR instanceを列挙する旧経路を通らない。

1. 非標準対は通常候補の隣接座標から追加する。canonical-neighborを外す場合は、
   塩基の位置索引を使い、左端・塩基対型ごとに最大 `--dd-nmr-pair-beam` 個を
   残りの右端範囲から等間隔に採る。未採用の全右端を走査しない。
2. スタックは通常候補と上記のseedから延長する。各深さの候補を
   `--dd-nmr-pattern-beam`、観測ごとの最終witnessを `--dd-nmr-witnesses` に制限。
   fallbackは相補的な右鎖wordの索引で配列全体のdirect instanceの有無を調べる。
   ビームがdirect候補を失っても、bulgeを誤って許可しない。
3. coaxialの端末対は左端・右端の固定幅索引から照会する。まず固定幅の端末contextを
   残し、そのcontextだけについて通常候補を走査して実際の第三ヘリックスを探す。
   topology blockerは採用したwitnessごとに候補を一度走査して構築する。
4. 上下レベルの交差支援は既存DDと同じactive beam・行ごとのwitness制限を使う。
   同一レベルの対同士の交差禁止行は作らず、Nussinovと回復時のplanarity検査で守る。
5. 回復の1状態における区間伝播を `--dd-constraint-passes` 回の全行走査に制限する。
   強制された対の内側の親を左右の端点に付け、異なる親を持つ自由端点の候補を除く。
   これは `O(K*n+候補数)` の走査で交差を除く。対同士の全比較はしない。
   DFSは分岐の決定だけを保存し、各状態で領域を再構築する。
   状態ごとにモデル全体のコピーを保持するメモリ増大を避ける。

ハッシュ索引で行の重複係数を統合し、長いcount行などの全項のソートも避ける。
幅内の順位付け・ソートのコストは固定幅の定数に含まれる。

| オプション | 既定 | 制限する量 |
|---|---:|---|
| `--dd-dp` | beam | beam探索または枝刈りなしのnussinov DP |
| `--dd-noe-solver` | relaxed | NOE内部も緩和する従来方式、またはNOE専用ILP |
| `--dd-beam` | 100 | beam DPの区間開始点。nussinov指定時は参照しない |
| `--dd-crossing-beam` | 100 | 各レベルの未閉鎖候補 |
| `--dd-witnesses` | 16 | 上位対・下位レベルごとの交差相手 |
| `--dd-nmr-pair-beam` | 64 | 左端・型ごとの追加候補、coaxial端点索引幅 |
| `--dd-nmr-witnesses` | 64 | 観測ごとのwitnessとcoaxial context |
| `--dd-nmr-pattern-beam` | 32 | seedの各深さで保持するpattern延長 |
| `--dd-max-iter` | 50 | 各閾値・refinementのDD反復上限 |
| `--dd-constraint-states` | 2048 | 初期・途中・終了時を合計した回復状態数 |
| `--dd-constraint-recovery-every` | 0 | 途中回復を再開するDD反復間隔。0は終了時のみ |
| `--dd-constraint-passes` | 8 | 1回復状態での伝播走査回数 |

### 計算量の条件と内訳

配列長 `n`、入力の物理的な疎BPP対数 `M`、レベル数 `K`、NMRスタック観測数 `Q`、
最大pattern長 `H`、要求する非標準対型数 `U`、追加候補幅 `B`、NMR witness幅 `W` を
用いる。レベル込みの対候補数 `A` は次のように制限される。

```text
A = O(K * (M + n*(1+U*B) + Q*W*H))
```

通常候補からの隣接対追加は `O(M)`、固定対は最大 `n/2` 個。
NMRから注入する対は観測ごとに最大 `W*H` 個。coaxialの第三ヘリックスは通常候補。
交差witness幅を `R` とすると、制約の非ゼロ項数 `E` は

```text
E = O(n*K + K*R*A + Q*W*K*A + Q*W*H*K)
```

と抑えられる。項数にはendpoint、stack、count、PK積、NMRリンク・blockerを含む。
pattern生成は固定した深さ・幅の分岐、coaxial/blockerは固定個数のwitnessと候補の
走査なので、`K,Q,H,U` と全ビーム幅を固定すると構築は期待 `O(n+M)` になる。

反復数を `T`、回復状態数を `S`、伝播回数を `P` とすると、固定DP幅で

```text
時間 = O(T*(n*K+A+E) + S*P*(n*K+A+E) + S^2)
記憶 = O(n*K + A + E + S)
```

となる。`S^2` は各分岐状態の決定パスの復元に由来する。
**`K,Q,H`、全幅、`T,S,P` を配列長と独立に固定した場合、全制約付きdecoderは
疎な入力規模 `n+M` に対して期待線形時間・線形記憶量である。**
ハッシュ検索・selectionの平均計算量を用いる。NMR観測数やpattern長も増やす場合は
上の依存関係を使い、無条件に全入力サイズに対して線形とは主張しない。
密なBPPは `M=Theta(n^2)`。配列長に対して線形にするには、固定幅LinearPartitionなどの
疎な確率生成器、固定個数の閾値候補・refinement回数と組み合わせる必要がある。
ViennaRNA/CONTRAfoldや密な補助BPP全体への線形保証ではない。

## beam探索を使わないNussinov DP

`--dd-dp nussinov` は各レベルの部分問題を、すべての区間開始点を保持する
Nussinov区間DPで厳密に解く。beamのselection・枝刈りを行わず、`--dd-beam`を参照しない。
`--dd-dp beam`は従来のbeam版。互換設定の`--dd-beam 0`でも枝刈りを外せる。
DPの選択は塩基対確率を計算するエンジンやNMR候補グラフの選択と独立である。

レベル込み候補数をAとすると、疎な再帰式の主要な結合走査はO(n*A)。
密な候補での最悪計算量はO(n^3)、DPとそのtraceはO(n^2+A)記憶である。
traceは各保持セルの最良遷移だけを保存し、比較途中の全改善候補の記録による
立方オーダーのメモリ増大を避ける。MXfold2などの計算コストとは別に評価する。

`--dd-constraints linear`もこのDPと併用できる。この場合linearは候補・制約構築と
回復探索の方式を指し、DPの時間・記憶量に線形保証はない。
NMR説明候補や交差支援の制限は残るため、DPの厳密性と全整数問題の厳密性は別である。

## NOE専用ILPを組み込む場合

`--dd-noe-solver ilp` で、各DD反復のNOE副問題をILPで解く。
リンクされたILPソルバーが必要。既定値 `relaxed` は従来のソルバー不要DDを保つ。
`linear` と `full` の両方、およびbeam DPとNussinov DPの両方に対応する。
NOE観測がない場合、専用副問題は作らない。
候補構築は同じだが、ILPの探索時間は線形とは限らない。

塩基対は引き続き各レベルのDPで選ぶ。NOE用ILPの変数は説明witness、
coaxial context、soft違反だけで、塩基対変数やそのコピーを含めない。
観測ごとの説明選択、face容量などのNOE内部の行をILPに残し、RNAとNOEを結ぶ行を緩和する。
共有禁止については、説明の使用量が物理的な対の選択量以下である行と、
元モデルの対容量からNOEだけの容量行を導く。
coaxial blockerとrequired-pairリンクが同じRNA式を参照する場合も、NOE同士の競合行を導く。
別レベルの部分和を同じ式として扱わない。いずれも元モデルの実行可能解を除外しない行である。

例えば `y_w <= sum_level x_ij,level` を緩和すると、最大化目的に
`lambda * (sum_level x_ij,level - y_w)` を加える。
DP側には塩基対への加点、NOE側には説明候補への減点として同じ乗数が入る。
NOEをまとめて最適化するため、説明候補同士の整数不整合を副問題内で解消できる。
DPとの不一致やレベル間の整数ギャップが全て解消する保証はない。

HiGHSでは副問題の行列を一度だけ登録し、各反復で目的係数を変更する。
同じ係数が続く場合は解と上界を再利用する。他の対応バックエンドは各回モデルを再構築する。
MIPの相対・絶対ギャップを0に設定し、最適性を確認した解だけを使う。
HiGHSではソルバーの最大化双対上界をDDの上界に加える。
非最適終了を実行可能目的値で代用せず、エラーにする。上界の対象は採用した候補モデルである。

実行可能解の回復では、DPが選んだRNA構造を固定し、元のNOE罰則を最適化する
別のNOE割当ILPを解く。RNAだけの行に違反する候補は先に除外する。
回復解は下界と最良解の更新に使い、DDの上界や劣勾配には使わない。
同じ固定構造が続く場合、この割当結果も再利用する。

通常のRNA定式化では回復DFSからNOE変数の分岐を外し、RNA側の整数決定について探索する。
途中の構造候補や葉でNOE割当をILPへ任せ、返った整数解を全ての元の行とplanarityで検証する。
`--dd-constraint-states` はこのDFSの訪問状態数であり、NOE ILP内の探索ノードは含まない。
初期・途中・最後の回復で同じ予算と最良解を共有する。
NOEと連続スコアを直接結ぶ汎用モデルでは、連続区間の最適化を守るため元の整数分岐を保つ。

トレースのproblem/summaryにNOE変数数、行数、呼出し数、キャッシュ数、時間を記録する。
`--dd-trace-state` では `noe_columns` と `noe_rows` により、残したNOE内部の制約も検査できる。

```sh
ipknot --decoder dd --dd-dp nussinov --dd-noe-solver ilp \
  --dd-max-iter 50 -e MXfold2 -t auto,auto -r 0 \
  --stack-constraint "GC AU" input.fa
```

L2の8配列で50/100反復を確認した結果は
[統合後の評価](benchmarks/dd-noe-ilp-integration-20261006/report.txt)に保存した。
一部の配列では改善するが、一律の精度向上・高速化はまだ得られていない。

## 分解、上界と近似の意味

`IP(IPModel&)` はソルバーを構築せずに式を記録する。列を `v`、不等式行を
`A_r*v <= b_r` とし、最大化のラグランジアンを次のように置く。

```text
L(v,q) = c*v + sum_r q_r*(b_r - A_r*v)
```

不等式の `q_r >= 0`、等式の `q_r` は符号自由。レベルごとのNussinovを
`c-A^T*q` で解く。既定では補助列は有限範囲の最良端点を選び、
`--dd-noe-solver ilp` ではNOE列を内部の制約の下でまとめて最適化する。
劣勾配は `b_r-A_r*v`、更新は `q <- q-eta*g`、不等式だけ非負に射影する。
stack条件も緩和行なので、部分問題のNussinov自体は孤立対を許す。

符号付きPK接触はbinary積 `z=x_upper*x_lower` と3本のAND行で表す。
正負どちらの係数でも元の接触スコアを正しく計算する。
DPがビームで枝刈りされた場合は各端点の正の係数最大値の和を2で割る上界を使う。
枝刈りがなければDPの最適部分問題値を使う。
`graph_upper_bound` は**保持した候補・交差辺・NMR witnessからなるモデル**の
目的値に対する上界。枝刈りされたoracle値を上界とみなさない。
浮動小数点計算なので、形式的な区間演算証明ではない。

候補/witnessの制限は実行可能な構造を落とすことがあり、PKスコアも保持した接触だけで
計算する。したがって、この上界は無制限のILPモデルに対する上界ではない。
整数最適性・hard制約の実行可能解の発見を一般には保証しない。
返した解については全hard制約を検査する。soft違反は設定したペナルティで記録する。

回復は元の目的係数で枝を評価する。初期は必要な場合だけ最大半分の状態予算を使い、
最初の実行可能解を探す。soft観測には空構造と違反変数も初期候補として検査する。
各状態では強制済みの対と未確定列の下限による完成解も検査し、不要な0への全分岐を避ける。

`--dd-constraint-recovery-every 10`ではDDの10反復ごとに、同じDFSの未探索部分から回復を再開する。
探索の分岐決定と未訪問のfrontierを保存し、最良の実行可能解を全期間で引き継ぐ。
既に選んだ分岐は保持し、新たな分岐にはその時点までのDD調整重みを積算した優先値を使う。
改善した最良解は下界・停止判定・最終出力に反映する。
既存の`--dd-recovery-target baseline`（既定）ではstep targetを初期解とDD部分問題から
直接得られた実行可能解に保ち、途中回復によるDDの軌道変更を抑える。
`--dd-recovery-target best`では途中回復の下界も直ちにstep targetへ反映する。
この選択は目的関数・F1・収束速度のいずれにも一律の優劣を保証しない。

`--dd-constraint-states S`は1回の回復ごとの予算ではなく、1閾値・1refinementの
初期・途中・終了時の訪問状態数の合計上限。各状態はDFSで取り出した部分的な変数割当てで、
制約伝播、完成解検査、目的値上界による枝刈り、必要なら次の整数分岐を行う。
塩基対数・候補数・DD反復数とは異なる。auto,autoの10閾値なら各閾値に独立なSがあり、
最大で10*S状態となる。制約伝播の走査回数は別の`--dd-constraint-passes`で制限する。

初期探索後の残り予算の約5%（整数切り捨て・最低1状態）ずつを最初の4checkpointで使う。
間隔10なら10,20,30,40反復で約20%を訪問し、50反復で残りを使う。
最後の区分は、後半のDD調整重みに基づく別の分岐順でrootから探索する。
最良整数解と累積訪問数は保持し、早い分岐順だけに残り予算を固定しない。
50より短く終了した場合は残りを終了時に使う。総予算が小さい場合は途中で使い切ることもある。
この配分はdd-max-iterと独立なので、同じ入力・他の設定・patience0の50/100/500反復では
50反復までの回復履歴を共有し、長い実行が50反復時の最良整数解を失わない。
総予算使用後にも、DD部分問題から得られた実行可能解の検査と上界の改善を続ける。
この性質は固定閾値の目的値についてであり、F1やautoで選ぶ最終構造の単調改善は保証しない。

`--dd-constraint-recovery-every 0`（既定）は途中回復を無効にする。初期と終了時の探索も
未探索部分を引き継ぐ。`--dd-constraints full --dd-constraint-states 0`は無制限探索で、
最初の回復機会に探索完了まで走る診断設定である。

- 実行可能解があれば予算切れでも保存解と上界を返す。
- 実行可能解がなく予算切れなら、探索不足を明記してエラーにする。
- 伝播・全探索で不可能と判明した場合も、**保持したgraphの不可能性**と明記する。
  幅を増やすかfullモードで診断できる。無制限モデルの不可能性と混同しない。

linearモードで交差の幅または回復状態予算に0を指定するとエラーにする。
交差グラフ・整数回復の無制限設定は `--dd-constraints full` を明示する。
DPだけの無制限設定は両モードで`--dd-dp nussinov`を使える。

## full診断モードと共通の挙動

`--dd-constraints full` は従来ILPと共有する候補・全交差・NMR制約行を記録してDDで解く。
NMR instance全列挙、全行、収束までの伝播を保持するため、線形保証はない。
`--dd-beam 0 --dd-constraint-states 0` は小規模の厳密比較向け。
PK積は既存の独立continuous score列で表し、整数列確定後にその区間を最適化する。

共通のstep、schedule、patience、projected-norm、traceを使える。
制約なし用の追加bound / joint / exchange / recoveryは制約付き経路で未対応。
回復には `--dd-constraint-states` と `--dd-constraint-recovery-every` を使い、targetはbaseline/bestを選べる。
自動閾値探索では、保持したモデルで証明された不可能な閾値組を飛ばす。
探索予算切れで閾値を飛ばさない。すべて不可能ならエラーになる。
FASTA・alignmentのrefinementでは元の固定制約を各回で再適用する。
auxiliary BPP入力にも固定構造・スタック観測を適用する。
alignment入力へのNMR観測の適用は既存CLIと同様に対象外で、固定構造だけを扱う。

traceの制約付きイベントは `constrained: true`。`--dd-trace-state` のproblemに
列係数・領域・対と列の対応・緩和した全行を保存する。
iterationのselectedは部分問題解、recoveredは保存された実行可能解。
recoveryイベントに再開時の反復、periodic/final、訪問数、累積状態数、下界、探索完了、分岐順のrestartを保存する。
最後のrepairに最終列値、summaryに目的値・上界・処理量カウンタを保存する。
summaryのrepair_calls/periodic_repair_callsは探索再開回数、repair_budget_exhaustedは
総状態予算に到達し未探索部分が残ったことを表す。dp_pruned_states/exact_dpでDP枝刈りを確認できる。
主問題解が未発見ならlower_boundはJSONのnull。
`nonzeros` は重複統合後の項数、`propagation_work` は伝播の項走査量、
`structural_work` はplanarity検査の配列と候補の処理量。
ログにはNMR seed数、幅で落とした候補・witness数も出す。

## 実行例

```sh
cmake -S . -B build-dd -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build-dd -j2
ctest --test-dir build-dd --output-on-failure

build-dd/ipknot --decoder dd -r 0 -t .01,.01 \
  --base-pairs UU=1 --without-canonical-neighbor tests/data/nmr_constraints.fa
build-dd/ipknot --decoder dd -r 0 -t .01,.01 \
  --stack-constraint 'GC AU' tests/data/nmr_single_bulge.fa
build-dd/ipknot --decoder dd -r 0 -t 0,0 --coaxial-stacking \
  --stack-constraint 'UU CC' --stack-constraint 'GC AU' \
  -c tests/data/nmr_coaxial_bulge_children.constraint \
  tests/data/nmr_coaxial_bulge_children.fa
```

## 検証と長さ倍増の測定

`dd_constrained_exhaustive` は符号付き行・重複係数・整数slackの100モデルと、
交差・共有端点を含む任意の対の120モデルを独立全列挙と比較する。
後者には交差禁止行を入れず、暗黙のplanarity判定を検証する。
1024項のcount行が固定回数の走査で処理されること、負の固定対、signed continuous PK積、
beam幅0/1/100の上界、探索不足と不可能性の区別も確認する。
`dd_constraints_cli` は固定構造、個数、soft、stack/bulge、3種類のcoaxial topology、
pair/face容量、aux、自動閾値・refinement、full診断モードを検証する。
fallbackのword索引はランダム40入力で独立した全配列走査と比較する。

`dd_constraint_scaling` は疎な長距離入れ子固定対＋count/stack、全U配列の非標準対＋soft、
複数multibranchの固定coaxialの3系列を使う。
候補・行・項数の密度と長さ倍増時の増加率、伝播回数の予算による上限を検査する。
返したBPSEQは入力の固定対または解析的に最適な空構造と独立に比較する。
時間の比率をCIの合否条件にしない。

2026-10-05、aarch64 Release、各点3回の中央値。
DP/交差/NMR幅8、反復3、回復状態32、伝播4、2レベルを全長で固定した。
[全設定・処理量・測定値](dd-constraint-scaling.json)を保存した。

| 系列 | 長さ | 非ゼロ項数 | 秒 |
|---|---:|---:|---:|
| 入れ子固定対＋count/stack | 2048 | 21,500 | 0.0099 |
| 同上 | 4096 | 43,004 | 0.0290 |
| 同上 | 8192 | 86,012 | 0.0334 |
| 非標準対＋soft | 2048 | 519,792 | 0.1949 |
| 同上 | 4096 | 1,040,817 | 0.4667 |
| 同上 | 8192 | 2,084,549 | 0.9816 |
| 固定coaxial | 2052 | 11,788 | 0.0090 |
| 同上 | 4104 | 23,452 | 0.0161 |
| 同上 | 8208 | 46,780 | 0.0306 |

```sh
python3 tests/dd_constraint_scaling_test.py build-dd/ipknot \
  --lengths 512,1024,2048,4096,8192 --repeat 3 --report scaling.json
```

ソルバーなし12テスト、HiGHS付き70テストが成功。
既存NMR回帰50ケースを既定linearモードでも実行し、成功・失敗とログ条件がすべて一致。
成功36ケースではHiGHSとの目的値も1e-6以内で一致した。
この一致は小規模fixtureの検証であり、一般入力で同じ予測・整数最適性を保証しない。
