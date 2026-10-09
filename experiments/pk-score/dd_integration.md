# DD改良beam100へのPKスコア統合

`codex/dual-decomposition` の `4b93a1e263c1420802ebf05cb7be8273ddef10ed` を、PK実験用の `codex/pk-score` にfast-forwardで統合した。元の `dev` とDD側のworktreeの作業ファイルは変更していない。

**指定された改良beam100・LPC100・最大50反復で、crossing型の形状、BPP、統合、DP、CCスコアは動作した。** 統合スコアの600配列の時間はILPの111.597秒に対してDDは7.958秒、約14倍速い。DDにスコアを加える時間増は約8.4%。ただしBPP計算を含まない比較であり、統合スコアによるPK対F1の上昇は両標本とも95%区間に0を含む。

## 指定された改良beamの再現

保存された [設定](../../docs/presentations/dd-overview-20261006/data/improved100-manifest.json) では、「DD（改良beam100）-LPC」はDDブランチ本体と異なる実験バイナリを使っていた。このバイナリの、beamを適用する前に支配される区間状態を除く処理と決定的なtracebackの同点処理を、`--dd-dp improved-beam` として取り込んだ。通常の `beam`、厳密な `nussinov` とは区別できる。実験用のheapやbeam幅の自動拡大は有効にしていない。

```
--decoder dd --dd-dp improved-beam --dd-beam 100
-e lpc --beam-size 100 -t auto,auto
--dd-max-iter 50 --dd-patience 0
--dd-crossing-beam 100 --dd-witnesses 16
```

最大50反復は、各閾値・各refinementの求解ごとに適用される。主比較は以前の依頼に合わせて `-r 0`。元のDD比較と同じ `-r 1` でも、別途、予測の一致とスコアの動作を確認した。

元の実験バイナリは `/tmp/ipknot-dd-heap-study/heap-build/ipknot`、SHA256は `698d762a25e570f4f3b13bb21a2405a65bde8b4d4622a338365af5a1d39f11ba`。保存された `src/dual_decomposition.cpp` は `1f185da41e6950284bb0bdbaaba87a3ea4b34f21114e34933bbaeb17bcc0e62f`。通常beamの既定値は変えていない。ただし、取り込んだDDブランチ自体がデコーダの既定値をDDに変更しているため、ILPを使う場合は `--decoder ilp` を明示する。

## スコアと比較条件

形状/BPP統合は [以前のモデル](models/hybrid-r0-v1-posterior.txt) と重みをそのまま使う。固定最大候補ステムブロックA、Bの長さと三つの区間長Lrについて、接点係数は

```
s = -0.05 + 0.025 min(A,B) - 0.0075 Σ log(1+Lr)
w_AB = s/(A B) + 0.05 (θ・f)/anchor_length
```

形状のみ、BPPのみ、統合を比較し、DPは以前のR0で選ばれた非ゼロ設定 `intercept=.05, scale=.002`、CCは以前の非ゼロ診断条件 `intercept=.1, scale=.005` とCC09表を使用する。以前のCCの選択結果は重み0だったため、このCC条件を推奨値とは扱わない。いずれも形状範囲はステム・ループ上限200、閾値選択へのスコア重みは0。再学習・DD向けの重み再調整は行っていない。

DDは保持された交差接点について、選択された上位対と下位対の両方が成立するときだけ `w_AB` を目的関数へ加える。負の係数も同時成立時に課金する。形状・DP・CC・BPP・統合それぞれで、保持された接点の係数がILP側の計算と一致することをテストした。

H型ブロック対の係数はキャッシュし、全ブロック対を総当たりしない。固定された交差beam・支持相手数・最大ステム長・反復数のもとでは、疎な候補数Mと配列長nに対して期待O(n+M)というDD側の計算量を保つ。スコア付き接点の局所最適化には定数倍の処理が加わる。CC09表のファイル読み込みは配列長とは別の固定費であり、CLIを配列ごとに起動する今回の時間にも含む。

DDの交差beam100と各支持行16相手という制限は、ILPの全交差接点を必ず保持する設定ではない。相手の選別にはBPP目的係数とPK係数が使われる。したがって、同じ重みでも保持される接点と実行可能集合が異なることがある。自動閾値選択のPK期待F値もDD側は有界な交差サンプルを使う。精度差を単にbeam DPの誤差だけに帰すことはできない。

この検証の基準コミット `2f08323` では `projected`, `supported`, `exact`, `rerank` はDD未対応だった。このページの測定は `crossing/blocks` を対象とする。その後追加した `projected/blocks` のDD対応と比較は [dd_projected.md](dd_projected.md) に記録する。`supported`, `exact`, `rerank` は引き続き `--decoder ilp` が必要になる。

## 精度・時間の比較

以前に評価した二つの300配列（それぞれPK陽性150・陰性150、12–200nt）を再使用する。両方とも今回の新しい独立検証データではない。全条件に同じ初期LPC100のBPPを与え、refinementを切る。時間はBPP計算を含まない壁時計時間。条件順をRNA・反復ごとに回し、逐次3反復のRNAごとの中央値を合計する。F1はmicro平均、PK対は少なくとも一つの対と交差する塩基対として評価する。

測定環境はaarch64、Release、HiGHS、ViennaRNA 2.7.2。CPU affinityは固定せず、許可された0–19番を使用した。各配列につきCLIを起動しており、起動・BPP読み込み・予測出力、CCでは表の読み込みも時間に含む。時間比は今回の同時期・同条件の比較から求め、以前の別測定の秒数とは直接比較しない。FASTAからのBPP計算込みや、一度に複数配列を処理する使い方で同じ時間比になるとは限らない。

測定値と対応比較の95%区間は [summary](measurements/dd-integration-summary.json)、配列ごとの値は [per RNA](measurements/dd-integration-per-rna.json)。重みを再調整しない回帰検証と探索的な精度比較であり、RNAファミリーを独立に分けた精度保証ではない。

以前の統合スコア検証に使った先の300配列。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK | 時間 |
|---|---:|---:|---:|---:|
| ILP、スコアなし | 0.624928 | 0.231560 | 64 | 35.720秒 |
| ILP、統合 | 0.628189 | 0.243064 | 62 | 55.691秒 |
| DD、スコアなし | 0.624238 | 0.228545 | 65 | 3.540秒 |
| DD、形状 | 0.624919 | 0.231174 | 66 | 3.830秒 |
| DD、BPP | 0.627369 | 0.233772 | 66 | 3.836秒 |
| DD、統合 | 0.627324 | 0.236665 | 65 | 3.844秒 |
| DD、DP | 0.625584 | 0.231579 | 64 | 3.828秒 |
| DD、CC | 0.624955 | 0.231190 | 65 | 16.363秒 |

以前のprojected比較で追加した後の300配列。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK | 時間 |
|---|---:|---:|---:|---:|
| ILP、スコアなし | 0.643634 | 0.223608 | 82 | 35.346秒 |
| ILP、統合 | 0.644409 | 0.228205 | 81 | 55.906秒 |
| DD、スコアなし | 0.645019 | 0.232143 | 81 | 3.802秒 |
| DD、形状 | 0.644011 | 0.233931 | 79 | 4.113秒 |
| DD、BPP | 0.646407 | 0.239565 | 79 | 4.118秒 |
| DD、統合 | 0.644453 | 0.234266 | 79 | 4.114秒 |
| DD、DP | 0.643492 | 0.232972 | 80 | 4.095秒 |
| DD、CC | 0.648970 | 0.232373 | 80 | 16.668秒 |

統合−DD基準のPK対F1差は先の300本で+0.008120、95%区間 `[-0.002616,+0.022742]`、後の300本で+0.002123、区間 `[-0.008454,+0.009983]`。全対F1は先で+0.003087、後で−0.000565だった。統合−ILP統合のPK対F1差は先で−0.006398、区間 `[-0.012229,-0.001284]`、後で+0.006061、区間 `[+0.000607,+0.012160]`。DDとILPの精度の優劣は標本によって変わった。

BPPのみは、後の300本でDD基準からのPK対F1差+0.007422、区間 `[+0.002086,+0.013832]` だった。一方、先の300本では統合のPK対F1がBPPのみより高い。これは既知の評価データで複数条件を比較した結果であり、DD向けの新しい重み選択や独立検証は行っていない。

DP/CCによるPK対F1の上昇は小さい。CCの時間増は毎回47,451行のCC09表を読み込む固定費を含む。CLI内では表を入力配列ループの前に一度読み込むため、複数FASTA配列を一つのプロセスに与える場合にはこの固定費を共有できる。

## 動作と回帰の確認

- 重み0のDDは、全600配列で基準DDと予測・各閾値の目的値・反復数・停止理由が完全一致。
- マージ後のILP基準・統合は、600配列×3反復すべてで以前の予測と完全一致。[回帰記録](measurements/dd-integration-ilp-regression.json)。
- 保存された改良beamバイナリと、24配列×統合/DP/CCの72比較で予測・目的値・PKスコア・反復数・停止理由が一致。
- R0の72条件とR1の8条件、合計854グラフ・15,183反復状態で、塩基の排他性、レベル内非交差、スタック支持、各下位レベルの交差支持を独立に検査。選択された接点の符号付きスコアから目的値を再計算し、最大誤差は `1.50e-13`。
- 小規模な符号付きDPの全列挙との一致、1024状態を超えるradix処理、ILPとDDの保持接点係数の一致、制約付き経路とCLIを含む79件のCTestがすべて成功。

長鎖は、以前の速度診断に使った562・1000・2883ntの3配列でLPC100、R1、auto/auto、改良DD100、最大50反復を確認した。基準DDの予測は保存された結果と一致し、統合スコアは新旧の改良beamバイナリで予測と各グラフの数値が一致した。562・2883ntでは非ゼロの接点スコアと実現PKスコアがあり、1000ntでは今回のH型スコアで採点される保持接点が0で、PKスコア0・基準DDとの完全一致を確認した。各配列に必ず非ゼロのH型スコアが入るとは限らない。これは1回ずつの動作確認であり、長鎖の一般的な精度・時間の評価ではない。[長鎖記録](measurements/dd-integration-long.json)。

## 再現方法

統合スコアを指定されたDD設定で使う例。`-r 1` にすればrefinementありになる。

```sh
build/ipknot --decoder dd --dd-dp improved-beam --dd-beam 100 \
  --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  -e lpc --beam-size 100 -r 0 -t auto,auto \
  --pk-h-formulation crossing --pk-h-allocation blocks \
  --pk-h-max-stem 200 --pk-h-max-loop 200 --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/hybrid-r0-v1-posterior.txt \
  --pk-learned-scale 0.05 --pk-hybrid-shape \
  --pk-h-intercept -0.05 --pk-h-stem-reward 0.025 \
  --pk-h-loop-penalty 0.0075 --pk-h-weight 1 input.fa
```

HiGHSをリンクしたReleaseビルドで実行する。以前の初期BPPキャッシュとモデル・CC表が必要。

```sh
python3 experiments/pk-score/dd_integration.py run
python3 experiments/pk-score/dd_integration.py audit
python3 experiments/pk-score/dd_integration.py summarize
python3 experiments/pk-score/dd_integration_regression.py
python3 experiments/pk-score/dd_integration_long.py
ctest --test-dir build --output-on-failure
```

生のコマンド、ログ、BPP・バイナリ・ソースのSHA256、失敗履歴は `results/dd-integration/` に保存する。失敗を除外したり、成功するまで再実行したりしない。長鎖のR1確認は [別記録](measurements/dd-integration-long.json) に保存する。
