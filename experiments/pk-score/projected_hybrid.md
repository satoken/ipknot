# 統合PKスコアのprojectedとcrossingの比較

**等倍projectedは、同じスコアのcrossingに近い精度で、約34%速かった。** 未使用300本ではPK対F1がcrossingの0.228205に対してprojectedは0.228037、差−0.000169の95%対応bootstrap区間は `[-0.004745,+0.006109]`。時間は31.364秒から20.620秒になり、基準版20.865秒とほぼ同じだった。精度差の区間は0を含み、小さな優劣は判定できない。

形状とBPP支持の統合スコアを、実際に選択された交差相手へ付けるcrossingと、最良の候補相手が全て成立すると仮定して上位対へ付けるprojectedで比較した。全条件でrefinementはオフ。以前のスコア・モデル・候補cut・固定最大ステムブロックを使い、学習し直していない。

## 同じスコアの二つの配分

[統合スコア](hybrid_r0.md)の接点係数をそのまま使う。候補ブロックの長さA、B、三つのループ長Lr、BPP特徴fと固定モデルθから、

```
s = -0.05 + 0.025 min(A,B) - 0.0075 Σ log(1+Lr)
w_AB = γ s/(A B) + λ (θ・f)/anchor_length
```

γ=1、λ=0.05が以前のcrossing主条件。高いBPP支持を持つ主導ステムをanchorと呼び、ILPの上位レベルとは区別する。モデルは [hybrid-r0-v1-posterior.txt](models/hybrid-r0-v1-posterior.txt)（SHA256 `f0cc53ede91372d35e632548f429a8d475fbf348f81415cfbcf963d0f6ca9578`）。

上位対uの選択変数をxu、下位対vをyvとすると、crossingは同時選択された接点にだけ加える。

```
Q_crossing = Σu Σlower<level(u) xu Σv∈lower w_uv yv
```

採点外の接点はw=0。今回のcrossing対照は、整数解の目的値を保つ `--pk-crossing-simplify --pk-crossing-hypograph` を使う。

projectedは、同じ候補H型ブロック対について「その下位レベルに存在する相手ブロックの全候補対が選ばれた時」の係数を計算する。上位対uと下位レベルごとに、最も高い候補ブロックの係数を残す。

```
a_u,lower = max over scored partner blocks B
            [w_AB × number of candidate pairs of B in lower]
Q_projected = Σu xu Σlower<level(u) a_u,lower
```

負の係数も符号付き最大値にする。全候補が負なら負の補正を残す。0係数・採点外の候補は候補形状の提案に含めず、その行に採点可能な相手がなければ補正0。これは以前の部分ステムprojectedと同じ扱いであり、実際に採点外の相手を選んで負の接点を避けられるcrossingとの差になる。

候補ブロックの全対が下位レベルに存在する場合、形状項は以前のprojectedの `γs/上位ステム長` と一致する。部分的な候補しかそのレベルに存在しない場合は存在する対数だけを使う。BPP項にも同じfull-partner投影を適用する。上位レベルの変数の目的係数を変えるだけで、変数と制約行は追加しない。

元の支持制約 `xu ≤ Σcrossing witnesses yv` は**各下位レベルについてそのまま残す**。そのため、選択された上位対には何らかの実際の交差相手が必ずある。ただし、スコアの根拠となった最良の相手ブロックが選ばれるとは限らず、選ばれた相手が一部だけでもfull-partner分の係数を使う。

例えば、候補相手の長さが4対で実際には1対だけ選ばれた場合、同じ候補ブロックからのcrossing補正は接点1個分、projectedは4個分になる。正の係数では加点を過大に、負の係数では減点を過大にする可能性がある。逆に、複数のブロックが同時に選べる場合のcrossingの合計を、最大値だけのprojectedが小さく見積もることもある。3レベル以上では下位レベルごとに最良相手を独立に仮定するため、その仮定が全て同時に実現可能とは限らない。

## 精度・時間の結果

未使用300本（PK陽性150本、陰性150本）での比較。時間はデコーダのみ、3回のRNAごとの中央値の合計。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK本数 | 時間 |
|---|---:|---:|---:|---:|
| 基準版 | 0.643634 | 0.223608 | 82 | 20.865秒 |
| crossing、同じγ=1／λ=0.05 | 0.644409 | 0.228205 | 81 | 31.364秒 |
| projected、等倍 | 0.643743 | 0.228037 | 82 | 20.620秒 |
| projected、調整用282本で選択した0.5倍 | 0.644714 | 0.226262 | 82 | 20.628秒 |

等倍projectedはcrossingより34.26%短時間（時間比0.65744、95%区間 `[0.58318,0.74897]`）、基準版に対して時間比0.98824（95%区間 `[0.96262,1.01193]`）。PK対のTP/FP/FNはcrossingが534/1828/1784、projectedが536/1847/1782。projectedは正解対が2対増えた一方で誤検出も19対増え、PK precisionは0.226080から0.224927、recallは0.230371から0.231234になった。最終予測が異なるRNAは14本。

主比較のPK対F1差 `projected−crossing` は−0.000169、95%対応区間 `[-0.004745,+0.006109]`。全対F1差は−0.000666、区間 `[-0.001690,+0.000214]`。基準版に対するprojectedのPK対F1差は+0.004429、区間 `[-0.001779,+0.011759]` で、この新しい標本単独では基準版からの改善も確定しない。

既知300本での対応比較は次の通り。この300本は以前の統合スコアの評価に使ったデータであり、独立検証とは区別する。

| 条件 | 全対F1 | PK対F1 | 陰性150本の偽PK本数 | 時間 |
|---|---:|---:|---:|---:|
| 基準版 | 0.624928 | 0.231560 | 64 | 20.175秒 |
| crossing | 0.628189 | 0.243064 | 62 | 30.222秒 |
| projected、等倍 | 0.627018 | 0.241096 | 62 | 19.967秒 |
| projected、0.5倍 | 0.625395 | 0.234052 | 62 | 19.906秒 |

既知300本の等倍projectedのPK対F1差は−0.001968（95%区間 `[-0.005854,+0.001526]`）。時間比0.66069（区間 `[0.56167,0.77601]`）。最終予測の変更は16本。

調整用282本では、PK対F1が基準0.186333、crossing0.192467、projected等倍0.189300、0.5倍0.189844だった。共通倍率の規則で選ばれたのは0.5倍だが、その後の二つの300本では等倍projectedよりPK対F1が低い。全体の倍率だけを調整する優位性は確認できなかった。主比較の等倍条件はこの結果にかかわらず固定している。

## なぜ速くなったか／相手はどれだけ崩れたか

新300本のcrossingと等倍projectedで、既存の求解回数はともに2973回。projectedは全ての閾値グラフで基準版と同じ総列数・総行数だった。LP反復数の合計はcrossing109761回、projected39265回、基準42395回。分枝ノード数はそれぞれ1369、1360、1362で、時間差は主にLP計算に表れている。スコア係数の構築はcrossing0.146秒、projected0.090秒で、壁時計時間の10.744秒差の大部分は求解時間にある。

projectedの最終構造には採点された上位対が612対あり、その最良候補相手が全て選ばれていたのは363対、部分的に選ばれたのは212対、一対も選ばれなかったのは37対だった。37対のうち17対には正の加点が入った。同点候補が複数ある場合は最も成立している相手を採用したので、同点の別相手だけが不成立であるケースは数えない。

最良相手が不成立でも、元の支持制約を満たす別の交差相手は選ばれている。したがって、これはPKとしての交差の消失ではなく、**採点した相手と実際の相手の対応が外れる**近似。中立な特徴出力と最終構造のレベルから、選択された閾値と投影目的値を独立に再構成し、native目的値との差は最大 `2.22e-16` だった。この近似は実際に起きているが、今回の精度差は小さい。構造ごとの誤りへの因果効果をこの集計だけで断定しない。

## 検証と保存データ

60 CTestが全て成功。固定割当の独立オラクルで、正負の係数、最良相手が不成立／部分成立／別レベルの相手、複数下位レベル、下位レベルに候補の一部しか存在しないブロック、BPPが欠ける場合の形状項の保持を検証した。CLIでは学習モデル＋projectedの既定ブロック割当、形状との統合、補正0の互換性、候補cutと不正なオプションを確認した。

調整2256回、既知評価4500回、新評価3600回、計10356回のスコア付き／基準デコーダ実行が全て成功し、コマンド・予測対数・ログとBPPのハッシュを監査した。projectedの候補と基準のILPサイズが一致するグラフは49649個。crossingは旧実行ファイル `45e60c1` と調整282本・既知300本の最終予測、閾値ごとの整数目的値（差1e−8以下）、採点対象が一致した。

- [集計と対応区間](measurements/projected-hybrid-summary.json)
- [RNA単位の精度・時間・solver計算量](measurements/projected-hybrid-per-rna.json)
- [未使用300本の参照構造・抽出履歴](measurements/projected-hybrid-fresh-references.json)
- [最良相手の成立状況](measurements/projected-hybrid-partner-diagnostic.json)

## 比較手順

- 調整用は以前のbpRNA282本（PK陽性141本／陰性141本）。学習係数θ、形状係数、二つの項の比率を固定し、projectedの共通倍率だけを0、0.125、0.25、0.5、1、2から選ぶ。全対F1が基準以上、PK陰性RNAでの偽PK本数が基準以下の条件でPK対F1を最大化する。倍率0も選択肢に含め、同点は0、次いで小さい倍率を優先する。
- 等倍projectedとcrossingが主な配分比較。倍率を調整したprojectedは追加の比較として扱う。基準版も同時に測る。
- 以前の300本は既知データの対応比較。新しい精度検証には、過去の全2082配列と一致するもの、およびglobal edit identityが80%以上の配列を除いた新しい300本を使う。PK陽性150本、長さの分布を合わせた陰性150本、長さ12–200塩基。
- 新しい配列も初期LinearPartition-C BPPを一回だけ計算し、9桁のfloat往復精度で保存する。FASTA入力とキャッシュ入力で基準予測・中立特徴出力が一致することを確認し、全スコア条件が同じキャッシュを読む。
- 全条件 `-r 0 -t auto,auto`、一スレッドHiGHS、元の候補cut、固定最大ブロック、maxstem/maxloop=200、閾値候補選択へのPKスコア直接加算0。追加のBPP計算やILP求解はしない。
- 評価は逐次3回、RNA／反復ごとに条件順を回転させる。時間はRNAごとの壁時計時間の中央値を合計する。モデル読み込みを含み、BPPキャッシュ作成・倍率選択・bootstrapは含めない。
- 旧実行ファイル `45e60c1` のcrossingも282本と既知300本で測り、予測・閾値ごとの整数目的値・候補ブロックの一致を確認する。projectedの候補ブロック／接点がcrossingと同じで、ILPの総列数／行数が基準版と同じであることもログから確認する。
- 対応identity-cluster bootstrapを5000回。新300本での「等倍projected−crossing」が主比較。調整倍率や基準版との差、既知300本の区間は探索的な比較として扱う。

短くPKを多く含むRNAの検証であり、identityでの分離はRNAファミリーの分離を保証しない。長鎖RNAや自然のPK頻度への一般化はこの実験から断定しない。

## 実行

```sh
cmake --build build -j4
ctest --test-dir build --output-on-failure
python3 experiments/pk-score/projected_hybrid.py tune
python3 experiments/pk-score/projected_hybrid.py freeze
python3 experiments/pk-score/projected_hybrid.py paired
python3 experiments/pk-score/projected_hybrid.py prepare
python3 experiments/pk-score/projected_hybrid.py cache
python3 experiments/pk-score/projected_hybrid.py fresh
python3 experiments/pk-score/summarize_projected_hybrid.py
```

元のcrossing実行ファイルは `/tmp/ipknot-crossing-45e60c1` に保存してある。初期BPPや既存の調整／評価対象は `results/hybrid-r0` から継承する。生の予測・コマンド・ログ・チェックサムは `results/projected-hybrid` に保存し、公開用の集計とRNA単位の結果は `measurements/projected-hybrid-*.json` に出力する。

```sh
build/ipknot -r 0 -t auto,auto \
  --pk-h-formulation projected --pk-h-allocation blocks \
  --pk-h-max-stem 200 --pk-h-max-loop 200 --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/hybrid-r0-v1-posterior.txt \
  --pk-learned-scale 0.05 --pk-hybrid-shape \
  --pk-h-intercept -0.05 --pk-h-stem-reward 0.025 \
  --pk-h-loop-penalty 0.0075 --pk-h-weight 1 input.fa
```

倍率0.5では `--pk-learned-scale 0.025 --pk-h-weight 0.5` にする。スコアは任意指定で、デフォルトは補正0。旧部分ステムprojectedもそのまま使用できる。特徴学習用の出力は従来通りcrossingのみで、projectedとは併用できない。
