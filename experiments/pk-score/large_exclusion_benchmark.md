# DP付き相互排他補正crossing：残余Rfam全件での精度評価

**結論：大きな未使用集合でも、このスコアによるPK予測精度の向上は確認できなかった。** 全てのスコア係数を固定し、過去のモデル調整・本評価で用いたmanifestと重ならないRfam RNAを集計した。残った条件適合RNAは4,445本で、PK陽性99本、陰性4,346本だった。全データの陰性率が高い影響を分けて見るため、全件、長さ帯を合わせたbalanced subset、過去の未使用300本を含む498本の比較を集計した。

## 条件

- スコアは `IPKNOT_PK_EXCLUSION_V1`、DPループ項、倍率 κ=0.4。切片1.8523222707481422、エネルギー係数0.043312241784810926、37℃を使う。調整後に係数を変更していない。
- 予測は refinement=0、同一のLPC100初期BPP、DD improved-beam100・最大50反復・patience 0・crossing-beam100・witnesses16。比較相手は同じ設定の補正なしbaseline。
- 過去の比較実験manifestのRNA IDと完全一致配列hashを除外。先行して試走した198本だけは本評価集合に含めて全体条件で再実行し、係数選択には使っていない。12〜200 ntの残余からPK陽性を全て採用し、陰性も全て採用した。長さ帯別に陰性を1本ずつ選んだbalanced subsetも、予測を始める前に固定した。
- 先行試走198本（この99陽性を全て含む）は、全体評価前に一度結果を確認している。係数調整には使っていないが、これらのRNAについて盲検の独立再現とは呼ばない。
- 4,445本全てで初期BPPを一度計算し、キャッシュからのbaseline予測が直接計算と一致することを確認。baseline・補正ありを各3回実行し、26,670回すべて成功した。信頼区間はaccession/URS接頭辞でまとめた5,000回のpaired bootstrap。
- スコアの特徴・学習には塩基配列自体を使わない。配列は既存LPCの入力と完全一致除外のhashだけに使った。近縁配列・同一familyの分離は保証できない。

## 全4,445本

この集合は陽性99本・陰性4,346本なので、PK対F1はクラス比に強く左右される。PK対precision、recall、陰性RNAでPKを出した数も併記する。

| 指標 | 補正なし | DP付き排他crossing |
|---|---:|---:|
| 全塩基対F1 | .73164 | .72895 |
| PK塩基対F1 | .03081 | .02794 |
| PK塩基対precision | .01649 | .01495 |
| PK塩基対recall | .23378 | .21354 |
| 陰性4,346本中、PKを出したRNA | 1,800 (41.42%) | 1,814 (41.74%) |
| 合計時間（RNAごとの中央値の和） | 23.26秒 | 24.97秒 |

スコアを加えたときの全塩基対F1差は−0.00269で、95%区間は[−0.00402, −0.00140]。PK塩基対F1差は−0.00287で、95%区間[−0.00653, +0.00015]は0をわずかに含む。陰性RNAでの偽PKは14本増え、合計時間は約7.4%増えた。

## 長さ帯を合わせた新規198本

全99 PK陽性と、各陽性と同じ50 nt長さ帯から事前に選んだ陰性99本を比較した。ここでもPK塩基対F1は.25273から.22692へ下がった。差−0.02582の95%区間[−0.05074, −0.00472]は0を含まない。precisionは.27504から.24209、recallは.23378から.21354、陰性の偽PKは37本から39本になった。全塩基対F1の差−0.00382の区間[−0.00949, +0.00190]は0を含む。

## 既評価300本を含むbalanced合計498本

既に係数を固定した後に測定した300本（PK陽性150、陰性150）に、新しいbalanced subset 198本を加えた参考集計。合計249陽性・249陰性で、固定係数のまま集計した。

| 指標 | 補正なし | DP付き排他crossing |
|---|---:|---:|
| 全塩基対F1 | .65273 | .65098 |
| PK塩基対F1 | .23331 | .23206 |
| PK塩基対precision | .24778 | .24126 |
| PK塩基対recall | .22044 | .22354 |
| 陰性RNAの偽PK | 103 | 108 |

PK塩基対F1差は−0.00125、95%区間[−0.01612, +0.01285]。全塩基対F1差−0.00175の区間[−0.00594, +0.00238]とともに0を含むため、合算でも精度向上は確認できなかった。この合算は既に見た300本を含む補足結果である。また、4,445本の全件評価にも先行試走の198本が含まれるので、それらを完全な独立再現としては数えない。どの評価後にも係数は変更していない。

全件の予測記録は[圧縮per-RNAデータ](measurements/large-exclusion-benchmark-all-per-rna.json.gz)（SHA-256 `a14fad5a88a0f15ffc2dc509449b36869f5c67df9702b7f8d43bdf7bb3784446`）、集計値とソース・モデル・manifestのhashは[機械可読summary](measurements/large-exclusion-benchmark-all.json)に保存した。評価コードは [large_exclusion_benchmark.py](large_exclusion_benchmark.py)。再実行には `/tmp/ipknot-dataset.zip` と既存の評価manifestが必要。

```sh
python3.12 experiments/pk-score/large_exclusion_benchmark.py prepare
python3.12 experiments/pk-score/large_exclusion_benchmark.py cache
python3.12 experiments/pk-score/large_exclusion_benchmark.py evaluate
python3.12 experiments/pk-score/large_exclusion_benchmark.py summary
```
