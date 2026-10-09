# 局所Boltzmann期待値補正：全評価データセット

**結論：残余4,445本ではcrossing・projectedの両方でPK塩基対F1がbaselineより有意に下がった。過去の4コホートを含む全5,645本でも、どちらもPK F1の改善は確認できない。** 係数は変更せず、局所Boltzmann補正をcrossingまたはprojectedとしてILPへ配分した。

## 固定した方式

使ったのは、過去の282本だけで調整した `IPKNOT_PK_DP_LOCAL_V1` モデルと倍率1。モデルは `intercept=1.8523222707481422`、`energy_scale=0.043312241784810926`、37℃。Dirks–Pierce (2003) のPKループ項を4状態の局所Boltzmannモデルに入れ、基準状態からの期待塩基対数の差を補正する。塩基配列の文字はスコアに使わない。

crossingは期待値補正を既存の交差接点に配分し、projectedは上位レベルの対合係数へ射影する。モデル係数、BPP、DD設定は共通で、この配分だけを比較した。refinementなし、LPC100、DD improved-beam100、最大50反復、patience 0、crossing-beam100、witnesses16。重みの再調整はしていない。

残余Rfamセット4,445本は過去の全実験IDと完全一致配列を除いて作成され、PK陽性99・陰性4,346。以前に評価した4コホート計1,200本も、過去ログの実行ファイルが今回と異なっていたため、既存LPC100キャッシュを使って現行binaryで再実行した。5,645本すべてについてbaseline、crossing、projectedを3回ずつ測り、計50,805回は全て成功した。残余4,445本ではbaseline予測hashが元の全件benchmarkと完全一致した。係数選択に使ったのは別の282本だけだが、歴史的な1,200本は以前に結果を確認済みで、残余セットの長さ整合198本も前回の試走で確認済み。よって、全5,645本を完全な盲検データとは呼ばない。

## コホート別PK塩基対F1

`Δ95%区間` はbaselineに対する差のaccession/URS接頭辞クラスタpaired bootstrap（5,000回）。過去4コホートでは全区間が0を含み、改善を確証できない。残余4,445本では両方式とも区間が負側にあり、PK F1が下がった。

| データセット | 本数（陽性/陰性） | baseline | local crossing | crossing Δ95%区間 | local projected | projected Δ95%区間 |
|---|---:|---:|---:|---:|---:|---:|
| previous300 | 300 (150/150) | .22855 | .23254 | [−.01100, +.01925] | .23077 | [−.01270, +.01795] |
| later300 | 300 (150/150) | .23214 | .23685 | [−.01496, +.02443] | .23134 | [−.01958, +.01763] |
| known_last300 | 300 (150/150) | .24583 | .25061 | [−.00984, +.01940] | .24544 | [−.01477, +.01472] |
| fresh300 | 300 (150/150) | .24494 | .24799 | [−.01357, +.01844] | .23857 | [−.02334, +.00921] |
| residual Rfam | 4,445 (99/4,346) | .03081 | .02817 | **[−.00555, −.00054]** | .02758 | **[−.00543, −.00152]** |
| 全5コホート | 5,645 (699/4,946) | .12480 | .12429 | [−.00483, +.00386] | .12114 | [−.00766, +.00051] |

全5,645本ではPK F1の差はcrossing −.00050、projected −.00366で、どちらの95%区間も0を含む。全塩基対F1はbaseline .70931に対しcrossing .70729、projected .70840。crossingの低下は有意で、projectedの差の区間は[−.00193, +.00011]だった。

## 残余4,445本の詳細

| 条件 | 全塩基対F1 | PK塩基対F1 | PK precision | PK recall | 陰性RNAでPKを予測 | 時間比 |
|---|---:|---:|---:|---:|---:|---:|
| baseline | .73164 | .03081 | .01649 | .23378 | 1,800 / 4,346 | 1.000 |
| local crossing | .72908 | .02817 | .01504 | .22121 | 1,837 / 4,346 | 1.077 |
| local projected | .73037 | .02758 | .01472 | .21842 | 1,852 / 4,346 | 1.051 |

PK F1の低下幅はcrossing −.00264（95%区間[−.00555, −.00054]）、projected −.00324（[−.00543, −.00152]）。誤ったPK予測を抑える効果もなく、陰性RNAでPKを出した数はcrossingで37本、projectedで52本増えた。長さ帯を合わせた198本でもbaseline .25273に対しcrossing .23352、projected .23177となり、両差の95%区間は負側だった。

残余セットの処理時間はBPP計算を含めず、baselineのRNAごとの中央値の合計が22.10秒、crossingが23.81秒（+7.7%）、projectedが23.23秒（+5.1%）。全5,645本では同じ指標でbaseline29.84秒、crossing32.37秒（+8.5%）、projected30.97秒（+3.8%）。

## 判断と再現

過去4コホートでは小さな正の点推定があったが、その区間は広く、残余全件では両配分ともPK F1が有意に悪化した。全コホートを合わせても改善は確認できず、今回の設定を既定値にする根拠はない。残余だけでの悪化が明確なため、次に重みを動かす場合もこの全件を使って最適化せず、新しい独立データで検証する必要がある。

集計とhashを含む[JSON](measurements/large-local-boltzmann-benchmark.json)、全5,645本の予測・反復別結果は[圧縮JSON](measurements/large-local-boltzmann-benchmark-per-rna.json.gz)。今回の実行コードは [large_local_boltzmann_benchmark.py](large_local_boltzmann_benchmark.py)。再実行には残余データと全コホートのBPPキャッシュが必要。

```sh
python3.12 experiments/pk-score/large_local_boltzmann_benchmark.py prepare
python3.12 experiments/pk-score/large_local_boltzmann_benchmark.py evaluate
python3.12 experiments/pk-score/large_local_boltzmann_benchmark.py summary
```
