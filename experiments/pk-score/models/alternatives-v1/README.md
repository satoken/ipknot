# 新しいPKスコアの凍結モデル

係数はold282だけで決め、新規300 RNAの評価で変更していない。
精度・偽PK・時間の比較は [評価記録](../../alternatives_evaluation.md)。既定値への推奨ではない。

| モデル | 用途 | 評価した倍率 |
|---|---|---:|
| exclusion_dp.txt | 排他性補正＋オリジナルDPループprior | crossing/projected kappa=.4 |
| exclusion_flat.txt | J=0の排他性補正 | crossing kappa=.4、projected kappa=1 |
| local.txt | 従来の独立4状態局所DP | 最良相手crossing kappa=2 |
| rank-local_crossing_c0_k2_best.txt | 最良相手＋局所DPの候補順位補正 | lambda=.25 |
| rank-local_projected_c0_k2.txt | 従来局所projected kappa=2の候補順位補正 | lambda=1 |
| rank-exclusion_flat_projected_c0_k1.txt | 排他flat projected kappa=1の候補順位補正 | lambda=.5 |
| rank-exclusion_dp_projected_c0_k0.4.txt | 排他DP projected kappa=.4の候補順位補正 | lambda=.5 |

DP付き排他性補正crossingの実行例：

```sh
build/ipknot -e lpc --beam-size 100 -r 0 -t auto,auto \
  --decoder dd --dd-dp improved-beam --dd-beam 100 \
  --dd-max-iter 50 --dd-patience 0 --dd-crossing-beam 100 --dd-witnesses 16 \
  --pk-h-formulation crossing --pk-h-allocation blocks \
  --pk-h-max-stem 200 --pk-h-max-loop 200 --pk-selection-weight 0 \
  --pk-learned-model experiments/pk-score/models/alternatives-v1/exclusion_dp.txt \
  --pk-learned-scale .4 INPUT.fa
```

最良相手＋順位学習には、上のモデルをlocal.txt、kappaを2へ変更し、
`--pk-best-partner`、`--pk-rank-model experiments/pk-score/models/alternatives-v1/rank-local_crossing_c0_k2_best.txt`、
`--pk-rank-scale .25` を加える。最良相手方式は無制約DD crossing専用。
projectedの実験ではformをprojectedへ変え、最良相手オプションは使わない。
短い区間の診断には `--pk-core-width 3` を指定する。
