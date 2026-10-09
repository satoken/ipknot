# PKスコア / DD実装の統合（2026-10-09）

「整数計画法にシュードノットスコアを追加」チャットの実装を、現在の
`dev` の作業ツリーへ統合した。既存の未コミット変更を保護するため、
Gitのmerge commitではなくファイル単位で三者統合し、通常の実装コミットとして保存する。

## 統合元と範囲

- 元worktree: `/tmp/ipknot-pk-score`
- 元ブランチ: `codex/pk-score`
- 元HEAD: `ef8ab8f2ccea92f4bf0092403910dedfc34bbd15`
- 統合前の統合先HEAD: `61900da938a346591aa228c6774982dad40c2dd6`
- ファイル単位の三者比較には、以前のNMR変更を保存した `7d560bb` を使用。

取り込んだものは `src/` のPKスコア・DD・DD/NOE実装、疎なBPP集約、
ビルド設定、関連テスト、CMakeが必要とする診断用C++プログラム、
PKモデル、READMEと説明資料。元worktreeの未コミットモデル
`experiments/pk-score/models/sequence-free-log-v1-T0.2.txt` も取り込んだ。
元のソースとworktreeは変更していない。

大量の入力RNA・計測JSON・実験スクリプト・発表資料は取り込まず、元worktreeに残した。
コピーした実験説明にはそれらへのリンクがあり、数値は元チャットの実験結果である。
今回の統合で精度評価をやり直したわけではない。

## こちら側で維持した変更

NOEの疎な配列/塩基対索引、選択不能候補の除去、制約の集約、
coaxial文脈の因数分解、HiGHSの行列バッファ移譲を維持した。
`--nmr-threshold-penalty-scale` と
`--nmr-coaxial-no-adjacent-bulge` も維持し、閾値選択のPKスコア項と組み合わせた。
既存のNMR索引ヘッダ3個と、その単体テスト・行列移譲テストはバックアップとの
SHA-256一致を確認した。既存のスライド・実験出力・ルートの `ipknot` は変更していない。

取り込み元に合わせ、既定デコーダは **DD** になった。
従来のILPを使うには `--decoder ilp` を指定する。
`--decoder auto` はソルバーをリンクしていればILP、なければDDを使う。
PKスコアは引き続き任意指定で、既定の特徴重みはゼロ。
DDの候補/反復制限は近似であり、ILPと同じ予測を保証するものではない。

統合時に、ソルバーなしビルドの `pk_learned_test` に不足していた
`pk_energy.cpp` のリンク依存を補った。
既存のNOEベンチマーク用モデルレコーダーは新しいIP APIに対応させ、
保存された旧バイナリのABIにも対応した。ベンチマークは明示的にILPを選び、
MakefilesとNinjaの両方のビルドから再リンクできる。

## 検証

| 検証 | 結果 |
|---|---|
| HiGHS + ViennaRNA + MXfold2、Release | ビルド成功、CTest **92/92** 成功 |
| `ENABLE_ILP=OFF`、Release | ビルド成功、CTest **20/20** 成功 |
| 統合前ILPとのランダム予測比較 | **64/64** 条件で目的値・構造が一致 |
| NOEモデルの順序付きハッシュ比較（固定構造なし） | **144/144** 条件で一致 |
| coaxial関連フィクスチャの実際のILP求解 | **29/29** 条件で終了状態・目的値・構造が一致 |
| MXfold2でのILP / DD実行 | 両方成功。DDはsoft NOE、bulge/coaxial、no-adjacent-bulgeも指定 |
| HiGHS行列移譲の大きいチェック | 10,000列・20,000行・10,000,000係数で統合前後とも成功 |
| LPC / LPVの疎なBPP経路 | beam 100、長さ2,000 / 4,000 / 8,000で実行成功 |
| `git diff --check` | 成功 |

coaxialフィクスチャの**モデルハッシュ**は5/29件のみ一致する。
固定対 `(i,j)` の逆向き `(j,i)` を別の変数として追加しない修正が統合元にあり、
候補フィルターの入力とNOEモデルのサイズが変わるためである。
これはモデルが完全一致したという意味ではない。実際の求解で上記29条件を別途比較し、
予測と目的値が一致することを確認した。
単発の疎経路チェックは線形計算量の証明や速度向上の測定ではない。

統合前ソース・差分・状態・検証ログは
`/tmp/ipknot-pk-merge-20261009/` に保存した。
バックアップは同ディレクトリの `before.tar.gz`、統合前比較バイナリは
`/tmp/ipknot-pre-pk-merge-build/ipknot`。
主要結果は `regression/summary.json`、`model-regression/random.json`、
`fixed-fixture-predictions/summary.json`、`sparse-check/measurements.json`。
`/tmp` のファイルは将来削除され得るため、永続的な実験アーカイブではない。

## この計算機でのビルド

HiGHSは `~/.local/app/highs-arm64` の **1.15.1**、
ViennaRNAは `~/.local/app/viennarna-2.7.2-arm64` の **2.7.2** を使用。
`ldd` でも当該HiGHS共有ライブラリのリンクを確認した。

```sh
PKG_CONFIG_PATH="$HOME/.local/app/viennarna-2.7.2-arm64/lib/pkgconfig" \
cmake -S . -B build-pk-merged-highs -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DENABLE_HIGHS=ON \
  -DHiGHS_INCLUDE_DIR="$HOME/.local/app/highs-arm64/include/highs" \
  -DHiGHS_LIBRARY="$HOME/.local/app/highs-arm64/lib/libhighs.so" \
  -DWITH_MXFOLD2=ON \
  -DPython_EXECUTABLE=/tmp/ipknot-nmr-mxfold2-native-venv/bin/python \
  -Dpybind11_DIR=/tmp/ipknot-nmr-mxfold2-native-venv/lib/python3.12/site-packages/pybind11/share/cmake/pybind11 \
  -DPython3_EXECUTABLE=/usr/bin/python3
cmake --build build-pk-merged-highs --parallel 8
ctest --test-dir build-pk-merged-highs --output-on-failure --parallel 4
```

通常のLPC/LPVはPython環境変数なしで実行できる。
現在の一時Python環境を使うMXfold2では次の環境を指定して動作確認した。

```sh
PYTHONHOME=/tmp/ipknot-nmr-python/cpython-3.12.14-linux-aarch64-gnu \
PYTHONPATH=/tmp/ipknot-nmr-mxfold2-native-venv/lib/python3.12/site-packages \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build-pk-merged-highs/ipknot --decoder ilp -e mxfold2 -r 0 \
  -t .2,.2 tests/data/nmr_constraints.fa
```

ソルバーなしDD版は `/tmp/ipknot-pk-merged-dd-build/ipknot` に作成した。
再構築には `-DENABLE_ILP=OFF -DWITH_MXFOLD2=OFF` を指定する。
PKスコアの利用方法と近似の違いは [PKスコアの説明](../experiments/pk-score/README.md)
を参照する。実験モデルを自動で有効にはしていない。
