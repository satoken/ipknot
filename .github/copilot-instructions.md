## 要約
このリポジトリは C++17 で実装された RNA pseudoknot 予測ツール `IPknot` です。主要な開発ワークフローは CMake を使ったビルド、オプションで複数の MIP ソルバ（GLPK/Gurobi/CPLEX/SCIP/HiGHS）や ViennaRNA / NUPACK / CONTRAfold 等の確率モデルとリンクすることです。AI エージェントは次の短いルールに従って実作業（修正・追加・テスト）を行ってください。

## 重要なポイント（すぐ役立つ観察）
- ビルド: `CMakeLists.txt` が単一エントリ。デフォルトは GLPK（または他の solver が指定されればそれを使用）。CMake フラグでソルバを切替（例: `-DENABLE_GUROBI`）。
- 実行フロー: 入力 FASTA/アラインメント -> `src/fold.*` のモデルで確率行列（posterior）を計算 -> `src/ipknot.*` の `IPknot::solve` が IP 変数と制約を作成 -> ソルバに渡して構造を出力。
- モデルの抽象化: `BPEngineSeq` / `BPEngineAln`（`src/fold.h`）の factory を使ってモデルを生成（例: LinearPartition, CONTRAfold, RNAfold, Nupack）。AI が新モデルを追加するときはこれらのインタフェースを実装する。
- 制約の表現: `BPConstraints` と `StackConstraints`（`src/ipknot.h`）が主要なデータ構造。スタック制約は具体インスタンスを `StackInstance` として管理する。
- インデックス慣習: コード中で 0-based と 1-based が混在する（例: `VSVF` の sbp は 1-based の位置を使う箇所が多い）。変更時は入出力のどちらが 0/1 ベースかを必ず確認する。

## 典型的な開発ワークフロー（最小の手順）
1. 依存を満たす（ViennaRNA, chosen solver, pkg-config, pybind11 が CMake で検出されること）
2. ビルド（作業ディレクトリ: `build/`）
```
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target ipknot
```
3. スモーク実行（用意済みのテストシーケンス）
```
./ipknot ../test.fa
```
4. Docker（隔離実行）が必要なら README の Docker セクションを参照。

## コードベースで注目すべきファイル（参照先）
- `CMakeLists.txt` — ビルド設定、ソルバ切替フラグ、pybind11 埋め込み、外部ライブラリの link オプション。
- `README.md` — 利用方法、モデル、オプションの説明（コマンド例を再利用可）。
- `src/main.cpp` — CLI（`cxxopts`）とログ（`spdlog`）の設定、入出力の取り回し。
- `src/fold.h / src/fold.cpp` — 確率モデルの実装と factory、重要なデータ変換（posterior -> sbp / bp matrix）。
- `src/ipknot.h / src/ipknot.cpp` — IP モデルの生成ロジック、制約、レベル分解（pk_level）や stacking constraints の実装。
- `src/ip.*` — 抽象化された IP ラッパー（ソルバ固有の実装に影響）。
- `src/mxfold2.*` — mxfold2 連携の入り口（現在ブランチ名に関連した作業がある場合はここを確認）。

## プロジェクト固有のコーディング規約・パターン
- ログ: `spdlog` を使用。デバッグ/情報/警告ログが豊富に書かれているので新しい処理にも同様のログを追加する。
- CLI: `cxxopts` を使う。オプションの追加は `src/main.cpp` を編集。
- 例外とエラー: 古い C スタイルのエラーハンドリングと `spdlog::warn`/`info` が混在する。中断すべき致命的エラーは `throw` するかログの後に return するパターンがある。
- テスト資産: ルートに多数の `test_*.fa` と制約ファイルがある。小さな修正はこれらを用いたスモークテストで確認する。

## 依存・統合ポイント（要チェック）
- MIP ソルバ: GLPK（デフォルト）、Gurobi、CPLEX、SCIP、HiGHS。CMake フラグで有効化。CMake は `find_package` を使う。リンク時の定義（-DWITH_GLPK 等）を利用するコードパスあり。
- 確率モデル: LinearPartition（`src/linearpartition`）、CONTRAfold（`src/contrafold`）、ViennaRNA（C ヘッダ経由）、NUPACK（`src/nupack`）。モデルの追加/設定は `BPEngine*::build()` を確認。
- pybind11: 埋め込み（`pybind11::embed` を link）。テストで Python 統合が必要な場合は pybind11 の有無に注意。

## 作業時のチェックリスト（AI エージェント用）
- 新しい C++ ファイルを追加する前に `CMakeLists.txt` の `add_executable` に追加されているか確認。
- 変更はビルドしてスモーク実行する（上記のビルド手順）。ビルドエラーが出たら最小限の修正でコンパイルを通す。
- インデックスの 0/1 ベース混在に注意。テストシーケンス（`test.fa` 系）で入出力を比較する。
- 外部 solver の挙動確認は、該当 solver の CMake フラグを付けてローカルで再現するか、GLPK（デフォルト）でまず動作確認する。

---
このファイルはリポジトリに新規追加しました。既存の AI 向け指示ファイルはリポジトリ内に見つかりませんでした（`AGENT.md` などのグロブ検索を実行済み）。不明点や追記してほしい項目（例: mxfold2 の詳細設定、CI のコマンド、特定のデバッグ出力例）があれば教えてください。
