---
name: AoR cutoff短縮
overview: AUTO_WALL_WITHDRAW を用いた短縮ランでクーロンカットオフ r_c を掃引し、AoR差±0.5°で収束する最小 r_c を採用して計算時間を短縮します。必要なら第2段としてクーロン用の粗いセル格子を追加し、nspanを大幅削減します。
todos:
  - id: setup-short-input
    content: AUTO_WALL_WITHDRAW=1 前提の短縮掃引用 input テンプレ（rcだけ差し替え可能）を作る
    status: completed
  - id: add-rc-array-job
    content: COULOMB_CUTOFF を 0.02/0.03/0.04/0.05 m で掃引する配列ジョブ（出力整理込み）を追加する
    status: completed
    dependencies:
      - setup-short-input
  - id: add-summarizer
    content: 各ケースの repose_angle_results.csv と timing_report.csv を集計して 1枚にまとめるスクリプトを追加する
    status: completed
    dependencies:
      - add-rc-array-job
  - id: decide-rc
    content: AoR差±0.5°基準で最小rcを決定し、本番ジョブ/入力へ反映する
    status: completed
    dependencies:
      - add-summarizer
  - id: optional-coulomb-grid
    content: rc最適化で不十分な場合に、クーロン専用粗セル格子（別cell list）を導入して nspan を削減する
    status: cancelled
    dependencies:
      - decide-rc
---

# 安息角計測の計算時間短縮（r_c最適化＋短縮ラン）

## 目的

- **安息角（AoR）の精度を保ったまま、クーロン力計算コストを下げて全体の計算時間を短縮**します。
- 判定基準はユーザー指定の **AoR差±0.5°以内**。

## 現状（コードから分かった前提）

- **クーロン力は `COULOMB_CUTOFF`（m）でカットオフ**され、`COULOMB_SHIFT_FORCE=1` なら **shifted-force**（カットオフで力が連続）です（[src/dem_valid.f90](/LARGE0/gr20001/b37581/DEM-valid/src/dem_valid.f90)）。
- セル法版クーロン力は `nspan = ceil(r_c / cell_size)` を使い、近傍セル範囲を広げて探索します（`coulomb_force_cell_cutoff_sub`、[src/dem_valid.f90](/LARGE0/gr20001/b37581/DEM-valid/src/dem_valid.f90)）。したがって **r_c を下げるほど nspan が下がり、クーロン部が強く高速化**します。
- **AUTO_WALL_WITHDRAW が既に実装済み**で、
- 充填が静止 → 自動引き抜き → 崩落開始検知 → 再静止で `goto 200` して終了

という短縮ラン向きの停止ができます（メインループのフェーズ制御、[src/dem_valid.f90](/LARGE0/gr20001/b37581/DEM-valid/src/dem_valid.f90)）。

- AoR解析はジョブ側で `src/analyze_repose_angle.py` を叩く構成です（[job_coulomb.sh](/LARGE0/gr20001/b37581/DEM-valid/job_coulomb.sh)、[job_coulomb_comparison.sh](/LARGE0/gr20001/b37581/DEM-valid/job_coulomb_comparison.sh)）。

## 方針（実装・数値実験）

### A. まずは「r_c最適化」だけでどこまで短縮できるかを確定させる（最優先）

- 一次見積りで候補を絞り（ユーザー提示の通り）、**短縮ランで AoR が収束する最小 r_c を採用**します。
- 掃引候補（まず4点）:
- **`r_c = 0.02, 0.03, 0.04, 0.05` m**
- 実行は **AUTO_WALL_WITHDRAW=1** を使い、再静止到達で自動終了させます。

### B. 追加の第2段（必要なら）：クーロン用の粗いセル格子を導入

- いまのセル法クーロンは `nspan = ceil(r_c / cell_size)` で近傍セルを広げます。
- **cell_size が接触計算向けに小さい場合**、r_c が数cmだと nspan が大きくなりやすいので、
- **クーロン専用に粗いセル格子（cell_size_coulomb ≈ r_c）**を別途持つ
- 近傍探索を `nspan≈1` に近づける

という改造で、rc最適化以上の高速化余地があります。

- ただし接触計算（`ncel_sub`/`pcont_sub`）とセルを共有しているため、**接触用セルを大きくするのは逆効果**になり得ます。別グリッド化は「必要なら」着手します。

## 具体的な実装タスク（ファイル単位）

### 1) 短縮掃引用の入力テンプレを用意

- ベースは [`inputs/input_coulomb.dat`](/LARGE0/gr20001/b37581/DEM-valid/inputs/input_coulomb.dat) を流用。
- 変更方針（テンプレ側）:
- `AUTO_WALL_WITHDRAW 1`
- `ENABLE_WALL_WITHDRAW 1`
- `WITHDRAW_*`（対象壁ID）を本番と同じに固定
- `MIN_STEPS_BEFORE_STATIC_CHECK` と `KINETIC_ENERGY_THRESHOLD` を短縮ラン用に調整（「早すぎる誤判定」を避けつつ最短化）
- `OUTPUT_INTERVAL` を粗くしてI/Oを削減（AoR解析に必要な最低限に）
- `PROFILING_SAMPLE_INTERVAL` は0（常時計測）か、一定間隔でサンプリング（軽量化）を選択

### 2) r_c=0.02–0.05 m を回す配列ジョブを追加

- [`job_coulomb_comparison.sh`](/LARGE0/gr20001/b37581/DEM-valid/job_coulomb_comparison.sh) を参考に、
- **rcを絶対値（`COULOMB_CUTOFF`）で直接上書き**するスイープを新設（倍数指定ではなくm指定）
- 各ケースの出力ディレクトリを `results/.../rc_0.03/...` のように明確化
- 各出力に `timing_report.csv`（Fortran側）と `repose_angle_results.csv`（Python側）を必ず残す

### 3) 収束判定と採用ルール（AoR差±0.5°）を自動集計

- 実装案:
- 追加の集計スクリプト（例: `src/summarize_cutoff_sweep.py`）を作り、各ケースの
    - AoR（`repose_angle_results.csv`）
    - 実時間/クーロン時間比（`timing_report.csv`）

を1枚のCSVにまとめる。

- 基準は `r_c=0.05` のAoR、各候補の **|ΔAoR| ≤ 0.5°** を満たす最小 r_c を採用。

## 数値実験プロトコル（短縮ラン）

- **固定するもの（超重要）**
- `RANDOM_SEED`
- 電荷分布（`CHARGE_DISTRIBUTION_TYPE` と `CHARGE_LIST_FILE` など）
- それ以外のDEMパラメータ（摩擦・転がり摩擦・壁条件・層数）
- **変えるもの**
- `COULOMB_CUTOFF` のみ（0.02/0.03/0.04/0.05）
- **終端条件**
- AUTO_WALL_WITHDRAW のフェーズ3で静止検出 → 自動終了（メインコードが既に `goto 200` で実装）
- **見る指標**
- AoR（最重要）
- `timing_report.csv` の `coulomb_force` の総時間と割合（最適化効果の可視化）

## 期待される短縮（概算）

- セル法クーロンは近傍探索が `nspan = ceil(r_c / cell_size)` に比例して重くなるため、**r_c を 0.05→0.03 m に落とせれば、クーロン部は概ね数倍短縮**が期待できます（実測は `timing_report.csv` で確認）。

## リスクと対策

- **静止判定が早すぎて誤停止**: `MIN_STEPS_BEFORE_STATIC_CHECK` を確保し、`KINETIC_ENERGY_THRESHOLD` を掃引前に1回だけ妥当化。
- **AoRのばらつき**: まずは同一seedで比較し、境界ケースだけ seed を2–3本追加して再確認。

## 受け入れ条件（Doneの定義）

- `r_c=0.02/0.03/0.04/0.05` の短縮掃引が回り、
- **AoR差±0.5°を満たす最小 r_c が決まり**、
- そのときの **実時間短縮（total と coulomb_force の内訳）がCSVで提示**できる。

---

## 実装Todo

- **setup-short-input**: 短縮掃引用入力テンプレ（AUTO_WALL_WITHDRAW前提）を作る
- **add-rc-array-job**: COULOMB_CUTOFF（m指定）を掃引する配列ジョブを追加
- **add-summarizer**: AoRとtiming_reportを横断集計するスクリプトを追加