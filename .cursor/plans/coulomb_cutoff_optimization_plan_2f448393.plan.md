---
name: Coulomb Cutoff Optimization Plan
overview: クーロン力のカットオフ距離の見積もりが物理的に妥当であることを確認し、推奨される範囲（0.006m - 0.020m）でのパラメータスイープを実行するための計画。
todos:
  - id: create_script
    content: Create new job script `job_aor_cutoff_sweep_low_charge.sh` with updated parameters
    status: completed
  - id: submit_job
    content: Submit the job using sbatch
    status: completed
    dependencies:
      - create_script
---

# クーロンカットオフ距離最適化計画

## 1. 見積もりの妥当性確認

ユーザーによる見積もり計算は物理的に正しく、妥当です。

-   **物理パラメータ**: 質量 $m \simeq 1.04 \times 10^{-5}\,\mathrm{kg}$, 電荷 $q \approx 0.1 \sim 0.16\,\mathrm{nC}$
-   **基準**: 遠方1ペアのクーロン力が重力の1% ($\varepsilon=0.01$) 以下となる距離
-   **計算結果**: $r_c \approx 15.4\,\mathrm{mm}$
-   **結論**: 推奨されるスイープ範囲 ($6\,\mathrm{mm} \sim 20\,\mathrm{mm}$) は、この理論値をカバーしており適切です。

## 2. 実装計画

既存の `job_aor_cutoff_sweep.sh` をベースに、低電荷量条件に合わせた新しいジョブスクリプトを作成します。

### 変更点

-   **スクリプト名**: `job_aor_cutoff_sweep_low_charge.sh` (新規作成)
-   **スイープパラメータ ($r_c$)**:
    -   Task 1: `0.006` (約 3d)
    -   Task 2: `0.010` (約 5d)
    -   Task 3: `0.015` (約 7.5d)
    -   Task 4: `0.020` (約 10d)
-   **ベース入力ファイル**: `inputs/input_coulomb_short.dat` (変更なし)

## 3. 実行手順

1.  `job_aor_cutoff_sweep.sh` をコピーして `job_aor_cutoff_sweep_low_charge.sh` を作成
2.  `case` 文内のパラメータを上記の値に書き換え