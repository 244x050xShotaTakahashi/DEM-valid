# DEM Charge Distribution & Angle of Repose Study

## 概要
このプロジェクトは、帯電した粉体の挙動解析、特に**帯電分布（Charge Distribution）**と**クーロン力の計算範囲（Cutoff Distance）**が**安息角（Angle of Repose, AoR）**に与える影響を評価するためのDEM（個別要素法）シミュレーションコードおよびワークフローです。

FortranによるDEMソルバーと、Slurmジョブスケジューラを用いた自動パラメータスタディ、Pythonによる解析ツールで構成されています。

## 主な機能
1. **充填シミュレーション (Phase 1)**
   - 異なる帯電分布（二峰性、正規分布、一様分布）を持つ粒子群を生成し、容器へ充填します。
   - 大きなカットオフ距離を用いて安定した初期状態を作成します。

2. **安息角計測シミュレーション (Phase 2)**
   - Phase 1で生成した充填状態から壁を取り払い、崩落後の安息角を形成させます。
   - クーロン力のカットオフ距離を段階的に変更し、計算精度（物理的妥当性）と計算コストのトレードオフを解析します。

3. **自動解析・集計**
   - 崩落後の粒子配置から安息角を自動計算。
   - 基準となる高精度計算と比較し、許容誤差範囲内で最適な（最小の）カットオフ距離を提案するツールを含みます。

## ディレクトリ構成

| ディレクトリ | 説明 |
| --- | --- |
| `src/` | Fortranソースコード (`dem_valid.f90`) および解析用Pythonツール |
| `inputs/` | シミュレーションの入力パラメータファイルおよびテンプレート |
| `scripts/` | 結果の集計、可視化を行う補助スクリプト |
| `results/` | シミュレーション結果の出力先（自動生成） |
| `log/` | ジョブ実行時の標準出力・エラーログ |
| `validation/` | 基本的な物理挙動（単一粒子、衝突、振動など）の検証ケース |

## 実行ワークフロー

シミュレーションは主に2段階のフェーズで実行されます。

### Phase 1: 初期充填 (Filling)
粒子の帯電分布を設定し、容器への充填計算を行います。

```bash
# 3種類の帯電分布(bimodal, normal, uniform)を並列計算
sbatch job_charge_dist_filling.sh
```
- **出力先**: `results/charge_distribution_study/phase1_filling/`
- **生成物**: `filled_particles.dat` (Phase 2の初期配置として使用)

### Phase 2: 安息角計測 (AoR Study)
Phase 1の結果を入力として、カットオフ距離を変えながら安息角を計測します。

```bash
# Phase 1の計算完了後に実行してください
sbatch job_charge_dist_aor.sh
```
- **出力先**: `results/charge_distribution_study/phase2_aor/`
- **パラメータ**:
    - 分布: `bimodal`, `normal`, `uniform`
    - カットオフ距離 ($r_c$): `0.02`, `0.015`, `0.01`, `0.006`, `0.003` [m]

### 解析と最適化 (Analysis)
計算終了後、結果を集計し、最適なカットオフ距離を見積もります。

**1. 結果の集約**
全ケースの安息角と計算時間をCSVにまとめます。
```bash
python3 scripts/summarize_charge_dist_aor.py
```
- 出力: `results/charge_distribution_study/summary/aor_vs_cutoff.csv`

**2. カットオフ距離の判定**
基準となる計算結果（例: $r_c=0.02\text{m}$）と比較し、許容誤差（例: $0.5^\circ$）を満たす最小のカットオフ距離を判定します。
```bash
python3 src/summarize_cutoff_sweep.py \
  --root results/charge_distribution_study/phase2_aor \
  --ref-rc 0.02 \
  --tol-deg 0.5
```

## 必要要件
- **コンパイラ**: Intel Fortran Compiler (`ifort`)
- **Python環境**: Python 3.x (`pandas`, `numpy`, `matplotlib` 等)
- **ジョブ管理**: Slurm Workload Manager

## ファイルの役割詳細

- **`src/dem_valid.f90`**: DEMソルバーのメインプログラム。
- **`src/analyze_repose_angle.py`**: 粒子座標から安息角を算出するスクリプト。
- **`src/summarize_cutoff_sweep.py`**: カットオフ距離の最適値を判定するツール。
- **`job_charge_dist_filling.sh`**: Phase 1 実行用スクリプト。
- **`job_charge_dist_aor.sh`**: Phase 2 実行用スクリプト。
