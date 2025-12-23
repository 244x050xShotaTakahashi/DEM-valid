---
name: ColorByCharge
overview: "`scripts/plot_snapshot.py` の着色を `radius` から `charge` に切り替え可能にし、`charge` は0中心の発散カラーマップで可視化します。"
todos:
  - id: add-color-by-arg
    content: "`scripts/plot_snapshot.py` に `--color-by radius|charge` を追加（デフォルト radius）し、help/epilog を更新"
    status: completed
  - id: refactor-color-scalar
    content: "`plot_snapshot()` の着色ロジックを `radius`/`charge` 切替可能にリファクタし、カラーバーラベルも切替"
    status: completed
    dependencies:
      - add-color-by-arg
  - id: charge-diverging-norm
    content: "`charge` 用に 0中心の発散表示（TwoSlopeNorm + 対称clim）と端ケース処理を追加"
    status: completed
    dependencies:
      - refactor-color-scalar
---

# `plot_snapshot.py` 帯電量カラー対応プラン

## 目的

- `particles.csv` の `charge` 列を使って粒子を色分けできるようにする（従来の `radius` 色分けも維持）。
- `charge` は **0中心の発散カラーマップ**で、正負が直感的に分かる表示にする。

## 変更対象

- [`scripts/plot_snapshot.py`](/LARGE0/gr20001/b37581/DEM-valid/scripts/plot_snapshot.py)

## 実装方針

- **CLIに `--color-by` を追加**: `radius` / `charge` を選択（デフォルトは `radius`）。
- **着色に使うスカラー列を抽象化**: `plot_snapshot()` 内で `radii` と別に `color_values` を作り、`radius` または `charge` を割り当てる。
- **`charge` の正負表現**:
- `matplotlib.colors.TwoSlopeNorm(vcenter=0.0, vmin=..., vmax=...)` を使用。
- `vmin/vmax` はデータの最小・最大から決めつつ、発散表示が安定するように **対称レンジ**（`v = max(abs(min), abs(max))` で `[-v, +v]`）をデフォルトにする。
- カラーマップは発散系（例: `coolwarm`）を使用。
- **カラーバー**: `--color-by` に応じてラベルを切替（`Radius [m]` / `Charge [unit]`）。
- **例外/端ケース**:
- `charge` 列が無いCSVに対しては、分かりやすいエラーメッセージで終了。
- 全粒子の `charge` が同一（vmin==vmax）などの場合は、ノーマライズが壊れないように微小幅を与えるか、単色扱いで警告を出す。

## 動作確認（手元で想定）

- 既存互換:
- `python scripts/plot_snapshot.py --time 0.01`（従来通り `radius` 着色）
- 新機能:
- `python scripts/plot_snapshot.py --time 0.01 --color-by charge`
- `python scripts/plot_snapshot.py --step 100000 --color-by charge --output snapshot_charge.png`

## 追加で入れると便利（任意・小さめ）

- `--cmap`（デフォルト: `viridis`/`coolwarm` を `color-by` に応じて自動）