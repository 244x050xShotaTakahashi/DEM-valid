---
name: ""
overview: ""
todos: []
---

# セル法クーロン力高速化プラン

## 方針

- 既存のセル法（`cell_head` / `particle_cell_next` を `ncel_sub` で構築）を流用し、クーロン力のペア列挙をカットオフ半径 `rc` 以内の近傍セルみに限定する。初期はセル流用（案1）で実装し、性能が不足すればクーロン専用の粗いセル（cell_size_coulomb >= rc, 近傍8セルのみ走査：案2）を追加検討。
- 力はソフトニング付き、必要ならシフトドフォースで連続化。旧 `coulomb_force_sub` はフォールバックとして保持。
- 入力パラメータを拡張し、挙動を入力ファイルから制御できるようにする。

## 設計の芯（ブレない決めごと）

1) **セル利用の方針**

- 案1（まず実装）：接触用セルを流用し `nspan = ceil(rc / cell_size)` で (2*nspan+1)^2 セルを走査。
- 案2（必要なら切替）：クーロン専用で `cell_size_coulomb >= rc` を用意し、近傍8セルのみを走査。

2) **ペア重複排除ルール（i<j）**

- 同一セル：`i` をセルの linked-list で回し、`j = particle_cell_next(i)` から開始。
- 異なるセル：近傍セルを「前方向のみ」（例：dz>0, または dz==0 なら dx>0）で走査。

3) **力の式（softening + shifted-force を固定）**

- `r_eff^2 = r^2 + delta^2` を必ず通す。
- `F = k * qi * qj * r_vec * (1/r_eff^3 - 1/rc^3)` for `r < rc`, else 0。`COULOMB_SHIFT_FORCE` が false なら `1/r_eff^3` のみでカットオフ。

4) **重なり発散対策の運用指針**

- 接触で中心距離が小さくなっても暴れないよう softening を常時適用。
- `COULOMB_SOFTENING` デフォルトは平均粒径の数%〜10% 目安。

5) **OpenMP 方針**

- 初期は現行と同じ `atomic` で安全優先。
- 追加高速化オプションとしてスレッドローカル力配列 + 最終合算（reduction）を検討。

## 主要変更点

1) **入力パラメータ拡張**

- `simulation_parameters_mod` に以下を追加し既存 `ENABLE_COULOMB_FORCE` と併用：`COULOMB_CUTOFF` (m), `COULOMB_SOFTENING`, `COULOMB_SHIFT_FORCE` (logical), `COULOMB_USE_CELL` (logical, デフォルト ON)。
- `read_input_file` と入力サンプルに項目を追加し、未指定時の安全デフォルトを設定。

2) **セル＋カットオフ版サブルーチン追加**

- `src/dem_valid.f90` に `coulomb_force_cell_cutoff_sub(rc, ...)` を新設。
- `nspan = ceiling(rc / cell_size)`（案1）で近傍セル幅を決定し、同一セルは `j=next(i)`、異セルは前方向のみを走査して二重カウントを排除。
- 距離判定に `rc^2` を用い、`r_eff^2 = r^2 + delta^2` のソフトニングを適用。`COULOMB_SHIFT_FORCE` が有効な場合は上記シフトドフォース式を使用。
- OpenMP は初回 `atomic` を踏襲し、必要ならローカル和版を追加可能。

3) **既存ループへの組込み**

- 時間積分ステップで `coulomb_force_sub` 呼び出しを `coulomb_force_cell_cutoff_sub` に置換し、`COULOMB_USE_CELL` が false なら旧実装を使用。
- セル再構築 `ncel_sub` の呼び出しタイミングは現状を維持。

4) **最小テスト/検証**

- 小規模粒子数で `COULOMB_CUTOFF` を十分大に設定し、旧実装と力総和が一致することを確認するチェックを追加（簡易比較ルーチンまたは一時的ログ）。

## 実装TODO

- add-params: COULOMB_* パラメータを `simulation_parameters_mod` と入力読み込みに追加
- add-subroutine: `coulomb_force_cell_cutoff_sub` を実装（セル走査＋カットオフ）
- wire-call: メインステップで新旧切替を組込み
- sanity-test: 大カットオフで旧実装と一致確認