---
name: animate_pem CLI+charge
overview: "`src/animate_pem.py` をフラグ中心のCLIに改修し、particles.csvのchargeで粒子を着色（ヒートマップ相当）してGIF/MP4を出力できるようにします。併せてgraph11.d(legacy)も自動判定で扱えるようにします。"
todos:
  - id: cli-flags
    content: "`src/animate_pem.py` に `--file/-f`, `--output/-o`, `--walls/-w`, `--format`, `--fps` 等のフラグを追加し、既存位置引数とも両立させる"
    status: completed
  - id: input-auto-detect
    content: 入力が `particles.csv` / `graph11.d` のどちらでも動くよう `--format auto|csv|legacy` と自動判定を実装する
    status: completed
    dependencies:
      - cli-flags
  - id: charge-coloring
    content: CSVから `charge` を読み込み、`--color-by charge` 時に `TwoSlopeNorm` + `coolwarm` で粒子を着色しカラーバーを追加する
    status: completed
    dependencies:
      - input-auto-detect
  - id: docs-example
    content: （任意）`scripts/README.md` に `animate_pem.py` のフラグ版使用例を追記する
    status: completed
    dependencies:
      - cli-flags
---

# animate_pem.py のCLI

拡張と帯電可視化

## ゴール

- `src/animate_pem.py` で **入力データをコマンドライン引数（フラグ）で指定**してGIF/MP4を生成できるようにする
- `particles.csv` の `charge` を使って **粒子を発散カラーマップで着色（TwoSlopeNorm, 0中心）**し、カラーバーも表示する（ユーザー要望の「帯電量ヒートマップ」）

## 方針

- 既存の位置引数（`data_file output_file walls_file`）は残しつつ、推奨は `scripts/plot_snapshot.py` と同系のフラグに統一
- 入力は `--format auto|csv|legacy` を追加し、`auto` では拡張子やヘッダで `particles.csv` / `graph11.d` を判定
- `--color-by` と `--charge-mode` を追加し、`--color-by charge` で帯電着色＋カラーバー
- カラーマップのスケールは **全フレーム共通**（読み込み時に charge の全体 min/max を追跡して対称レンジへ）にして、動画中で色がブレないようにする

## 変更内容（主要ファイル）

- [`/LARGE0/gr20001/b37581/DEM-valid/src/animate_pem.py`](/LARGE0/gr20001/b37581/DEM-valid/src/animate_pem.py)
- `parse_arguments()` をフラグ中心に拡張（例: `--file/-f`, `--output/-o`, `--walls/-w`, `--fps`, `--format`, `--color-by`, `--charge-mode`）
- `read_simulation_data()` で `charge` を読み込み（存在しない場合は0）
- `animate()` に
    - `color_by` / `charge_mode` を受け渡し
    - `--color-by charge` 時に `TwoSlopeNorm(vcenter=0)` と `coolwarm` を使って粒子の `facecolor` を毎フレーム更新
    - 静的なカラーバーを追加
- `--format legacy` / `auto` 時は `animate_pem_legacy.py` の読み込みロジック相当を `animate_pem.py` 側へ取り込むか、関数として同居させて一本化（運用上 `run_friction_study.sh` が `graph11.d` を渡すため）
- （任意・軽微）[`/LARGE0/gr20001/b37581/DEM-valid/scripts/README.md`](/LARGE0/gr20001/b37581/DEM-valid/scripts/README.md)
- `animate_pem.py` の使用例（フラグ版）を追記して、今後迷わないようにする

## CLI 仕様（案）

- 入力/出力
- `--file/-f`: `particles.csv` または `graph11.d`
- `--output/-o`: `*.gif` / `*.mp4` など
- `--walls/-w`: `walls.dat`
- `--format`: `auto|csv|legacy`（デフォルト `auto`）
- 時間間引き
- `--frame-step`: ステップ間引き
- `--max-frames`: 最大フレーム数
- `--fps`: 出力fps
- 帯電可視化
- `--color-by`: `none|radius|charge`（デフォルト `none` で現状維持）