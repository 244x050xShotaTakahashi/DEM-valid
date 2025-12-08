# 粒子分布生成・可視化ツール

DEMシミュレーション用の粒子半径・電荷の分布リストを生成し、可視化するためのPythonツールです。

## ファイル構成

| ファイル | 説明 |
|---------|------|
| `generate_distribution.py` | 任意の分布関数から半径・電荷リストを生成 |
| `plot_distribution.py` | 生成したリストファイルを読み込んで分布をプロット |
| `radii.dat` | 半径（直径）リストファイル（サンプル） |
| `charges.dat` | 電荷リストファイル（サンプル） |

---

## 1. generate_distribution.py - 分布リスト生成

### 基本的な使い方

```bash
python3 generate_distribution.py --type <radius|charge> --dist <分布タイプ> [オプション]
```

### 対応している分布タイプ

| 分布タイプ | 説明 | 必要なパラメータ |
|-----------|------|-----------------|
| `normal` | 正規分布 | `--mean`, `--std` |
| `uniform` | 一様分布 | `--min`, `--max` |
| `lognormal` | 対数正規分布（粒径分布に推奨） | `--mean`, `--std` |
| `exponential` | 指数分布 | `--mean` |
| `gamma` | ガンマ分布 | `--shape`, `--scale` |
| `weibull` | ワイブル分布 | `--shape`, `--scale` |
| `bimodal` | 二峰性分布（正負電荷に便利） | `--mean1`, `--std1`, `--mean2`, `--std2`, `--ratio` |

### オプション一覧

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--type` | 生成タイプ (`radius` or `charge`) | 必須 |
| `--dist` | 分布タイプ | 必須 |
| `-n`, `--num` | 生成する粒子数 | 100 |
| `-o`, `--output` | 出力ファイル名 | `radii.dat` or `charges.dat` |
| `--seed` | 乱数シード（再現性のため） | なし |
| `--clip-min` | 出力値の下限 | なし |
| `--clip-max` | 出力値の上限 | なし |

### 使用例

#### 半径リストの生成

```bash
# 正規分布（平均7.5mm, 標準偏差1mm, 100個）
python3 generate_distribution.py --type radius --dist normal --mean 0.0075 --std 0.001 -n 100 -o radii.dat

# 対数正規分布（粒径分布によく使われる）
python3 generate_distribution.py --type radius --dist lognormal --mean 0.008 --std 0.002 -n 100 -o radii.dat

# 一様分布（5mm〜10mm）
python3 generate_distribution.py --type radius --dist uniform --min 0.005 --max 0.010 -n 100 -o radii.dat

# 再現性のためシード指定
python3 generate_distribution.py --type radius --dist normal --mean 0.0075 --std 0.001 -n 100 --seed 42 -o radii.dat
```

#### 電荷リストの生成

```bash
# 正規分布（平均1nC, 標準偏差0.2nC）
python3 generate_distribution.py --type charge --dist normal --mean 1e-9 --std 0.2e-9 -n 100 -o charges.dat

# 一様分布（-1nC〜1nC）
# 注意: 負の値は = を使って指定
python3 generate_distribution.py --type charge --dist uniform --min=-1e-9 --max=1e-9 -n 100 -o charges.dat

# 二峰性分布（正負の電荷が混在）
python3 generate_distribution.py --type charge --dist bimodal --mean1=1e-9 --std1=0.2e-9 --mean2=-1e-9 --std2=0.2e-9 --ratio 0.5 -n 100 -o charges.dat
```

> **注意**: 負の値を指定する場合は `--mean2=-1e-9` のように `=` を使ってください。スペースで区切ると正しく認識されません。

---

## 2. plot_distribution.py - 分布可視化

### 基本的な使い方

```bash
python3 plot_distribution.py [オプション]
```

### オプション一覧

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `-r`, `--radius-file` | 半径リストファイル | `radii.dat` |
| `-c`, `--charge-file` | 電荷リストファイル | `charges.dat` |
| `-o`, `--output` | 出力画像ファイル名 | `distribution_plot.png` |
| `--dpi` | 出力画像のDPI | 150 |
| `--no-show` | プロット表示をスキップ | False |

### 使用例

```bash
# デフォルト（radii.dat と charges.dat をプロット）
python3 plot_distribution.py

# 出力ファイル名を指定
python3 plot_distribution.py -o my_distribution.png

# 高解像度で出力
python3 plot_distribution.py --dpi 300 -o high_res_plot.png

# 半径のみプロット（存在しないファイルを指定）
python3 plot_distribution.py --charge-file nonexistent.dat

# 別のファイルを指定
python3 plot_distribution.py -r my_radii.dat -c my_charges.dat -o custom_plot.png
```

### プロットの特徴

- **ヒストグラム**: 分布の形状を可視化
- **KDE曲線**: カーネル密度推定による滑らかな分布曲線（赤線）
- **統計情報ボックス**: 平均、標準偏差、最小/最大値、粒子数を表示
- **電荷の色分け**: 正電荷（オレンジ）、負電荷（青）で視覚的に区別
- **ゼロライン**: 電荷プロットでは0の位置に点線を表示

---

## 3. DEMシミュレーションでの使用方法

生成したファイルをDEMシミュレーションで使用するには、入力ファイル（例: `input_distribution.dat`）に以下を追加します：

```
# 半径分布: ファイルから読み込み
RADIUS_DISTRIBUTION_TYPE  file
RADIUS_LIST_FILE  inputs/radii.dat

# 電荷分布: ファイルから読み込み
CHARGE_DISTRIBUTION_TYPE  file
CHARGE_LIST_FILE  inputs/charges.dat
```

### その他の分布指定方法

ファイルを使わず、直接統計パラメータを指定することもできます：

```
# 半径: 一様分布
RADIUS_DISTRIBUTION_TYPE  uniform
RADIUS_UNIFORM_MIN  0.005
RADIUS_UNIFORM_MAX  0.010

# 電荷: 正規分布
CHARGE_DISTRIBUTION_TYPE  normal
CHARGE_NORMAL_MEAN  1.0e-9
CHARGE_NORMAL_STD   0.2e-9
```

---

## 4. ワークフロー例

```bash
# 1. 対数正規分布で半径リストを生成
python3 generate_distribution.py --type radius --dist lognormal --mean 0.008 --std 0.002 -n 200 --seed 42 -o radii.dat

# 2. 二峰性分布で電荷リストを生成
python3 generate_distribution.py --type charge --dist bimodal --mean1=1e-9 --std1=0.2e-9 --mean2=-1e-9 --std2=0.2e-9 -n 200 --seed 42 -o charges.dat

# 3. 分布をプロットして確認
python3 plot_distribution.py -o distribution_check.png

# 4. DEMシミュレーションを実行
cd ..
./build/dem_valid inputs/input_distribution.dat data
```

---

## 依存ライブラリ

```
numpy
matplotlib
scipy
```

インストール:
```bash
pip install numpy matplotlib scipy
```





