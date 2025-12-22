---
name: Add Distribution Types to DEM
overview: generate_distribution.pyにある分布タイプのうち、lognormal（対数正規分布）、bimodal（二峰性分布）、exponential（指数分布）の3種類を dem_valid.f90 に追加し、半径・電荷の両方で使用可能にする。
todos:
  - id: add-params
    content: simulation_parameters_mod に exponential/bimodal 用パラメータ変数を追加
    status: completed
  - id: add-input-parsing
    content: read_input_file に新パラメータの読み込み case 文を追加
    status: completed
  - id: add-lognormal-func
    content: generate_lognormal_random サブルーチンを追加
    status: completed
  - id: add-exponential-func
    content: generate_exponential_random サブルーチンを追加
    status: completed
  - id: add-bimodal-func
    content: generate_bimodal_random サブルーチンを追加
    status: completed
  - id: update-rmax-rmin
    content: rmax/rmin 推定の select case に3分布を追加
    status: completed
  - id: update-radius-assign
    content: 半径割り当ての select case に3分布を追加
    status: completed
  - id: update-charge-assign
    content: 電荷割り当ての select case に3分布を追加
    status: completed
  - id: test-compile
    content: コンパイル確認
    status: completed
---

# 粒子分布タイプ追加実装プラン

## 追加する分布

- `lognormal`: 対数正規分布（粒径分布に適用）
- `bimodal`: 二峰性分布（正負電荷の混在に有用）
- `exponential`: 指数分布

## 変更ファイル

`src/dem_valid.f90` のみ

---

## 1. パラメータ変数の追加

`simulation_parameters_mod` 内（68-81行目付近）に以下を追加:

```fortran
! 対数正規分布パラメータ（mean/stdは既存normalと共有可能）

! 指数分布パラメータ
real(8) :: radius_exponential_mean = 0.0075d0   ! [m]
real(8) :: charge_exponential_mean = 1.0d-9     ! [C]

! 二峰性分布パラメータ
real(8) :: radius_bimodal_mean1 = 0.005d0, radius_bimodal_std1 = 0.001d0
real(8) :: radius_bimodal_mean2 = 0.010d0, radius_bimodal_std2 = 0.001d0
real(8) :: radius_bimodal_ratio = 0.5d0

real(8) :: charge_bimodal_mean1 = 1.0d-9, charge_bimodal_std1 = 0.2d-9
real(8) :: charge_bimodal_mean2 = -1.0d-9, charge_bimodal_std2 = 0.2d-9
real(8) :: charge_bimodal_ratio = 0.5d0
```

---

## 2. 入力ファイル読み込み処理の追加

`read_input_file` サブルーチン内（1041-1079行目付近）の `select case` に追加:

```fortran
! 指数分布
case ('RADIUS_EXPONENTIAL_MEAN')
case ('CHARGE_EXPONENTIAL_MEAN')

! 二峰性分布（半径）
case ('RADIUS_BIMODAL_MEAN1'), ('RADIUS_BIMODAL_STD1')
case ('RADIUS_BIMODAL_MEAN2'), ('RADIUS_BIMODAL_STD2')
case ('RADIUS_BIMODAL_RATIO')

! 二峰性分布（電荷）
case ('CHARGE_BIMODAL_MEAN1'), ('CHARGE_BIMODAL_STD1')
case ('CHARGE_BIMODAL_MEAN2'), ('CHARGE_BIMODAL_STD2')
case ('CHARGE_BIMODAL_RATIO')
```

---

## 3. 乱数生成関数の追加

プログラム末尾（2939行目付近）に3つのサブルーチンを追加:

```fortran
subroutine generate_lognormal_random(seed_io, mean_val, std_val, result_out)
  ! 正規分布を生成し、exp()で変換
  ! パラメータ変換: μ_ln = ln(mean) - 0.5*ln(1+(std/mean)^2)
  !                σ_ln = sqrt(ln(1+(std/mean)^2))
end subroutine

subroutine generate_exponential_random(seed_io, mean_val, result_out)
  ! 逆関数法: -mean * log(u)
end subroutine

subroutine generate_bimodal_random(seed_io, mean1, std1, mean2, std2, ratio, result_out)
  ! 確率ratioで2つの正規分布を切り替え
end subroutine
```

---

## 4. 粒子生成ロジックへの追加

### 4.1 rmax/rmin推定（1318-1336行目付近）

```fortran
case ('lognormal')
  rmax_out = radius_normal_mean + 3.0d0 * radius_normal_std
  rmin_val = max(1.0d-6, radius_normal_mean - 3.0d0 * radius_normal_std)
case ('exponential')
  rmax_out = radius_exponential_mean * 5.0d0  ! 5σ相当
  rmin_val = 1.0d-6
case ('bimodal')
  rmax_out = max(radius_bimodal_mean1, radius_bimodal_mean2) + 3.0d0 * max(std1, std2)
  rmin_val = max(1.0d-6, min(mean1, mean2) - 3.0d0 * max(std1, std2))
```

### 4.2 半径割り当て（1398-1427行目付近）

```fortran
case ('lognormal')
  call generate_lognormal_random(...)
case ('exponential')
  call generate_exponential_random(...)
case ('bimodal')
  call generate_bimodal_random(...)
```

### 4.3 電荷割り当て（1430-1446行目付近）

同様にcase文を追加

---

## 入力ファイル例（追加後）

```dat
# 対数正規分布（半径に適用）
RADIUS_DISTRIBUTION_TYPE    lognormal
RADIUS_NORMAL_MEAN          0.0075   # mean と stdを共有
RADIUS_NORMAL_STD           0.001

# 二峰性分布（電荷に適用）
CHARGE_DISTRIBUTION_TYPE    bimodal
CHARGE_BIMODAL_MEAN1        1.0e-9
CHARGE_BIMODAL_STD1         0.2e-9
CHARGE_BIMODAL_MEAN2        -1.0e-9
CHARGE_BIMODAL_STD2         0.2e-9
CHARGE_BIMODAL_RATIO        0.5
```