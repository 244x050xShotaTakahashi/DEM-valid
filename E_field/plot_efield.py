#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
電場データの可視化プログラム
electro_field_2D.npzから電場とポテンシャルを読み込んでプロットする
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.font_manager as fm

# 日本語フォントの設定
# 利用可能な日本語フォントを探して設定
japanese_fonts = [
    'IPAGothic', 'IPAPGothic', 'IPAMincho', 'IPAPMincho',
    'Noto Sans CJK JP', 'Noto Serif CJK JP',
    'TakaoPGothic', 'TakaoGothic', 'VL Gothic',
    'DejaVu Sans'  # フォールバック
]

font_found = False
for font_name in japanese_fonts:
    font_list = [f.name for f in fm.fontManager.ttflist]
    if font_name in font_list:
        plt.rcParams['font.family'] = font_name
        font_found = True
        print(f"日本語フォントを設定: {font_name}")
        break

if not font_found:
    print("日本語フォントが見つかりません。英語ラベルを使用します。")
    plt.rcParams['font.family'] = 'DejaVu Sans'

# マイナス記号が文字化けしないようにする
plt.rcParams['axes.unicode_minus'] = False

# データの読み込み
data = np.load('data/electro_field_2D.npz')
x = data['x']
y = data['y']
phi = data['phi']  # ポテンシャル
Ex = data['Ex']    # x方向電場
Ey = data['Ey']    # y方向電場

# メッシュグリッドの作成
X, Y = np.meshgrid(x, y)

# 電場の大きさを計算
E_magnitude = np.sqrt(Ex**2 + Ey**2)

# 電場の2乗を計算
E_squared = Ex**2 + Ey**2

# 電場の2乗の勾配を計算: ∇(E²)
grad_E2_y, grad_E2_x = np.gradient(E_squared, y, x)
grad_E2_magnitude = np.sqrt(grad_E2_x**2 + grad_E2_y**2)

# プロットの作成
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# --- 左側：ポテンシャルのヒートマップ + ベクトル場 ---
ax1 = axes[0]
# ポテンシャルをヒートマップで表示
levels = np.linspace(np.min(phi), np.max(phi), 50)
contourf = ax1.contourf(X, Y, phi, levels=levels, cmap='viridis')
cbar1 = plt.colorbar(contourf, ax=ax1, label='ポテンシャル [V]')

# 電場ベクトルを矢印で表示（間引いて表示）
step = max(1, len(x) // 20)
ax1.quiver(X[::step, ::step], Y[::step, ::step], 
           Ex[::step, ::step], Ey[::step, ::step], 
           scale=3e4, color='white', alpha=0.7, width=0.003)

ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_title('ポテンシャル分布と電場ベクトル')
ax1.set_aspect('equal')
ax1.grid(True, alpha=0.3)

# --- 右側：電場の2乗の勾配のヒートマップ + ベクトル ---
ax2 = axes[1]
# 電場の2乗の勾配をヒートマップで表示
contourf2 = ax2.contourf(X, Y, grad_E2_magnitude, levels=50, cmap='plasma')
cbar2 = plt.colorbar(contourf2, ax=ax2, label='∇(E²)の大きさ [(V/m)²/m]')

# ポテンシャルの等高線を追加
contour_lines = ax2.contour(X, Y, phi, levels=10, colors='cyan', 
                             linewidths=0.8, alpha=0.6)
ax2.clabel(contour_lines, inline=True, fontsize=8, fmt='%1.0f V')

# ∇(E²)のベクトルを表示
# スケーリングを調整
grad_scale = np.max(grad_E2_magnitude) / np.max(E_magnitude) * 3e4
ax2.quiver(X[::step, ::step], Y[::step, ::step], 
           grad_E2_x[::step, ::step], grad_E2_y[::step, ::step], 
           scale=grad_scale, color='white', alpha=0.6, width=0.003)

ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_title('電場の2乗の勾配 ∇(E²)')
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('data/efield_visualization.png', dpi=300, bbox_inches='tight')
print("図を data/efield_visualization.png に保存しました")
plt.show()

# --- 統計情報の表示 ---
print("\n=== データ統計 ===")
print(f"グリッドサイズ: {len(x)} x {len(y)}")
print(f"領域サイズ: x=[{x.min():.4f}, {x.max():.4f}] m, y=[{y.min():.4f}, {y.max():.4f}] m")
print(f"\nポテンシャル:")
print(f"  最小値: {phi.min():.2f} V")
print(f"  最大値: {phi.max():.2f} V")
print(f"\n電場 Ex:")
print(f"  最小値: {Ex.min():.2e} V/m")
print(f"  最大値: {Ex.max():.2e} V/m")
print(f"\n電場 Ey:")
print(f"  最小値: {Ey.min():.2e} V/m")
print(f"  最大値: {Ey.max():.2e} V/m")
print(f"\n電場の大きさ:")
print(f"  最小値: {E_magnitude.min():.2e} V/m")
print(f"  最大値: {E_magnitude.max():.2e} V/m")
print(f"  平均値: {E_magnitude.mean():.2e} V/m")
print(f"\n∇(E²)の大きさ:")
print(f"  最小値: {grad_E2_magnitude.min():.2e} (V/m)²/m")
print(f"  最大値: {grad_E2_magnitude.max():.2e} (V/m)²/m")
print(f"  平均値: {grad_E2_magnitude.mean():.2e} (V/m)²/m")

