#!/usr/bin/env python3
"""
粒子の半径・電荷リストファイルを読み込んで分布をプロットするスクリプト

使用方法:
    python plot_distribution.py [オプション]

例:
    # デフォルト（radii.dat と charges.dat をプロット）
    python plot_distribution.py
    
    # 半径のみプロット
    python plot_distribution.py --radius-file radii.dat
    
    # 電荷のみプロット
    python plot_distribution.py --charge-file charges.dat
    
    # 両方を別ファイルで指定
    python plot_distribution.py --radius-file my_radii.dat --charge-file my_charges.dat
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple, Optional
import os

# フォント設定
import warnings
warnings.filterwarnings('ignore', message='findfont')
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

def load_data_file(filename: str) -> Tuple[np.ndarray, str]:
    """
    データファイルを読み込む
    
    Returns:
        values: 読み込んだ値の配列
        dist_info: 分布情報（ファイルヘッダから取得）
    """
    values = []
    dist_info = ""
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # 分布情報を抽出
                if '分布:' in line:
                    dist_info = line.split('分布:')[1].strip()
                continue
            try:
                values.append(float(line))
            except ValueError:
                continue
    
    return np.array(values), dist_info

def plot_radius_distribution(values: np.ndarray, dist_info: str, ax: plt.Axes, 
                             is_diameter: bool = True):
    """半径分布をプロット"""
    if is_diameter:
        # 直径から半径に変換
        radii = values / 2.0
        title_prefix = "半径"
    else:
        radii = values
        title_prefix = "半径"
    
    # mm単位に変換
    radii_mm = radii * 1000
    
    weights = np.ones_like(radii_mm) / len(radii_mm)
    # ヒストグラム
    n, bins, patches = ax.hist(radii_mm, bins=20, weights=weights, alpha=0.7,
                                color='steelblue', edgecolor='white', linewidth=0.5)

    # KDE（カーネル密度推定）のプロット
    from scipy import stats
    if len(radii_mm) > 1:
        weights = np.ones_like(radii_mm) / len(radii_mm)
        kde = stats.gaussian_kde(radii_mm, weights=weights)
        x_range = np.linspace(radii_mm.min() * 0.9, radii_mm.max() * 1.1, 200)
        # KDEを確率に変換（ヒストグラムのビン幅を掛ける）
        kde_density = kde(x_range)
        bin_width = (bins[-1] - bins[0]) / len(bins[:-1])  # ヒストグラムのビン幅
        kde_probability = kde_density * bin_width
        ax.plot(x_range, kde_probability, 'r-', linewidth=2, label='KDE')
    
    ax.set_xlabel('Radius [mm]', fontsize=11)
    ax.set_ylabel('Probability', fontsize=11) 
    ax.set_title(f'Particle Radius Distribution\n({dist_info})', fontsize=12)
    
    # 統計情報をテキストで表示
    stats_text = f'Mean: {np.mean(radii_mm):.3f} mm\n'
    stats_text += f'Std: {np.std(radii_mm):.3f} mm\n'
    stats_text += f'Min: {np.min(radii_mm):.3f} mm\n'
    stats_text += f'Max: {np.max(radii_mm):.3f} mm\n'
    stats_text += f'N: {len(radii_mm)}'
    
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

def plot_charge_distribution(values: np.ndarray, dist_info: str, ax: plt.Axes):
    """電荷分布をプロット"""
    # nC単位に変換
    charges_nC = values * 1e9
    
    weights = np.ones_like(charges_nC) / len(charges_nC)
    # ヒストグラム
    n, bins, patches = ax.hist(charges_nC, bins=20, weights=weights, alpha=0.7,
                                color='coral', edgecolor='white', linewidth=0.5)

    # 正負で色分け
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i+1]) / 2
        if bin_center < 0:
            patch.set_facecolor('cornflowerblue')
        else:
            patch.set_facecolor('coral')

    # KDE（カーネル密度推定）のプロット
    from scipy import stats
    if len(charges_nC) > 1:
        weights = np.ones_like(charges_nC) / len(charges_nC)
        kde = stats.gaussian_kde(charges_nC, weights=weights)
        x_range = np.linspace(charges_nC.min() * 1.1, charges_nC.max() * 1.1, 200)
        # KDEを確率に変換（ヒストグラムのビン幅を掛ける）
        kde_density = kde(x_range)
        bin_width = (bins[-1] - bins[0]) / len(bins[:-1])  # ヒストグラムのビン幅
        kde_probability = kde_density * bin_width
        ax.plot(x_range, kde_probability, 'darkred', linewidth=2, label='KDE')
    
    ax.set_xlabel('Charge [nC]', fontsize=11)
    ax.set_ylabel('Probability', fontsize=11) 
    ax.set_title(f'Particle Charge Distribution\n({dist_info})', fontsize=12)
    
    # ゼロラインを追加
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    # 統計情報をテキストで表示
    n_positive = np.sum(charges_nC > 0)
    n_negative = np.sum(charges_nC < 0)
    stats_text = f'Mean: {np.mean(charges_nC):.3f} nC\n'
    stats_text += f'Std: {np.std(charges_nC):.3f} nC\n'
    stats_text += f'Min: {np.min(charges_nC):.3f} nC\n'
    stats_text += f'Max: {np.max(charges_nC):.3f} nC\n'
    stats_text += f'N: {len(charges_nC)} (+:{n_positive}, -:{n_negative})'
    
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

def main():
    parser = argparse.ArgumentParser(
        description='粒子の半径・電荷リストファイルを読み込んで分布をプロット',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--radius-file', '-r', default='radii.dat',
                        help='半径（直径）リストファイル (default: radii.dat)')
    parser.add_argument('--charge-file', '-c', default='charges.dat',
                        help='電荷リストファイル (default: charges.dat)')
    parser.add_argument('--output', '-o', default='distribution_plot.png',
                        help='出力画像ファイル名 (default: distribution_plot.png)')
    parser.add_argument('--no-show', action='store_true',
                        help='プロットを表示しない（ファイル保存のみ）')
    parser.add_argument('--dpi', type=int, default=150,
                        help='出力画像のDPI (default: 150)')
    
    args = parser.parse_args()
    
    # ファイルの存在チェック
    radius_exists = Path(args.radius_file).exists()
    charge_exists = Path(args.charge_file).exists()
    
    if not radius_exists and not charge_exists:
        print(f"エラー: ファイルが見つかりません")
        print(f"  半径ファイル: {args.radius_file}")
        print(f"  電荷ファイル: {args.charge_file}")
        return 1
    
    # プロット数に応じてレイアウトを決定
    n_plots = sum([radius_exists, charge_exists])
    
    if n_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        ax_radius = axes[0]
        ax_charge = axes[1]
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        if radius_exists:
            ax_radius = ax
            ax_charge = None
        else:
            ax_radius = None
            ax_charge = ax
    
    # 半径分布のプロット
    if radius_exists:
        radii, dist_info = load_data_file(args.radius_file)
        if len(radii) > 0:
            plot_radius_distribution(radii, dist_info, ax_radius)
            print(f"半径データ読み込み: {len(radii)} 個 ({args.radius_file})")
        else:
            print(f"警告: 半径ファイルにデータがありません: {args.radius_file}")
    
    # 電荷分布のプロット
    if charge_exists:
        charges, dist_info = load_data_file(args.charge_file)
        if len(charges) > 0:
            plot_charge_distribution(charges, dist_info, ax_charge)
            print(f"電荷データ読み込み: {len(charges)} 個 ({args.charge_file})")
        else:
            print(f"警告: 電荷ファイルにデータがありません: {args.charge_file}")
    
    plt.tight_layout()
    
    # 保存
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight', facecolor='white')
    print(f"プロット保存: {args.output}")
    
    # 表示
    if not args.no_show:
        # X11が利用できない環境では表示をスキップ
        if os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'):
            plt.show()
        else:
            print("(ディスプレイが利用できないため、表示をスキップしました)")
    
    return 0

if __name__ == '__main__':
    exit(main())

