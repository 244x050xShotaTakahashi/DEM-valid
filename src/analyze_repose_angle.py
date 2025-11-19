#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
安息角測定スクリプト

DEMシミュレーションの出力データから粒子堆積物の安息角を測定します。
最終フレームの粒子座標を読み込み、表面粒子を検出して線形回帰により斜面角度を計算します。
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Union, Tuple, List, Dict
import sys
from scipy.stats import linregress


def read_final_frame(filename: Union[str, Path] = "data/graph11.d") -> Dict:
    """
    graph11.dファイルから最終フレームの粒子データを読み込む
    
    Args:
        filename: 入力ファイルパス
        
    Returns:
        粒子データを含む辞書 (x, z, r, num_particles, container_width)
    """
    filename = Path(filename)
    
    if not filename.exists():
        raise FileNotFoundError(f"データファイルが見つかりません: {filename}")
    
    last_frame = None
    
    with open(filename, "r") as f:
        tokens = []
        
        def read_floats(n):
            """n個の浮動小数点数を読み込む"""
            result = []
            while len(result) < n:
                if not tokens:
                    line = f.readline()
                    if not line:
                        raise EOFError("データが不完全です")
                    tokens.extend(line.split())
                
                result.append(float(tokens.pop(0)))
            return result
        
        while True:
            try:
                # ヘッダー: num_particles, time, container_width, container_height, rmax
                header = read_floats(5)
                num_particles = int(header[0])
                time_val = header[1]
                container_width = header[2]
                container_height = header[3]
                rmax_val = header[4]
                
                if num_particles < 0:
                    continue
                
                # 粒子データ (x, z, r) × num_particles
                particles_data_to_read = num_particles * 3
                particle_values = []
                if num_particles > 0:
                    particle_values = read_floats(particles_data_to_read)
                
                # 速度データ (読み飛ばす)
                velocity_values = []
                if num_particles > 0:
                    velocity_values = read_floats(particles_data_to_read)
                
                # 回転角度データ (読み飛ばす)
                rotation_angles = []
                if num_particles > 0:
                    rotation_angles = read_floats(num_particles)
                
                # このフレームのデータを保存
                x_coords = [particle_values[i * 3] for i in range(num_particles)]
                z_coords = [particle_values[i * 3 + 1] for i in range(num_particles)]
                radii = [particle_values[i * 3 + 2] for i in range(num_particles)]
                
                last_frame = {
                    'x': np.array(x_coords),
                    'z': np.array(z_coords),
                    'r': np.array(radii),
                    'num_particles': num_particles,
                    'time': time_val,
                    'container_width': container_width,
                    'container_height': container_height
                }
                
            except EOFError:
                break
    
    if last_frame is None:
        raise ValueError("有効なフレームデータが見つかりませんでした")
    
    return last_frame


def detect_surface_particles(x: np.ndarray, z: np.ndarray, r: np.ndarray, 
                             container_width: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    堆積物の表面粒子を検出
    
    Args:
        x: 粒子のx座標配列
        z: 粒子のz座標配列
        r: 粒子の半径配列
        container_width: 容器の幅
        
    Returns:
        左斜面の表面粒子(x, z)、右斜面の表面粒子(x, z)
    """
    # 容器の中央
    center_x = container_width / 2.0
    
    # 左右の領域に分割（中央の20%は除外）
    left_margin = 0.1
    right_margin = 0.1
    center_margin = 0.1
    
    left_region = (x < center_x - center_margin * container_width) & (x > left_margin * container_width)
    right_region = (x > center_x + center_margin * container_width) & (x < (1 - right_margin) * container_width)
    
    # 各x位置での最も高い粒子を検出
    def get_surface_points(region_mask):
        region_x = x[region_mask]
        region_z = z[region_mask]
        region_r = r[region_mask]
        
        if len(region_x) == 0:
            return np.array([]), np.array([])
        
        # x座標でビンを作成
        num_bins = 20
        x_min, x_max = region_x.min(), region_x.max()
        bins = np.linspace(x_min, x_max, num_bins)
        
        surface_x = []
        surface_z = []
        
        for i in range(len(bins) - 1):
            bin_mask = (region_x >= bins[i]) & (region_x < bins[i + 1])
            if bin_mask.sum() > 0:
                # このビン内で最も高い粒子（粒子上端）
                z_top = region_z[bin_mask] + region_r[bin_mask]
                max_idx = np.argmax(z_top)
                indices = np.where(bin_mask)[0]
                particle_idx = indices[max_idx]
                
                surface_x.append(region_x[particle_idx])
                surface_z.append(region_z[particle_idx])
        
        return np.array(surface_x), np.array(surface_z)
    
    left_surface_x, left_surface_z = get_surface_points(left_region)
    right_surface_x, right_surface_z = get_surface_points(right_region)
    
    return (left_surface_x, left_surface_z), (right_surface_x, right_surface_z)


def calculate_repose_angle(surface_x: np.ndarray, surface_z: np.ndarray, 
                           side: str = 'left') -> Tuple[float, float, float, float]:
    """
    表面粒子から安息角を計算
    
    Args:
        surface_x: 表面粒子のx座標
        surface_z: 表面粒子のz座標
        side: 'left' または 'right'
        
    Returns:
        angle_deg: 安息角（度）
        slope: 斜面の傾き
        intercept: 切片
        r_squared: 決定係数
    """
    if len(surface_x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    
    # 線形回帰
    slope, intercept, r_value, p_value, std_err = linregress(surface_x, surface_z)
    
    # 傾きから角度を計算
    angle_rad = np.arctan(abs(slope))
    angle_deg = np.degrees(angle_rad)
    
    # 左斜面の場合は正の傾き、右斜面の場合は負の傾きになるはず
    if side == 'right':
        angle_rad = np.arctan(abs(slope))
        angle_deg = np.degrees(angle_rad)
    
    r_squared = r_value ** 2
    
    return angle_deg, slope, intercept, r_squared


def plot_repose_angle(data: Dict, left_surface: Tuple, right_surface: Tuple,
                     left_fit: Tuple, right_fit: Tuple, output_file: str = None):
    """
    安息角の測定結果を可視化
    
    Args:
        data: 粒子データ
        left_surface: 左斜面の表面粒子座標
        right_surface: 右斜面の表面粒子座標
        left_fit: 左斜面のフィッティング結果 (angle, slope, intercept, r2)
        right_fit: 右斜面のフィッティング結果 (angle, slope, intercept, r2)
        output_file: 出力ファイル名
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 全粒子を描画
    x, z, r = data['x'], data['z'], data['r']
    for i in range(len(x)):
        circle = plt.Circle((x[i], z[i]), r[i], color='lightgray', ec='gray', linewidth=0.5)
        ax.add_patch(circle)
    
    # 表面粒子を強調
    left_x, left_z = left_surface
    right_x, right_z = right_surface
    
    if len(left_x) > 0:
        ax.scatter(left_x, left_z, c='red', s=50, marker='o', label='左斜面表面粒子', zorder=5)
    
    if len(right_x) > 0:
        ax.scatter(right_x, right_z, c='blue', s=50, marker='o', label='右斜面表面粒子', zorder=5)
    
    # フィッティング直線を描画
    left_angle, left_slope, left_intercept, left_r2 = left_fit
    right_angle, right_slope, right_intercept, right_r2 = right_fit
    
    if not np.isnan(left_angle) and len(left_x) > 0:
        x_fit = np.array([left_x.min(), left_x.max()])
        z_fit = left_slope * x_fit + left_intercept
        ax.plot(x_fit, z_fit, 'r--', linewidth=2, 
                label=f'左斜面: {left_angle:.1f}° (R²={left_r2:.3f})')
    
    if not np.isnan(right_angle) and len(right_x) > 0:
        x_fit = np.array([right_x.min(), right_x.max()])
        z_fit = right_slope * x_fit + right_intercept
        ax.plot(x_fit, z_fit, 'b--', linewidth=2, 
                label=f'右斜面: {right_angle:.1f}° (R²={right_r2:.3f})')
    
    # 平均安息角を計算
    angles = [a for a in [left_angle, right_angle] if not np.isnan(a)]
    if angles:
        avg_angle = np.mean(angles)
        ax.text(0.5, 0.98, f'平均安息角: {avg_angle:.2f}°', 
                transform=ax.transAxes, fontsize=14, weight='bold',
                ha='center', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel('x [m]', fontsize=12)
    ax.set_ylabel('z [m]', fontsize=12)
    ax.set_title('安息角測定結果', fontsize=14, weight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    # 容器の底を描画
    ax.axhline(y=0, color='black', linewidth=2)
    ax.axvline(x=0, color='black', linewidth=2)
    ax.axvline(x=data['container_width'], color='black', linewidth=2)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"プロットを保存しました: {output_file}")
    
    plt.close()


def analyze_repose_angle(data_file: str = "data/graph11.d", 
                        output_dir: str = "results",
                        case_name: str = "repose_angle") -> Dict:
    """
    安息角を測定してCSVに保存
    
    Args:
        data_file: 入力データファイル
        output_dir: 出力ディレクトリ
        case_name: ケース名
        
    Returns:
        測定結果を含む辞書
    """
    # データ読み込み
    print(f"データを読み込んでいます: {data_file}")
    data = read_final_frame(data_file)
    print(f"粒子数: {data['num_particles']}, 時刻: {data['time']:.4f} s")
    
    # 表面粒子を検出
    print("表面粒子を検出中...")
    left_surface, right_surface = detect_surface_particles(
        data['x'], data['z'], data['r'], data['container_width']
    )
    
    left_x, left_z = left_surface
    right_x, right_z = right_surface
    print(f"左斜面表面粒子数: {len(left_x)}, 右斜面表面粒子数: {len(right_x)}")
    
    # 安息角を計算
    print("安息角を計算中...")
    left_fit = calculate_repose_angle(left_x, left_z, side='left')
    right_fit = calculate_repose_angle(right_x, right_z, side='right')
    
    left_angle, left_slope, left_intercept, left_r2 = left_fit
    right_angle, right_slope, right_intercept, right_r2 = right_fit
    
    # 平均安息角
    angles = [a for a in [left_angle, right_angle] if not np.isnan(a)]
    avg_angle = np.mean(angles) if angles else np.nan
    std_angle = np.std(angles) if len(angles) > 1 else 0.0
    
    print(f"左斜面安息角: {left_angle:.2f}° (R²={left_r2:.3f})")
    print(f"右斜面安息角: {right_angle:.2f}° (R²={right_r2:.3f})")
    print(f"平均安息角: {avg_angle:.2f}° ± {std_angle:.2f}°")
    
    # 出力ディレクトリを作成
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # プロット作成
    plot_file = output_path / f"{case_name}_repose_angle.png"
    plot_repose_angle(data, left_surface, right_surface, left_fit, right_fit, 
                     output_file=str(plot_file))
    
    # CSVに保存
    csv_file = output_path / f"{case_name}_results.csv"
    with open(csv_file, 'w') as f:
        f.write("parameter,value\n")
        f.write(f"case_name,{case_name}\n")
        f.write(f"num_particles,{data['num_particles']}\n")
        f.write(f"simulation_time,{data['time']}\n")
        f.write(f"left_angle_deg,{left_angle}\n")
        f.write(f"left_r_squared,{left_r2}\n")
        f.write(f"right_angle_deg,{right_angle}\n")
        f.write(f"right_r_squared,{right_r2}\n")
        f.write(f"average_angle_deg,{avg_angle}\n")
        f.write(f"std_angle_deg,{std_angle}\n")
    
    print(f"結果を保存しました: {csv_file}")
    
    return {
        'left_angle': left_angle,
        'right_angle': right_angle,
        'average_angle': avg_angle,
        'std_angle': std_angle,
        'left_r2': left_r2,
        'right_r2': right_r2
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='DEMシミュレーションから安息角を測定')
    parser.add_argument('--data', '-d', default='data/graph11.d', 
                       help='入力データファイル (デフォルト: data/graph11.d)')
    parser.add_argument('--output', '-o', default='results',
                       help='出力ディレクトリ (デフォルト: results)')
    parser.add_argument('--name', '-n', default='repose_angle',
                       help='ケース名 (デフォルト: repose_angle)')
    
    args = parser.parse_args()
    
    try:
        results = analyze_repose_angle(args.data, args.output, args.name)
        print("\n=== 測定完了 ===")
        print(f"平均安息角: {results['average_angle']:.2f}°")
    except Exception as e:
        print(f"エラーが発生しました: {e}", file=sys.stderr)
        sys.exit(1)



