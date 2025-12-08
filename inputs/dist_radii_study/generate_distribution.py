#!/usr/bin/env python3
"""
粒子の半径・電荷リストを任意の分布関数から生成するスクリプト

使用方法:
    python generate_distribution.py [オプション]

例:
    # 正規分布で半径リスト生成（平均10mm, 標準偏差2mm, 100個）
    python generate_distribution.py --type radius --dist normal --mean 0.010 --std 0.002 -n 100 -o radii.dat
    
    # 一様分布で電荷リスト生成（-1nC〜1nC, 50個）
    python generate_distribution.py --type charge --dist uniform --min -1e-9 --max 1e-9 -n 50 -o charges.dat
    
    # 対数正規分布で半径リスト生成
    python generate_distribution.py --type radius --dist lognormal --mean 0.010 --std 0.003 -n 100 -o radii.dat
"""

import argparse
import numpy as np
from typing import Callable, Dict

def create_distribution_generator(dist_type: str, params: Dict) -> Callable:
    """分布タイプに応じた乱数生成関数を返す"""
    
    generators = {
        'normal': lambda n: np.random.normal(params['mean'], params['std'], n),
        'uniform': lambda n: np.random.uniform(params['min'], params['max'], n),
        'lognormal': lambda n: np.random.lognormal(
            np.log(params['mean']) - 0.5 * np.log(1 + (params['std']/params['mean'])**2),
            np.sqrt(np.log(1 + (params['std']/params['mean'])**2)),
            n
        ),
        'exponential': lambda n: np.random.exponential(params['mean'], n),
        'gamma': lambda n: np.random.gamma(params['shape'], params['scale'], n),
        'weibull': lambda n: params['scale'] * np.random.weibull(params['shape'], n),
        'bimodal': lambda n: np.where(
            np.random.random(n) < params.get('ratio', 0.5),
            np.random.normal(params['mean1'], params['std1'], n),
            np.random.normal(params['mean2'], params['std2'], n)
        ),
    }
    
    if dist_type not in generators:
        raise ValueError(f"未対応の分布タイプ: {dist_type}\n対応: {list(generators.keys())}")
    
    return generators[dist_type]

def generate_values(dist_type: str, params: Dict, n: int, 
                    min_val: float = None, max_val: float = None) -> np.ndarray:
    """指定した分布から値を生成"""
    generator = create_distribution_generator(dist_type, params)
    values = generator(n)
    
    # 範囲制限（オプション）
    if min_val is not None:
        values = np.maximum(values, min_val)
    if max_val is not None:
        values = np.minimum(values, max_val)
    
    return values

def save_to_file(values: np.ndarray, filename: str, value_type: str, 
                 dist_info: str, is_diameter: bool = False):
    """値をファイルに保存"""
    with open(filename, 'w') as f:
        if value_type == 'radius':
            f.write(f"# 粒子の直径リスト [m]\n")
            f.write(f"# 分布: {dist_info}\n")
            f.write(f"# 粒子数: {len(values)}\n")
            f.write(f"# 注意: 直径で指定（プログラム内で半径に変換）\n\n")
            # 半径を直径に変換して出力
            for v in values:
                f.write(f"{v * 2:.6e}\n")
        else:
            f.write(f"# 粒子の電荷リスト [C]\n")
            f.write(f"# 分布: {dist_info}\n")
            f.write(f"# 粒子数: {len(values)}\n\n")
            for v in values:
                f.write(f"{v:.6e}\n")
    
    print(f"生成完了: {filename} ({len(values)}個)")

def main():
    parser = argparse.ArgumentParser(
        description='粒子の半径・電荷リストを任意の分布から生成',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
分布タイプ:
  normal      正規分布 (--mean, --std)
  uniform     一様分布 (--min, --max)
  lognormal   対数正規分布 (--mean, --std)
  exponential 指数分布 (--mean)
  gamma       ガンマ分布 (--shape, --scale)
  weibull     ワイブル分布 (--shape, --scale)
  bimodal     二峰性分布 (--mean1, --std1, --mean2, --std2, --ratio)

例:
  # 正規分布の半径（平均7.5mm, 標準偏差1mm）
  %(prog)s --type radius --dist normal --mean 0.0075 --std 0.001 -n 100

  # 対数正規分布の半径（粒径分布によく使われる）
  %(prog)s --type radius --dist lognormal --mean 0.008 --std 0.002 -n 100

  # 二峰性の電荷分布（正負の電荷が混在）
  %(prog)s --type charge --dist bimodal --mean1 1e-9 --std1 0.2e-9 --mean2 -1e-9 --std2 0.2e-9 -n 100
        """
    )
    
    parser.add_argument('--type', choices=['radius', 'charge'], required=True,
                        help='生成するタイプ (radius/charge)')
    parser.add_argument('--dist', required=True,
                        help='分布タイプ (normal/uniform/lognormal/exponential/gamma/weibull/bimodal)')
    parser.add_argument('-n', '--num', type=int, default=100,
                        help='生成する粒子数 (default: 100)')
    parser.add_argument('-o', '--output', 
                        help='出力ファイル名 (default: radii.dat or charges.dat)')
    parser.add_argument('--seed', type=int, help='乱数シード（再現性のため）')
    
    # 分布パラメータ
    parser.add_argument('--mean', type=float, help='平均値')
    parser.add_argument('--std', type=float, help='標準偏差')
    parser.add_argument('--min', type=float, help='最小値（uniform分布）')
    parser.add_argument('--max', type=float, help='最大値（uniform分布）')
    parser.add_argument('--shape', type=float, help='形状パラメータ（gamma/weibull）')
    parser.add_argument('--scale', type=float, help='スケールパラメータ（gamma/weibull）')
    
    # 二峰性分布用
    parser.add_argument('--mean1', type=float, help='第1ピークの平均（bimodal）')
    parser.add_argument('--std1', type=float, help='第1ピークの標準偏差（bimodal）')
    parser.add_argument('--mean2', type=float, help='第2ピークの平均（bimodal）')
    parser.add_argument('--std2', type=float, help='第2ピークの標準偏差（bimodal）')
    parser.add_argument('--ratio', type=float, default=0.5, 
                        help='第1ピークの割合 (default: 0.5)')
    
    # 範囲制限
    parser.add_argument('--clip-min', type=float, help='出力値の下限')
    parser.add_argument('--clip-max', type=float, help='出力値の上限')
    
    args = parser.parse_args()
    
    # 乱数シード設定
    if args.seed is not None:
        np.random.seed(args.seed)
    
    # パラメータ辞書作成
    params = {}
    for key in ['mean', 'std', 'min', 'max', 'shape', 'scale', 
                'mean1', 'std1', 'mean2', 'std2', 'ratio']:
        val = getattr(args, key, None)
        if val is not None:
            params[key] = val
    
    # デフォルト値設定
    if args.type == 'radius':
        if 'mean' not in params and args.dist in ['normal', 'lognormal', 'exponential']:
            params['mean'] = 0.0075  # 7.5mm
        if 'std' not in params and args.dist in ['normal', 'lognormal']:
            params['std'] = 0.001  # 1mm
        if 'min' not in params and args.dist == 'uniform':
            params['min'] = 0.005  # 5mm
        if 'max' not in params and args.dist == 'uniform':
            params['max'] = 0.010  # 10mm
        clip_min = args.clip_min if args.clip_min else 1e-6
    else:  # charge
        if 'mean' not in params and args.dist in ['normal', 'lognormal', 'exponential']:
            params['mean'] = 1e-9  # 1nC
        if 'std' not in params and args.dist in ['normal', 'lognormal']:
            params['std'] = 0.2e-9  # 0.2nC
        if 'min' not in params and args.dist == 'uniform':
            params['min'] = -1e-9  # -1nC
        if 'max' not in params and args.dist == 'uniform':
            params['max'] = 1e-9  # 1nC
        clip_min = args.clip_min
    
    clip_max = args.clip_max
    
    # 値生成
    try:
        values = generate_values(args.dist, params, args.num, clip_min, clip_max)
    except Exception as e:
        print(f"エラー: {e}")
        return 1
    
    # 出力ファイル名
    if args.output:
        output_file = args.output
    else:
        output_file = 'radii.dat' if args.type == 'radius' else 'charges.dat'
    
    # 分布情報文字列
    dist_info = f"{args.dist} ({', '.join(f'{k}={v}' for k, v in params.items())})"
    
    # ファイル保存
    save_to_file(values, output_file, args.type, dist_info)
    
    # 統計表示
    print(f"\n統計情報:")
    if args.type == 'radius':
        print(f"  半径 平均: {np.mean(values)*1000:.3f} mm")
        print(f"  半径 標準偏差: {np.std(values)*1000:.3f} mm")
        print(f"  半径 範囲: {np.min(values)*1000:.3f} - {np.max(values)*1000:.3f} mm")
    else:
        print(f"  電荷 平均: {np.mean(values)*1e9:.3f} nC")
        print(f"  電荷 標準偏差: {np.std(values)*1e9:.3f} nC")
        print(f"  電荷 範囲: {np.min(values)*1e9:.3f} - {np.max(values)*1e9:.3f} nC")
    
    return 0

if __name__ == '__main__':
    exit(main())





