import argparse
import shutil
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

DATA_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "data" / "graph11.d"
WALLS_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "inputs" / "walls.dat"

def read_walls_data(filename: Union[str, Path] = WALLS_FILE_DEFAULT):
    """walls.dat を読み込み、斜面壁のリストを返す。
    
    各行は x_start z_start x_end z_end の形式。
    コメント行（#や!で始まる）や空行はスキップされる。
    """
    walls = []
    try:
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                # コメント行や空行をスキップ
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                
                # コメント部分を削除
                if '#' in line:
                    line = line.split('#')[0].strip()
                if '!' in line:
                    line = line.split('!')[0].strip()
                
                try:
                    values = [float(x) for x in line.split()]
                    if len(values) >= 4:
                        walls.append({
                            'x_start': values[0],
                            'z_start': values[1],
                            'x_end': values[2],
                            'z_end': values[3]
                        })
                except ValueError:
                    continue  # 解析できない行はスキップ
    except FileNotFoundError:
        # ファイルが見つからない場合は空のリストを返す（エラーにしない）
        pass
    
    return walls

def read_simulation_data(
    filename: Union[str, Path] = DATA_FILE_DEFAULT,
    frame_step: int = 1,
    max_frames: Optional[int] = None,
):
    """graph11.d を読み込み、各タイムステップのデータを辞書のリストとして返す。

    今回はヘッダーが 2 行に分かれているケース（例: 3 つの値 + 次行に rmax のみ）にも対応する。
    読み込みはトークンストリームとして処理し、必要な個数だけ値を順次取り出す方式に変更した。
    """

    if frame_step < 1:
        raise ValueError("frame_step must be >= 1")
    if max_frames is not None and max_frames < 1:
        raise ValueError("max_frames must be >= 1 when provided")

    def _read_floats(token_buffer, fh, n):
        """token_buffer (list[str]) から n 個の float を取り出して返す。
        不足している場合はファイルから行を読み、トークンを補充する。"""
        vals = []
        while len(vals) < n:
            if not token_buffer:
                line = fh.readline()
                if not line:
                    raise EOFError("予期せぬ EOF")
                token_buffer.extend(line.split())
            # token_buffer が空でなければ pop
            vals.append(float(token_buffer.pop(0)))
        return vals

    def _skip_floats(token_buffer, fh, n):
        """token_buffer から n 個の float を取り出すが保持しない。"""
        count = 0
        while count < n:
            if not token_buffer:
                line = fh.readline()
                if not line:
                    raise EOFError("予期せぬ EOF")
                token_buffer.extend(line.split())
            token_buffer.pop(0)
            count += 1

    frames_data = []
    print("データファイルを読み込み中...")
    try:
        with open(filename, "r") as fh:
            tokens: list[str] = []  # 行を跨いで残ったトークンを一時保持
            frame_count = 0
            stored_frames = 0
            while True:
                try:
                    # ヘッダー: num_particles, time, container_width, container_height, rmax
                    num_particles_f, time_val, container_width, container_height, rmax_val = _read_floats(tokens, fh, 5)
                    num_particles = int(num_particles_f)
                except EOFError:
                    break  # 正常終了

                if num_particles < 0:
                    print("警告: num_particles が負の値です。スキップします。")
                    continue

                take_frame = (frame_count % frame_step == 0)
                particles_data_to_read = num_particles * 3

                if num_particles > 0:
                    if take_frame:
                        particle_values = _read_floats(tokens, fh, particles_data_to_read)
                    else:
                        _skip_floats(tokens, fh, particles_data_to_read)
                else:
                    particle_values = []

                # 速度データも同数だけ存在する。描画で使わないため常に読み飛ばす。
                if num_particles > 0:
                    _skip_floats(tokens, fh, particles_data_to_read)

                if num_particles > 0:
                    if take_frame:
                        rotation_angles = _read_floats(tokens, fh, num_particles)
                    else:
                        _skip_floats(tokens, fh, num_particles)
                else:
                    rotation_angles = []

                if take_frame:
                    particles = [
                        {
                            "x": particle_values[i * 3 + 0],
                            "z": particle_values[i * 3 + 1],
                            "r": particle_values[i * 3 + 2],
                            "rotation_angle": rotation_angles[i] if i < len(rotation_angles) else 0.0,
                        }
                        for i in range(num_particles)
                    ]

                    frames_data.append(
                        {
                            "time": time_val,
                            "num_particles": num_particles,
                            "container_width": container_width,
                            "container_height": container_height,
                            "particles": particles,
                        }
                    )
                    stored_frames += 1
                    if max_frames is not None and stored_frames >= max_frames:
                        print(f"\r  {len(frames_data)} フレームの読み込み完了 (制限到達)!        ")
                        break

                frame_count += 1
                if frame_count % 100 == 0:
                    print(f"  {frame_count} フレーム読み込み完了...", end='\r')
    except FileNotFoundError:
        print(f"エラー: ファイル '{filename}' が見つかりません。")
        return None
    except Exception as e:
        print(f"ファイル読み込み中にエラーが発生しました: {e}")
        return None
    
    print(f"\r  {len(frames_data)} フレームの読み込み完了!        ")
    return frames_data

def animate(frames_data, output_filename="pem_animation.mp4", walls_data=None, fps: int = 10):
    if not frames_data:
        print("アニメーションするデータがありません。")
        return

    # フォント設定を無効化して英語のみ使用
    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Z coordinate")
    fig.suptitle("PEM Validation", fontsize=12)

    if walls_data is None:
        walls_data = []

    is_validation_mode = frames_data[0]['num_particles'] in (1, 2)
    is_slope_validation = False
    if is_validation_mode and frames_data[0]['num_particles'] == 1 and len(frames_data) > 10:
        if frames_data[0]['particles']:
            initial_x = frames_data[0]['particles'][0]['x']
            later_frame_idx = min(10, len(frames_data) - 1)
            if frames_data[later_frame_idx]['particles']:
                later_x = frames_data[later_frame_idx]['particles'][0]['x']
                x_movement = abs(later_x - initial_x)
                if x_movement > 0.5:
                    is_slope_validation = True

    all_widths = [frame['container_width'] for frame in frames_data]
    max_container_width = max(all_widths) if all_widths else 1.0
    min_container_width = min(all_widths) if all_widths else 0.0

    container_heights = [frame.get('container_height', 0.0) for frame in frames_data]
    max_container_height = max(container_heights) if container_heights else 0.0

    max_z_overall = 0.0
    for frame_data in frames_data:
        if frame_data['particles']:
            current_max_z = max(p['z'] + p['r'] for p in frame_data['particles'])
            if current_max_z > max_z_overall:
                max_z_overall = current_max_z

    if max_z_overall <= 0:
        reference_width = frames_data[0]['container_width']
        max_z_overall = reference_width * 0.5 if reference_width > 0 else 1.0

    display_height = max(max_z_overall, max_container_height)
    if is_validation_mode and max_container_height <= 0.0:
        display_height = max(display_height, 12.0)

    x_margin = 0.1 * max_container_width if max_container_width > 0 else 1.0
    ax.set_xlim(min_container_width - x_margin, max_container_width + x_margin)

    if max_container_height > 0.0:
        y_min = -0.05 * display_height
        y_max = 1.05 * display_height
    else:
        y_min = -0.1 * display_height
        y_max = 1.2 * display_height
    ax.set_ylim(y_min, y_max)

    max_wall_height_fallback = display_height * 1.15

    time_text_obj = ax.text(0.05, 0.95, '', transform=ax.transAxes, ha="left", va="top", fontsize=10)

    left_wall = Line2D([], [], color='black', lw=2)
    bottom_wall = Line2D([], [], color='black', lw=2)
    right_wall = Line2D([], [], color='black', lw=2)
    top_wall = Line2D([], [], color='black', lw=2)
    for wall_line in (left_wall, bottom_wall, right_wall, top_wall):
        ax.add_line(wall_line)
    top_wall.set_visible(False)

    if walls_data:
        for idx, wall in enumerate(walls_data):
            line = Line2D(
                [wall['x_start'], wall['x_end']],
                [wall['z_start'], wall['z_end']],
                color='black',
                lw=2.5,
                alpha=0.8,
            )
            if idx == 0:
                line.set_label(f"Slope Walls ({len(walls_data)})")
            ax.add_line(line)
        ax.legend(loc='upper right')

    slope_line = None
    if not walls_data and is_validation_mode and is_slope_validation:
        slope_line = Line2D([], [], color='red', lw=2, alpha=0.7, label='Slope Wall')
        slope_line.set_visible(False)
        ax.add_line(slope_line)
        ax.legend(loc='upper right')

    max_particles = max((frame['num_particles'] for frame in frames_data), default=0)
    particle_patches: List[Circle] = []
    rotation_lines: List[Line2D] = []
    for _ in range(max_particles):
        circle = Circle((0.0, 0.0), 0.0, facecolor='white', edgecolor='black', linewidth=1.5, alpha=0.9)
        circle.set_visible(False)
        ax.add_patch(circle)
        particle_patches.append(circle)

        rotation_line = Line2D([], [], color='black', linewidth=1, alpha=0.8)
        rotation_line.set_visible(False)
        ax.add_line(rotation_line)
        rotation_lines.append(rotation_line)

    artists_to_update: List[object] = [time_text_obj, left_wall, bottom_wall, right_wall, top_wall]
    artists_to_update.extend(particle_patches)
    artists_to_update.extend(rotation_lines)
    if slope_line is not None:
        artists_to_update.append(slope_line)

    frame_counter = {'count': 0}
    total_frames = len(frames_data)

    def update_frame(frame_idx):
        frame_counter['count'] += 1
        if frame_counter['count'] % 10 == 0 or frame_counter['count'] == total_frames:
            progress = (frame_counter['count'] / total_frames) * 100
            print(f"  フレーム処理中: {frame_counter['count']}/{total_frames} ({progress:.1f}%)", end='\r')

        data = frames_data[frame_idx]
        current_container_width = data['container_width']
        current_container_height = data.get('container_height', 0.0)

        wall_height = current_container_height if current_container_height > 0.0 else max_wall_height_fallback
        left_wall.set_data([0, 0], [0, wall_height])
        bottom_wall.set_data([0, current_container_width], [0, 0])
        right_wall.set_data([current_container_width, current_container_width], [0, wall_height])

        if current_container_height > 0.0:
            top_wall.set_data([0, current_container_width], [current_container_height, current_container_height])
            top_wall.set_visible(True)
        else:
            top_wall.set_visible(False)

        if slope_line is not None and data['num_particles'] == 1:
            slope_angle = np.deg2rad(30.0)
            slope_x_end = min(current_container_width, max_wall_height_fallback / np.tan(slope_angle))
            slope_z_end = slope_x_end * np.tan(slope_angle)
            slope_line.set_data([0, slope_x_end], [0, slope_z_end])
            slope_line.set_visible(True)
        elif slope_line is not None:
            slope_line.set_visible(False)

        particles = data['particles']
        for idx, circle in enumerate(particle_patches):
            if idx < len(particles):
                p_data = particles[idx]
                circle.center = (p_data['x'], p_data['z'])
                circle.radius = p_data['r']
                circle.set_visible(True)

                x_center = p_data['x']
                z_center = p_data['z']
                radius = p_data['r']
                rotation_angle = p_data['rotation_angle']
                x_end = x_center + radius * np.cos(rotation_angle)
                z_end = z_center + radius * np.sin(rotation_angle)

                rotation_lines[idx].set_data([x_center, x_end], [z_center, z_end])
                rotation_lines[idx].set_visible(True)
            else:
                circle.set_visible(False)
                rotation_lines[idx].set_visible(False)

        time_text_obj.set_text(f"Time: {data['time']:.6f} s")
        return artists_to_update

    def init_frame():
        return artists_to_update

    print("\nアニメーション作成を開始します...")
    print(f"  使用フレーム数: {len(frames_data)}")
    print(f"  出力ファイル: {output_filename}")

    ani = animation.FuncAnimation(
        fig,
        update_frame,
        frames=len(frames_data),
        init_func=init_frame,
        blit=True,
        interval=150,
    )

    try:
        print("\nアニメーション保存中...")
        output_path = Path(output_filename)
        suffix = output_path.suffix.lower()
        final_output = str(output_path)
        ffmpeg_path = shutil.which("ffmpeg")
        used_ffmpeg_writer = False

        if suffix == '.gif':
            writer = animation.PillowWriter(fps=fps)
        elif ffmpeg_path:
            matplotlib.rcParams['animation.ffmpeg_path'] = ffmpeg_path
            try:
                writer = animation.FFMpegWriter(fps=fps, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])
                used_ffmpeg_writer = True
            except Exception:
                print("FFmpeg writer の初期化に失敗しました。PillowWriter へフォールバックします。")
                writer = animation.PillowWriter(fps=fps)
                final_output = str(output_path.with_suffix('.gif'))
        else:
            print("FFmpeg 実行ファイルが見つかりません。GIF 出力へフォールバックします。")
            writer = animation.PillowWriter(fps=fps)
            final_output = str(output_path.with_suffix('.gif'))

        try:
            ani.save(final_output, writer=writer, dpi=80)
        except FileNotFoundError:
            if used_ffmpeg_writer:
                print("FFmpeg 実行時にエラーが発生しました。GIF 出力へフォールバックします。")
                fallback_output = str(output_path.with_suffix('.gif'))
                fallback_writer = animation.PillowWriter(fps=fps)
                ani.save(fallback_output, writer=fallback_writer, dpi=80)
                final_output = fallback_output
            else:
                raise

        print(f"\n\nアニメーション保存完了: {final_output}")
    except KeyboardInterrupt:
        print("Animation creation was interrupted.")
        return
    except Exception as e:
        print(f"Error during animation creation: {e}")
        print("Please check if the writer dependencies are installed correctly.")

def parse_arguments():
    parser = argparse.ArgumentParser(description="PEM アニメーション作成ツール")
    parser.add_argument("data_file", nargs="?", default=str(DATA_FILE_DEFAULT), help="解析するデータファイル")
    parser.add_argument("output_file", nargs="?", default="pem_animation.mp4", help="出力ファイル名（拡張子で writer を判定）")
    parser.add_argument("walls_file", nargs="?", default=str(WALLS_FILE_DEFAULT), help="斜面壁データファイル")
    parser.add_argument("--frame-step", type=int, default=1, help="読み込み時にこのステップ間隔でフレームを抽出")
    parser.add_argument("--max-frames", type=int, default=200, help="保持する最大フレーム数（0 で無制限）")
    return parser.parse_args()


def main():
    args = parse_arguments()

    data_file = Path(args.data_file).expanduser()
    output_file = Path(args.output_file).expanduser()
    walls_file = Path(args.walls_file).expanduser()

    frame_step = max(1, args.frame_step)
    max_frames = args.max_frames if args.max_frames and args.max_frames > 0 else None

    print("=" * 60)
    print("PEM アニメーション作成ツール")
    print("=" * 60)
    print(f"データファイル: {data_file}")
    print(f"出力ファイル: {output_file}")
    print(f"斜面壁ファイル: {walls_file}")
    print(f"読み込みステップ: {frame_step}")
    print(f"最大フレーム数: {max_frames if max_frames is not None else '無制限'}")
    print("=" * 60)

    walls = read_walls_data(walls_file)
    if walls:
        print(f"斜面壁を読み込みました: {len(walls)} 本")
        for i, wall in enumerate(walls, 1):
            print(f"  壁{i}: ({wall['x_start']:.3f}, {wall['z_start']:.3f}) -> ({wall['x_end']:.3f}, {wall['z_end']:.3f})")
    else:
        print("斜面壁ファイルが見つからないか、壁が定義されていません。")
    print("=" * 60)

    all_frames_data = read_simulation_data(
        data_file,
        frame_step=frame_step,
        max_frames=max_frames,
    )

    if all_frames_data:
        animate(all_frames_data, str(output_file), walls_data=walls)
        print("\n" + "=" * 60)
        print("処理が完了しました!")
        print("=" * 60)
    else:
        print(f"\nエラー: {data_file} からデータを読み込めませんでした。")


if __name__ == "__main__":
    main()
