! メインプログラムおよびモジュール: 2次元粒子要素法 (球体モデル)

! モジュール: シミュレーション定数 (配列サイズ、数学定数)
module simulation_constants_mod
    implicit none
    integer, parameter :: ni_max = 1000  ! ni: 最大粒子数
    integer, parameter :: nj_max = 40    ! nj: 粒子ごとの最大接触点数 (粒子間10 + 壁(軸/斜面)30)
    integer, parameter :: nc_max = 20000 ! nc: グリッド内の最大セル数
    real(8), parameter :: PI_VAL = 3.141592653589793d0 ! pi: 円周率
    real(8), parameter :: GRAVITY_ACCEL = 9.80665d0    ! g: 重力加速度
end module simulation_constants_mod

! モジュール: シミュレーション制御パラメータと材料物性値
module simulation_parameters_mod
    use simulation_constants_mod, only: ni_max
    implicit none

    ! シミュレーション制御パラメータ
    real(8) :: time_step                      ! dt: 時間刻み
    real(8) :: friction_coeff_particle        ! fri: 粒子間摩擦係数
    real(8) :: friction_coeff_wall            ! frw: 壁-粒子間摩擦係数
    real(8) :: young_modulus_particle         ! e: 粒子のヤング率
    real(8) :: young_modulus_wall             ! ew: 壁のヤング率
    real(8) :: poisson_ratio_particle         ! po: 粒子のポアソン比
    real(8) :: poisson_ratio_wall             ! pow: 壁のポアソン比
    real(8) :: shear_to_normal_stiffness_ratio ! so: せん断弾性係数と法線方向弾性係数の比
    real(8) :: particle_density               ! de: 粒子の密度
    real(8) :: reference_overlap            ! 参照食い込み量 δ_ref（平均半径の5%）
    logical :: stop_when_static             ! 静止検出時に計算を停止するか（1=停止, 0=続行）

    ! 粒子生成パラメータ
    real(8) :: particle_radius_large          ! r1: 大きな粒子の半径
    real(8) :: particle_radius_small          ! r2: 小さな粒子の半径
    real(8) :: container_width                ! w: 容器の幅
    real(8) :: container_height               ! h: 容器の高さ (0.0=制限なし)
    integer :: particle_gen_layers            ! ipz: 初期粒子生成層数
    integer :: random_seed                    ! 乱数シード
    
    ! セル法アルゴリズム制御パラメータ
    logical :: disable_cell_algorithm         ! セル法アルゴリズムを無効化するフラグ
    real(8) :: cell_size_override            ! セルサイズの手動設定値 (0.0=自動計算)
    
    ! 出力制御パラメータ
    integer :: output_interval                ! 出力間隔 [ステップ]
    integer :: max_calculation_steps          ! 最大計算ステップ数
    
    ! 明示座標入力の制御
    logical :: use_explicit_positions         ! 明示座標ファイルの有無で切替
    character(len=256) :: positions_file      ! 明示座標ファイルパス
    
    ! クーロン力関連パラメータ
    real(8) :: coulomb_constant               ! k: クーロン定数 [N⋅m²/C²]
    real(8) :: default_charge                 ! デフォルトの粒子電荷 [C]
    logical :: enable_coulomb_force           ! クーロン力の有効化フラグ

    save
end module simulation_parameters_mod

! モジュール: 粒子固有データ (物理特性、運動学、力)
module particle_data_mod
    use simulation_constants_mod, only: ni_max, nj_max
    implicit none

    ! 物理特性
    real(8), dimension(ni_max) :: radius         ! rr(ni): 粒子半径
    real(8), dimension(ni_max) :: mass           ! wei(ni): 粒子質量
    real(8), dimension(ni_max) :: moment_inertia ! pmi(ni): 粒子の慣性モーメント
    real(8), dimension(ni_max) :: charge         ! q(ni): 粒子の電荷 [C]

    ! 位置と向き
    real(8), dimension(ni_max) :: x_coord        ! x0(ni): 粒子中心のx座標
    real(8), dimension(ni_max) :: z_coord        ! z0(ni): 粒子中心のz座標 (原文ではy、コードではz)
    real(8), dimension(ni_max) :: rotation_angle ! qq(ni): 粒子の回転変位 (角度)

    ! 速度 (並進および回転)
    real(8), dimension(ni_max) :: x_vel          ! u0(ni): 粒子のx方向速度
    real(8), dimension(ni_max) :: z_vel          ! v0(ni): 粒子のz方向速度
    real(8), dimension(ni_max) :: rotation_vel   ! f0(ni): 粒子の回転速度

    ! 合力とモーメント
    real(8), dimension(ni_max) :: x_force_sum    ! xf(ni): 粒子に働くx方向の合力
    real(8), dimension(ni_max) :: z_force_sum    ! zf(ni): 粒子に働くz方向の合力
    real(8), dimension(ni_max) :: moment_sum     ! mf(ni): 粒子に働くモーメント

    ! 接触力の成分と接触相手のインデックス
    real(8), dimension(ni_max, nj_max) :: normal_force_contact  ! en(ni,nj): 法線方向接触力
    real(8), dimension(ni_max, nj_max) :: shear_force_contact   ! es(ni,nj): せん断方向接触力
    integer, dimension(ni_max, nj_max) :: contact_partner_idx ! je(ni,nj): 接触点番号配列 (接触している粒子/壁のインデックスを格納)
    real(8), dimension(ni_max, nj_max) :: previous_overlap      ! 前ステップのオーバーラップ（接触開始検出用）

    ! 現時間ステップにおける増分変位 (common/dpm/ より)
    real(8), dimension(ni_max) :: x_disp_incr    ! u(ni): x方向変位増分
    real(8), dimension(ni_max) :: z_disp_incr    ! v(ni): z方向変位増分
    real(8), dimension(ni_max) :: rot_disp_incr  ! f(ni): 回転変位増分
    
    save
end module particle_data_mod

! モジュール: セル格子システムデータ
module cell_system_mod
    use simulation_constants_mod, only: ni_max, nc_max
    implicit none

    integer :: num_particles          ! n: 粒子数
    integer :: cells_x_dir            ! idx: x方向のセル数
    integer :: cells_z_dir            ! idz: z方向のセル数 (使用状況から推測)

    real(8) :: cell_size              ! c: セルの幅/サイズ

    ! 各セルに属する粒子を連結リストで保持する（セル先頭 index → next → ...）
    integer, dimension(nc_max) :: cell_head        ! そのセルで最初に登録された粒子インデックス (空=0)
    integer, dimension(ni_max) :: particle_cell_next ! 次の粒子インデックス (0 ならリスト終端)

    ! 後方互換用に「最後に登録された粒子」を保持（デバッグ用途）
    integer, dimension(nc_max) :: cell_particle_map

    integer, dimension(ni_max) :: particle_cell_idx ! 粒子iが格納されているセル番号

    save
end module cell_system_mod

! モジュール: 斜面壁データ
module wall_data_mod
    implicit none
    integer, parameter :: nw_max = 128
    integer :: num_walls = 0
    real(8), dimension(nw_max) :: wall_x_start
    real(8), dimension(nw_max) :: wall_z_start
    real(8), dimension(nw_max) :: wall_x_end
    real(8), dimension(nw_max) :: wall_z_end
    real(8), dimension(nw_max) :: wall_length
    real(8), dimension(nw_max) :: wall_tangent_x
    real(8), dimension(nw_max) :: wall_tangent_z
    real(8), dimension(nw_max) :: wall_normal_x
    real(8), dimension(nw_max) :: wall_normal_z
    character(len=256) :: walls_file = 'inputs/walls.dat'
    logical :: walls_file_exists = .false.
    save
end module wall_data_mod

! メインプログラム
program two_dimensional_pem
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    use wall_data_mod
    implicit none

    integer :: it_step, static_judge_flag          ! static_judge_flag: 静止判定フラグ
    integer :: i ! ループカウンタ用にiを宣言
    real(8) :: current_time                        ! 現在時刻
    real(8) :: rmax_particle_radius                ! fpositから返される最大粒子半径
    logical :: leapfrog_initialized                ! 蛙飛び法の初期化フラグ
    
    ! 計算時間計測用変数
    integer :: start_time, end_time, clock_rate
    real(8) :: elapsed_time
    
    ! 蛙飛び法用速度記録変数（位置更新時の速度v(t+Δt/2)）
    real(8) :: z_vel_at_position_update

    ! 計算時間計測開始
    call system_clock(start_time, clock_rate)
    
    ! inputファイルからパラメータを読み込み
    call read_input_file
    
    ! 初期位置と初期条件の設定
    call fposit_sub(rmax_particle_radius)
    call inmat_sub
    call init_sub

    current_time = 0.0d0
    leapfrog_initialized = .false.

    ! 各ステップの繰り返し計算
    do it_step = 1, max_calculation_steps
        current_time = current_time + time_step
       
        if (.not. leapfrog_initialized) then
            ! 初回のみ: v(0) → v(Δt/2) への変換
            call ncel_sub
            
            ! 全粒子の合力をクリア
            do i = 1, num_particles
                x_force_sum(i) = 0.0d0
                z_force_sum(i) = 0.0d0
                moment_sum(i) = 0.0d0
            end do
            
            do i = 1, num_particles
                ! 粒子と壁との接触力計算
                call wcont_sub(i)
                ! 粒子間の接触力計算
                call pcont_sub(i, rmax_particle_radius)
            end do
            
            ! クーロン力の計算
            call coulomb_force_sub()
            
            call nposit_leapfrog_sub(static_judge_flag, 0)
            leapfrog_initialized = .true.
        else
            ! 通常ループ: 蛙飛び法のメインステップ
            ! フェーズ1: 位置更新
            call nposit_leapfrog_sub(static_judge_flag, 1)
            
            ! 新しい位置で力を計算
            call ncel_sub
            
            ! 全粒子の合力をクリア
            do i = 1, num_particles
                x_force_sum(i) = 0.0d0
                z_force_sum(i) = 0.0d0
                moment_sum(i) = 0.0d0
            end do
            
            do i = 1, num_particles
                ! 粒子と壁との接触力計算
                call wcont_sub(i)
                ! 粒子間の接触力計算
                call pcont_sub(i, rmax_particle_radius)
            end do
            
            ! クーロン力の計算
            call coulomb_force_sub()
            
            ! フェーズ2: 速度更新（新しい位置での力を使用）
            call nposit_leapfrog_sub(static_judge_flag, 2)
        end if

        ! 静止状態の判定
        if (static_judge_flag == 1) then
            if (stop_when_static) then
                write(*,*) '静止状態に到達しました。時刻: ', current_time
                goto 200 ! シミュレーションループを抜ける
            end if
        end if

        ! 計算状況の出力
        if (mod(it_step, 10000) == 0) then
            write(*, '(A,F10.6,A,F12.6,A,F12.6)') 'Time= ', current_time, &
                                                 ' Z0(N)= ', z_coord(num_particles), &
                                                 ' V0(N)= ', z_vel(num_particles)
        end if

        ! グラフィック用データの出力
        if (it_step == 1 .or. mod(it_step, output_interval) == 0) then
            call gfout_sub(it_step, current_time, rmax_particle_radius)
        end if
    end do

200 continue ! シミュレーションループ脱出用のラベル

    ! バックアップデータの出力
    call bfout_sub

    close(10) ! data/graph11.d (グラフィックデータファイル1)
    close(11) ! data/graph21.d (グラフィックデータファイル2)
    close(13) ! data/backl.d (バックアップファイル、bfout_subで開かれていれば)

    ! 計算時間計測終了
    call system_clock(end_time)
    elapsed_time = real(end_time - start_time) / real(clock_rate)
    
    write(*,*) '================================='
    write(*,*) 'シミュレーション実行結果'
    write(*,*) '================================='
    write(*,*) '粒子数: ', num_particles
    write(*,*) '計算ステップ数: ', it_step
    write(*,*) '実行時間: ', elapsed_time, ' 秒'
    write(*,*) '1ステップあたりの平均時間: ', elapsed_time / real(it_step), ' 秒'
    write(*,*) 'コンテナ幅: ', container_width
    if (container_height > 0.0d0) then
        write(*,*) 'コンテナ高さ: ', container_height, ' (上壁あり)'
    else
        write(*,*) 'コンテナ高さ: 制限なし (上壁なし)'
    end if
    
    if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
        write(*,*) 'セル法アルゴリズム: 無効化'
        write(*,*) 'セルサイズ: ', cell_size
    else
        write(*,*) 'セル法アルゴリズム: 有効'
        write(*,*) 'セルサイズ: ', cell_size
        write(*,*) 'セル数 (X方向): ', cells_x_dir
        write(*,*) 'セル数 (Z方向): ', cells_z_dir
        write(*,*) '総セル数: ', cells_x_dir * cells_z_dir
    end if

    write(*,*) '================================='
    
    stop
contains

    !> inputファイルからパラメータを読み込むサブルーチン
    subroutine read_input_file
        use wall_data_mod
        implicit none
        character(len=256) :: line, keyword
        character(len=256) :: input_filename
        integer :: ios, unit_num
        real(8) :: value
        
        ! inputファイル名の決定（コマンドライン引数または固定名）
        if (command_argument_count() > 0) then
            call get_command_argument(1, input_filename)
            ! 相対パスの場合、inputsフォルダを追加
            if (input_filename(1:1) /= '/' .and. input_filename(1:6) /= 'inputs') then
                input_filename = 'inputs/' // trim(input_filename)
            end if
        else
            input_filename = "inputs/input_valid.dat"
        end if
        
        write(*,*) 'inputファイルを読み込み中: ', trim(input_filename)
        
        unit_num = 20
        open(unit=unit_num, file=input_filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: inputファイルを開けません: ', trim(input_filename)
            stop
        end if
        
        ! デフォルト値の設定
        time_step = 5.0d-7
        friction_coeff_particle = 0.25d0
        friction_coeff_wall = 0.17d0
        young_modulus_particle = 4.9d9
        young_modulus_wall = 3.9d9
        poisson_ratio_particle = 0.23d0
        poisson_ratio_wall = 0.25d0
        particle_density = 2.48d3
        particle_radius_large = 1.0d-2
        particle_radius_small = 5.0d-3
        container_width = 5.0d-1
        container_height = 0.0d0  ! 0.0=制限なし（上壁なし）
        particle_gen_layers = 30
        random_seed = 584287
        disable_cell_algorithm = .false.
        stop_when_static = .true.
        cell_size_override = 0.0d0
        output_interval = 50000
        max_calculation_steps = 2000000
        enable_coulomb_force = .false.
        coulomb_constant = 8.99d9  ! クーロン定数 k = 1/(4πε₀) [N⋅m²/C²]
        default_charge = 0.0d0      ! デフォルト電荷 [C]
        
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! コメント行と空行をスキップ
            if (line(1:1) == '#' .or. line(1:1) == '!' .or. len_trim(line) == 0) cycle
            
            ! キーワードと値を分離
            read(line, *, iostat=ios) keyword, value
            if (ios /= 0) cycle
            
            select case (trim(keyword))
                case ('TIME_STEP')
                    time_step = value
                case ('FRICTION_COEFF_PARTICLE')
                    friction_coeff_particle = value
                case ('FRICTION_COEFF_WALL')
                    friction_coeff_wall = value
                case ('YOUNG_MODULUS_PARTICLE')
                    young_modulus_particle = value
                case ('YOUNG_MODULUS_WALL')
                    young_modulus_wall = value
                case ('POISSON_RATIO_PARTICLE')
                    poisson_ratio_particle = value
                case ('POISSON_RATIO_WALL')
                    poisson_ratio_wall = value
                case ('PARTICLE_DENSITY')
                    particle_density = value
                case ('PARTICLE_RADIUS_LARGE')
                    particle_radius_large = value
                case ('PARTICLE_RADIUS_SMALL')
                    particle_radius_small = value
                case ('CONTAINER_WIDTH')
                    container_width = value
                case ('CONTAINER_HEIGHT')
                    container_height = value
                case ('PARTICLE_GEN_LAYERS')
                    particle_gen_layers = int(value)
                case ('RANDOM_SEED')
                    random_seed = int(value)
                case ('DISABLE_CELL_ALGORITHM')
                    disable_cell_algorithm = (int(value) == 1)
                case ('STOP_WHEN_STATIC')
                    stop_when_static = (int(value) == 1)
                case ('CELL_SIZE_OVERRIDE')
                    cell_size_override = value
                case ('MAX_CALCULATION_STEPS')
                    max_calculation_steps = int(value)
                case ('ENABLE_COULOMB_FORCE')
                    enable_coulomb_force = (int(value) == 1)
                case ('COULOMB_CONSTANT')
                    coulomb_constant = value
                case ('DEFAULT_CHARGE')
                    default_charge = value
                case default
                    write(*,*) '警告: 不明なキーワード: ', trim(keyword)
            end select
        end do
        
        close(unit_num)
        
        write(*,*) 'inputファイルの読み込み完了'
        
        ! 数値積分法の表示
        write(*,*) '数値積分法: 蛙飛び法'
        
        ! 明示座標ファイルの存在チェック（固定パス）
        positions_file = 'inputs/particle_positions.dat'
        use_explicit_positions = .false.
        inquire(file=trim(positions_file), exist=use_explicit_positions)
        if (use_explicit_positions) then
            write(*,*) '粒子配置: 明示座標ファイルを使用: ', trim(positions_file)
        else
            write(*,*) '粒子配置: 乱数生成（明示座標ファイルなし）'
        end if

        call read_walls_file
        
    end subroutine read_input_file

    !> 斜面壁ファイルを読み込むサブルーチン
    subroutine read_walls_file
        use wall_data_mod
        implicit none
        integer :: unit_num, ios, line_no
        character(len=256) :: line
        real(8) :: x1, z1, x2, z2
        real(8) :: dx, dz, length_val
        logical :: file_exists
        
        walls_file_exists = .false.
        num_walls = 0
        
        inquire(file=trim(walls_file), exist=file_exists)
        if (.not. file_exists) then
            write(*,*) '斜面壁ファイルなし: ', trim(walls_file)
            return
        end if
        
        unit_num = 22
        open(unit=unit_num, file=trim(walls_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: 斜面壁ファイルを開けません: ', trim(walls_file)
            stop 'read_walls_file: open failed'
        end if
        
        line_no = 0
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_no = line_no + 1
            
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            
            read(line, *, iostat=ios) x1, z1, x2, z2
            if (ios /= 0) then
                write(*,*) '警告: 斜面壁ファイル解析失敗 (行', line_no, '): ', trim(line)
                cycle
            end if
            
            dx = x2 - x1
            dz = z2 - z1
            length_val = sqrt(dx*dx + dz*dz)
            if (length_val < 1.0d-10) then
                write(*,*) '警告: 壁長さが短すぎるためスキップ (行', line_no, ')'
                cycle
            end if
            
            if (num_walls >= nw_max) then
                write(*,*) 'エラー: 壁数が上限 (', nw_max, ') を超過しました'
                stop 'read_walls_file: too many walls'
            end if
            
            num_walls = num_walls + 1
            wall_x_start(num_walls) = x1
            wall_z_start(num_walls) = z1
            wall_x_end(num_walls) = x2
            wall_z_end(num_walls) = z2
            wall_length(num_walls) = length_val
            wall_tangent_x(num_walls) = dx / length_val
            wall_tangent_z(num_walls) = dz / length_val
            wall_normal_x(num_walls) = -wall_tangent_z(num_walls)
            wall_normal_z(num_walls) =  wall_tangent_x(num_walls)
        end do
        
        close(unit_num)
        
        if (num_walls > 0) then
            walls_file_exists = .true.
            write(*,*) '斜面壁を読み込み: ', num_walls, ' 本 (', trim(walls_file), ')'
        else
            write(*,*) '斜面壁ファイルに有効な壁がありません: ', trim(walls_file)
        end if
        
    end subroutine read_walls_file

    !> 初期粒子配置と構成を設定するサブルーチン
    subroutine fposit_sub(rmax_out)
        
        use simulation_constants_mod, only: ni_max, PI_VAL
        use simulation_parameters_mod
        use particle_data_mod, only: radius, x_coord, z_coord, rotation_angle
        use cell_system_mod ! モジュールからセルシステム関連変数を取得
        implicit none

        real(8), intent(out) :: rmax_out ! 出力: 最大粒子半径

        integer :: i_layer, j_particle_in_layer, ipx_calc, current_particle_count
        real(8) :: r1_val, r2_val, rn_val, dx_offset, random_uniform_val
        real(8) :: rmin_val ! 宣言をここに移動
        integer :: particles_this_row ! 宣言をここに移動
        ! 明示座標入力用
        integer :: ios
        character(len=256) :: line
        real(8) :: xin, zin, rin, qin
        integer :: read_count

        r1_val = particle_radius_large
        r2_val = particle_radius_small
        
        if (use_explicit_positions) then
            ! 明示座標ファイルから読み込み
            num_particles = 0
            rmax_out = 0.0d0
            rmin_val = 1.0d99
            open(unit=21, file=trim(positions_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) 'エラー: 明示座標ファイルを開けません: ', trim(positions_file)
                stop 'fposit_sub: 位置ファイルopen失敗'
            end if
            do
                read(21,'(A)', iostat=ios) line
                if (ios /= 0) exit
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#' .or. line(1:1) == '!') cycle
                ! 電荷を含めて読み込み (4列目がなければデフォルト値を使用)
                qin = default_charge
                read(line, *, iostat=ios) xin, zin, rin, qin
                if (ios /= 0) then
                    ! 4列目がない場合は3列のみ読み込み
                    read(line, *, iostat=ios) xin, zin, rin
                    if (ios /= 0) cycle
                    qin = default_charge
                end if
                if (rin <= 0.0d0) cycle
                num_particles = num_particles + 1
                if (num_particles > ni_max) then
                    write(*,*) 'エラー: 粒子数がni_maxを超過: ', ni_max
                    stop 'fposit_sub: 粒子が多すぎます'
                end if
                x_coord(num_particles) = xin
                z_coord(num_particles) = zin
                radius(num_particles)  = rin
                charge(num_particles)  = qin
                rotation_angle(num_particles) = 0.0d0
                if (rin > rmax_out) rmax_out = rin
                if (rin < rmin_val) rmin_val = rin
            end do
            close(21)
            if (num_particles <= 0) then
                write(*,*) 'エラー: 明示座標ファイルに有効な粒子がありません: ', trim(positions_file)
                stop 'fposit_sub: 有効粒子なし'
            end if
        else
            ! 既存の乱数配置
            rmax_out = r1_val             ! 最大半径をr1_valとする
            rmin_val = r2_val             ! 最小半径をr2_valとする
            rn_val = rmax_out + 1.0d-5    ! パッキングのための有効半径
            ipx_calc = idint(container_width / (2.0d0 * rn_val)) ! 1行あたりの粒子数 (概算)
            current_particle_count = 0
            do i_layer = 1, particle_gen_layers
                if (mod(i_layer, 2) == 0) then  ! 偶数層
                    dx_offset = 2.0d0 * rn_val
                    particles_this_row = ipx_calc - 1
                else                            ! 奇数層
                    dx_offset = rn_val
                    particles_this_row = ipx_calc
                end if
                do j_particle_in_layer = 1, particles_this_row
                    call custom_random(random_seed, random_uniform_val)
                    if (random_uniform_val < 2.0d-1) cycle ! 一部の位置をスキップ
                    current_particle_count = current_particle_count + 1
                    if (current_particle_count > ni_max) then
                        write(*,*) '粒子数がni_maxを超えました: ', ni_max
                        stop 'fposit_sub: 粒子が多すぎます'
                    end if
                    num_particles = current_particle_count ! グローバルな粒子数を更新
                    x_coord(num_particles) = 2.0d0 * rn_val * (j_particle_in_layer - 1) + dx_offset
                    z_coord(num_particles) = 2.0d0 * rn_val * (i_layer - 1) + rn_val
                    rotation_angle(num_particles) = 0.0d0 ! 回転角を初期化
                    charge(num_particles) = default_charge ! 電荷を初期化
                    call custom_random(random_seed, random_uniform_val)
                    if (random_uniform_val < 0.5d0) then
                        radius(num_particles) = r1_val
                    else
                        radius(num_particles) = r2_val
                    end if
                end do
            end do
        end if
        write(*,*) '生成された粒子数: ', num_particles

        ! セルサイズ計算 (原文PDF p.35 eq 3.25: C < sqrt(2)*rmin)
        if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
            ! セル法アルゴリズムを無効化する場合
            if (cell_size_override > 0.0d0) then
                cell_size = cell_size_override
            else
                ! 計算領域全体を1つのセルとして設定
                cell_size = max(container_width, 2.0d0 * rmax_out * particle_gen_layers) * 2.0d0
            end if
            write(*,*) 'セル法アルゴリズムを無効化: cell_size = ', cell_size
        else
            ! 通常のセルサイズ計算
            if (rmin_val > 0.0d0) then
                 cell_size = rmin_val * 1.30d0 ! または入力から。原文ではrmin*1.35d0はコメントアウト
            else
                 cell_size = rmax_out * 1.30d0 ! rminが適切に定義されない場合のフォールバック
            end if
        end if
        
        if (cell_size <= 0.0d0) then
            write(*,*) "エラー: fposit_subでcell_sizeが正ではありません。"
            stop
        endif

        cells_x_dir = idint(container_width / cell_size) + 1
        
        ! cells_z_dirは粒子が到達しうる最大高さをカバーする必要がある
        ! container_heightが指定されている場合はそれを使用、そうでなければ従来の推定方法を使用
        if (container_height > 0.0d0 .and. cell_size > 0.0d0) then
            ! 上壁が指定されている場合は、その高さでセル数を計算
            cells_z_dir = idint(container_height / cell_size) + 1
        else if (num_particles > 0 .and. cell_size > 0.0d0) then
            ! 従来の方法: 生成時の最上部粒子のz座標を使用
            cells_z_dir = idint(z_coord(num_particles) / cell_size) + 10 
        else if (particle_gen_layers > 0 .and. cell_size > 0.0d0 .and. rn_val > 0.0d0) then ! 粒子がない場合でも推定
             if (particle_gen_layers > 0 .and. rn_val > 0 .and. cell_size > 0) then
                cells_z_dir = idint( (2.0d0 * rn_val * (real(particle_gen_layers) -1.0d0) + rn_val) / cell_size) + 10
             else
                cells_z_dir = 20 ! Fallback if values are still problematic
             end if
        else
            cells_z_dir = 20 ! デフォルト値 (粒子も層もない、またはセルサイズが0の場合)
        end if


        if (cells_x_dir * cells_z_dir > nc_max) then
            write(*,*) 'ncl (cell_particle_map)がオーバーフローしました!! 要求セル数: ', cells_x_dir * cells_z_dir
            stop 'fposit_sub: セル配列が小さすぎます'
        end if

    end subroutine fposit_sub

    !> 材料物性値を初期化し、定数を計算するサブルーチン
    subroutine inmat_sub
        use simulation_constants_mod, only: PI_VAL
        use simulation_parameters_mod, only: time_step, reference_overlap
        use particle_data_mod, only: radius, mass, moment_inertia 
        use cell_system_mod, only: num_particles 
        implicit none
        integer :: i
        
        ! time_step, particle_densityなどの値はモジュールsimulation_parameters_modで設定されていると仮定
        ! 粒子のポアソン比に基づいてせん断弾性係数と法線方向弾性係数の比(so)を計算
        shear_to_normal_stiffness_ratio = 1.0d0 / (2.0d0 * (1.0d0 + poisson_ratio_particle))

        do i = 1, num_particles
            ! 質量: 3D球体 V = 4/3 pi r^3。2Dディスク (面積 pi r^2)の場合、deが面密度ならば。
            ! 元のコードは2Dシミュレーションの文脈でも3D球体の体積で質量を計算しているように見える。
            mass(i) = (4.0d0 / 3.0d0) * PI_VAL * radius(i)**3 * particle_density
            
            ! 慣性モーメント: 3D球体 I = 2/5 m r^2。
            ! 元のコード: pmi(i)=8.d0/15.d0*de*pi*(rr(i)**5)
            ! これは (2/5) * (4/3 pi r^3 de) * r^2 = (2/5) * mass * r^2。球体として正しい。
            moment_inertia(i) = (8.0d0 / 15.0d0) * particle_density * PI_VAL * (radius(i)**5)

            ! 参照食い込み量 δ_ref を平均半径の5%で設定
            if (num_particles > 0) then
                reference_overlap = 0.05d0 * sum(radius(1:num_particles)) / real(num_particles, 8)
            else
                reference_overlap = 0.0d0
            end if  

        end do
    end subroutine inmat_sub

    !> 接触力関連の配列を初期化するサブルーチン
    subroutine init_sub
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer :: i, j

        ! 実際の粒子数まで繰り返す
        if (num_particles > 0) then
            do i = 1, num_particles 
                do j = 1, nj_max
                    normal_force_contact(i, j) = 0.0d0
                    shear_force_contact(i, j) = 0.0d0
                    contact_partner_idx(i, j) = 0
                    previous_overlap(i, j) = -1.0d0  ! 負値で「接触なし」を表現
                end do
            end do
            ! 増分変位を最初にゼロで初期化
            x_disp_incr(1:num_particles) = 0.0d0
            z_disp_incr(1:num_particles) = 0.0d0
            rot_disp_incr(1:num_particles) = 0.0d0
        end if
    end subroutine init_sub

    !> 近傍探索のために粒子をセルに割り当てるサブルーチン
    subroutine ncel_sub
        use simulation_constants_mod, only: nc_max
        use particle_data_mod, only: x_coord, z_coord
        use cell_system_mod
        implicit none
        integer :: i, cell_block_idx
        integer :: ix_cell, iz_cell ! 宣言をここに移動
    
        ! 連結リストをクリア
        if (nc_max > 0) cell_head(1:nc_max) = 0
        if (nc_max > 0) cell_particle_map(1:nc_max) = 0  ! デバッグ用途
        if (num_particles > 0) particle_cell_next(1:num_particles) = 0

        do i = 1, num_particles
            particle_cell_idx(i) = 0 ! 初期化
            if (cell_size <= 0.0d0) then
                write(*,*) "エラー: ncel_subでcell_sizeが0または負です。"
                stop
            endif
            if (cells_x_dir <= 0) then
                 write(*,*) "エラー: ncel_subでcells_x_dirが0または負です。"
                 stop
            endif

            ! 粒子iを含むセルの1次元インデックスを計算
            ix_cell = idint(x_coord(i) / cell_size) ! 座標が負にならないように注意
            iz_cell = idint(z_coord(i) / cell_size)

            ! インデックスが有効範囲 [0, cells_x_dir-1] および [0, cells_z_dir-1] 内にあることを保証
            ix_cell = max(0, min(ix_cell, cells_x_dir - 1))
            iz_cell = max(0, min(iz_cell, cells_z_dir - 1))
            
            cell_block_idx = iz_cell * cells_x_dir + ix_cell + 1 ! Fortranの1ベースインデックス

            if (cell_block_idx > 0 .and. cell_block_idx <= nc_max) then
                 ! 連結リストの先頭に追加
                 particle_cell_next(i) = cell_head(cell_block_idx)
                 cell_head(cell_block_idx) = i

                 ! 旧 single-map も更新（最後に登録された粒子）
                 cell_particle_map(cell_block_idx) = i

                 particle_cell_idx(i) = cell_block_idx
            else
                write(*,*) 'エラー: ncel_subで粒子', i, 'のcell_block_idxが範囲外です。'
                write(*,*) 'x0, z0, c: ', x_coord(i), z_coord(i), cell_size
                write(*,*) 'ix_cell, iz_cell, idx, computed_block_idx: ', ix_cell, iz_cell, cells_x_dir, cell_block_idx
                stop 'ncel_sub: 粒子が無効なセルインデックスにマッピングされました'
            end if
        end do
    end subroutine ncel_sub

    !> 粒子iと壁との接触力を計算するサブルーチン
    subroutine wcont_sub(particle_idx)
        use simulation_parameters_mod, only: container_width, container_height
        use particle_data_mod
        use cell_system_mod, only: num_particles
        use wall_data_mod
        use simulation_constants_mod, only: nj_max
        implicit none
        integer, intent(in) :: particle_idx ! 対象の粒子インデックス
    
        real(8) :: xi, zi, ri_particle ! 粒子iのx座標, z座標, 半径
        real(8) :: wall_angle_sin, wall_angle_cos, overlap_gap ! 壁の法線ベクトル成分, 重なり量
        integer :: wall_contact_slot_idx, wall_partner_id ! 壁の接触スロットインデックス, 壁の相手粒子インデックス
        integer :: wall_idx, slot_idx, first_sloped_slot
        real(8) :: proj_len, closest_x, closest_z
        real(8) :: vec_x, vec_z, dist_sq, dist_val
        real(8) :: dynamic_angle_sin, dynamic_angle_cos
        logical, save :: wall_slot_warning_emitted = .false.
        real(8) :: en_coeff, et_coeff ! 反発係数 (normal, tangential)

        xi = x_coord(particle_idx)
        zi = z_coord(particle_idx)
        ri_particle = radius(particle_idx)
        first_sloped_slot = 14
    
        ! 左壁 (contact_partner_idx = num_particles + 1)
        wall_contact_slot_idx = 11 ! 元のコードでの左壁用の固定スロット
        wall_partner_id = num_particles + 1
        if (xi < ri_particle) then  ! 左壁と接触
            en_coeff = 0.5d0
            et_coeff = 0.5d0
            wall_angle_sin = 0.0d0  ! 法線ベクトル成分 sin(alpha_ij) (粒子中心から壁中心へ向かうベクトル)
            wall_angle_cos = -1.0d0 ! 法線ベクトル成分 cos(alpha_ij)
            overlap_gap = ri_particle - xi ! 元のコードでは dabs(xi)、ここでは重なり量を正とする
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 下壁 (contact_partner_idx = num_particles + 2)
        wall_contact_slot_idx = 12 ! 元のコードでの下壁用の固定スロット
        wall_partner_id = num_particles + 2
        if (zi < ri_particle) then  ! 下壁と接触
            en_coeff = 0.5d0
            et_coeff = 0.5d0
            wall_angle_sin = -1.0d0
            wall_angle_cos = 0.0d0
            overlap_gap = ri_particle - zi 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 右壁 (contact_partner_idx = num_particles + 3)
        wall_contact_slot_idx = 13 ! 元のコードでの右壁用の固定スロット
        wall_partner_id = num_particles + 3
        if (xi + ri_particle > container_width) then ! 右壁と接触
            en_coeff = 0.5d0
            et_coeff = 0.5d0
            wall_angle_sin = 0.0d0
            wall_angle_cos = 1.0d0
            overlap_gap = (xi + ri_particle) - container_width 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                                        ! 接触なし
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 上壁 (contact_partner_idx = num_particles + 4)
        ! container_height > 0の場合のみ上壁を有効化
        if (container_height > 0.0d0) then
            wall_contact_slot_idx = 10 ! 上壁用の固定スロット
            wall_partner_id = num_particles + 4
            if (zi + ri_particle > container_height) then ! 上壁と接触
                en_coeff = 0.5d0
                et_coeff = 0.5d0
                wall_angle_sin = 1.0d0
                wall_angle_cos = 0.0d0
                overlap_gap = (zi + ri_particle) - container_height
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
            else                                            ! 接触なし
                normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
            end if
        end if
    
        ! 斜面壁 (任意本数、ファイルで定義)
        if (num_walls > 0 .and. first_sloped_slot <= nj_max) then
            do wall_idx = 1, num_walls
                wall_partner_id = num_particles + 4 + wall_idx
                en_coeff = 0.5d0
                et_coeff = 0.5d0
                proj_len = (xi - wall_x_start(wall_idx)) * wall_tangent_x(wall_idx) + &
                           (zi - wall_z_start(wall_idx)) * wall_tangent_z(wall_idx)
                proj_len = max(0.0d0, min(proj_len, wall_length(wall_idx)))
                closest_x = wall_x_start(wall_idx) + wall_tangent_x(wall_idx) * proj_len
                closest_z = wall_z_start(wall_idx) + wall_tangent_z(wall_idx) * proj_len
    
                vec_x = closest_x - xi
                vec_z = closest_z - zi
                dist_sq = vec_x * vec_x + vec_z * vec_z
    
                if (dist_sq > 1.0d-20) then
                    dist_val = sqrt(dist_sq)
                    dynamic_angle_cos = vec_x / dist_val
                    dynamic_angle_sin = vec_z / dist_val
                else
                    dist_val = 0.0d0
                    dynamic_angle_cos = -wall_normal_x(wall_idx)
                    dynamic_angle_sin = -wall_normal_z(wall_idx)
                end if
    
                overlap_gap = ri_particle - dist_val
    
                wall_contact_slot_idx = 0
                do slot_idx = first_sloped_slot, nj_max
                    if (contact_partner_idx(particle_idx, slot_idx) == wall_partner_id) then
                        wall_contact_slot_idx = slot_idx
                        exit
                    end if
                end do
    
                if (overlap_gap > 0.0d0) then
                    if (wall_contact_slot_idx == 0) then
                        do slot_idx = first_sloped_slot, nj_max
                            if (contact_partner_idx(particle_idx, slot_idx) == 0) then
                                wall_contact_slot_idx = slot_idx
                                exit
                            end if
                        end do
                    end if
    
                    if (wall_contact_slot_idx == 0) then
                        if (.not. wall_slot_warning_emitted) then
                            write(*,*) '警告: 斜面壁との接触を格納するスロットが不足しています (粒子', particle_idx, ')'
                            wall_slot_warning_emitted = .true.
                        end if
                        cycle
                    end if
    
                    contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
                    call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, dynamic_angle_sin, dynamic_angle_cos, &
                                  overlap_gap, en_coeff, et_coeff)
                else
                    if (wall_contact_slot_idx /= 0) then
                        normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                        shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                        contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
                    end if
                end if
            end do
        end if
    end subroutine wcont_sub

    !> 粒子iと他の粒子との接触力を計算するサブルーチン
    subroutine pcont_sub(particle_i_idx, rmax_val)
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod
        implicit none

        integer, intent(in) :: particle_i_idx     ! 対象の粒子iのインデックス
        real(8), intent(in) :: rmax_val           ! 最大粒子半径 (探索範囲に使用)

        real(8) :: xi, zi, ri_particle_i
        real(8) :: xj, zj, rj_particle_j
        real(8) :: center_dist, overlap_gap
        real(8) :: contact_angle_sin, contact_angle_cos ! 粒子iから粒子jへの法線ベクトル成分
        integer :: particle_j_idx                 ! 接触相手の粒子jのインデックス
        integer :: contact_slot_for_i_j, contact_slot_for_j_i ! 接触リストのスロット
        integer :: iz_cell_min, iz_cell_max, ix_cell_min, ix_cell_max ! セル探索範囲
        integer :: current_iz_cell, current_ix_cell, cell_block_idx
        integer :: jj, max_particle_contacts_check
        real(8) :: dx, dz               ! 位置差計算用
        integer :: curr_j_idx
        real(8) :: search_extent        ! 近傍探索半径

        max_particle_contacts_check = 10 ! 元のコードのループ(do 11 jj=1,10)から、粒子間接触は最大10個と仮定

        xi = x_coord(particle_i_idx)
        zi = z_coord(particle_i_idx)
        ri_particle_i = radius(particle_i_idx)

        ! セル格子内での探索範囲を決定
        if (cell_size <= 0.0d0) stop "pcont_sub: cell_sizeが正ではありません。"

        ! 以前は ±(2*rmax) で探索していたが、cell_size が大きい場合に隣接セル2つ分を
        ! またぐ衝突を取りこぼすことがあった。そこで探索半径を (2*rmax + cell_size) に拡張する。
        search_extent = 2.0d0 * rmax_val + cell_size

        iz_cell_max = idint((zi + search_extent) / cell_size)
        iz_cell_min = idint((zi - search_extent) / cell_size)
        ix_cell_min = idint((xi - search_extent) / cell_size)
        ix_cell_max = idint((xi + search_extent) / cell_size)

        ! セルインデックスを有効なグリッド境界内に収める
        iz_cell_max = min(iz_cell_max, cells_z_dir - 1)
        iz_cell_min = max(iz_cell_min, 0)
        ix_cell_min = max(ix_cell_min, 0)
        ix_cell_max = min(ix_cell_max, cells_x_dir - 1)
        
        if (iz_cell_max < iz_cell_min .or. ix_cell_max < ix_cell_min) then
             ! 粒子が通常の領域外にある場合や、rmax_valが小さすぎて探索範囲が無効になる場合に発生しうる
             return
        end if

        ! デバッグ出力 (検証モードで粒子1のみ)
        ! if (validation_mode .and. particle_i_idx == 1) then
        !     write(*,'(A,I2,A,I3,A,I3,A,I3,A,I3)') 'Search p',particle_i_idx,': cells x[',ix_cell_min,':',ix_cell_max,'] z[',iz_cell_min,':',iz_cell_max,']'
        !     write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') 'Position: x=',xi,', z=',zi,', search_extent=',search_extent
        ! end if

        do current_iz_cell = iz_cell_min, iz_cell_max      ! z方向のセルループ
            do current_ix_cell = ix_cell_min, ix_cell_max  ! x方向のセルループ
                if (cells_x_dir <=0) stop "pcont_sub: cells_x_dirが正ではありません。"
                cell_block_idx = current_iz_cell * cells_x_dir + current_ix_cell + 1

                if (cell_block_idx <= 0 .or. cell_block_idx > (cells_x_dir * cells_z_dir) ) cycle ! セルが範囲外ならスキップ

                curr_j_idx = cell_head(cell_block_idx) ! セルに属する最初の粒子

                do while (curr_j_idx > 0)
                    particle_j_idx = curr_j_idx

                    if (particle_j_idx == particle_i_idx) then
                        curr_j_idx = particle_cell_next(curr_j_idx)
                        cycle
                    end if

                    xj = x_coord(particle_j_idx)
                    zj = z_coord(particle_j_idx)
                    rj_particle_j = radius(particle_j_idx)

                    dx = xi - xj
                    dz = zi - zj
                    center_dist = sqrt(dx*dx + dz*dz) 
                    overlap_gap = (ri_particle_i + rj_particle_j) - center_dist ! 重なり量 (正なら接触)

                    if (overlap_gap > 0.0d0) then ! 粒子が接触している (重なっている)
                        ! デバッグ出力 (検証モードのみ)
                        ! if (validation_mode) then
                        !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        !         'Contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        ! end if

                        if (center_dist < 1.0d-12) then ! 粒子中心が完全に一致する場合のゼロ除算を回避
                            contact_angle_cos = 1.0d0   ! 暫定的にx軸方向とする
                            contact_angle_sin = 0.0d0
                        else
                            ! 法線ベクトル (i から j へ向かう方向) の成分
                            contact_angle_cos = (xj - xi) / center_dist 
                            contact_angle_sin = (zj - zi) / center_dist
                        end if
                        
                        ! ---- 連絡スロット確保 (i -> j) ----------------
                        contact_slot_for_i_j = 0
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                contact_slot_for_i_j = jj
                                exit
                            end if
                        end do
                        if (contact_slot_for_i_j == 0) then
                            do jj = 1, max_particle_contacts_check
                                if (contact_partner_idx(particle_i_idx, jj) == 0) then
                                    contact_slot_for_i_j = jj
                                    contact_partner_idx(particle_i_idx, jj) = particle_j_idx
                                    exit
                                end if
                            end do
                        end if
                        if (contact_slot_for_i_j == 0) then
                            curr_j_idx = particle_cell_next(curr_j_idx)
                            cycle  ! スロット不足
                        end if

                        ! -------------------------------------------------

                        ! 実際の力計算
                        call actf_sub(particle_i_idx, particle_j_idx, contact_slot_for_i_j, &
                                      contact_angle_sin, contact_angle_cos, overlap_gap, 0.5d0, 0.5d0)
                    else ! 幾何学的な接触なし / 粒子が離れた
                        ! デバッグ出力 (検証モードのみ、最初の数回のみ)
                        ! if (validation_mode .and. particle_i_idx == 1 .and. particle_j_idx == 2) then
                        !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        !         'No contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        ! end if

                        ! particle_i_idx の particle_j_idx に関する接触情報をクリア
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                normal_force_contact(particle_i_idx, jj) = 0.0d0
                                shear_force_contact(particle_i_idx, jj) = 0.0d0
                                contact_partner_idx(particle_i_idx, jj) = 0
                                exit
                            end if
                        end do
                    end if

                    curr_j_idx = particle_cell_next(curr_j_idx) ! セル内の次の粒子へ
                end do !! セル内連結リスト走査
            end do ! x方向セルループ
        end do     ! z方向セルループ
    end subroutine pcont_sub

    !> 全粒子間のクーロン力を計算するサブルーチン
    subroutine coulomb_force_sub
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        
        integer :: i, j
        real(8) :: dx, dz, dist, dist_sq, dist_cubed
        real(8) :: force_magnitude, fx, fz
        real(8) :: qi, qj
        
        ! クーロン力が無効化されている場合は何もしない
        if (.not. enable_coulomb_force) return
        
        ! 全粒子ペアについてクーロン力を計算
        do i = 1, num_particles - 1
            qi = charge(i)
            if (abs(qi) < 1.0d-20) cycle ! 電荷がゼロならスキップ
            
            do j = i + 1, num_particles
                qj = charge(j)
                if (abs(qj) < 1.0d-20) cycle ! 電荷がゼロならスキップ
                
                ! 粒子間の距離ベクトルと距離
                dx = x_coord(j) - x_coord(i)
                dz = z_coord(j) - z_coord(i)
                dist_sq = dx*dx + dz*dz
                
                ! ゼロ除算を回避
                if (dist_sq < 1.0d-20) cycle
                
                dist = sqrt(dist_sq)
                dist_cubed = dist * dist_sq
                
                ! クーロン力の大きさ F = k * q1 * q2 / r²
                ! 力のベクトル成分 F_vec = F * (r_vec / |r|) = k * q1 * q2 * r_vec / r³
                force_magnitude = coulomb_constant * qi * qj / dist_cubed
                
                fx = force_magnitude * dx
                fz = force_magnitude * dz
                
                ! 粒子iに力を加算（粒子jへ向かう力）
                x_force_sum(i) = x_force_sum(i) + fx
                z_force_sum(i) = z_force_sum(i) + fz
                
                ! 粒子jに反作用力を加算（粒子iへ向かう力）
                x_force_sum(j) = x_force_sum(j) - fx
                z_force_sum(j) = z_force_sum(j) - fz
            end do
        end do
    end subroutine coulomb_force_sub

    !> 蛙飛び法による粒子の位置と速度を更新するサブルーチン
    subroutine nposit_leapfrog_sub(judge_static, phase)
        ! phase: 0=初期化（v(0)→v(Δt/2)）, 1=位置更新, 2=速度更新
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use simulation_parameters_mod, only: time_step
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        
        integer, intent(in) :: phase
        integer, intent(out) :: judge_static
        
        integer :: i
        real(8) :: sum_abs_disp, avg_abs_disp
        real(8) :: grav, dt
        
        dt = time_step
        grav = GRAVITY_ACCEL
        
        if (phase == 0) then
            ! 初回のみ: v(0) → v(Δt/2) への変換
            do i = 1, num_particles
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * (dt * 0.5d0)
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - grav) * (dt * 0.5d0)
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * (dt * 0.5d0)
            end do
            judge_static = 0
            
        else if (phase == 1) then
            ! フェーズ1: 位置更新のみ
            sum_abs_disp = 0.0d0
            do i = 1, num_particles
                ! 位置更新: x(t+Δt) = x(t) + v(t+Δt/2) * Δt
                x_disp_incr(i) = x_vel(i) * dt
                z_disp_incr(i) = z_vel(i) * dt
                rot_disp_incr(i) = rotation_vel(i) * dt
                
                x_coord(i) = x_coord(i) + x_disp_incr(i)
                z_coord(i) = z_coord(i) + z_disp_incr(i)
                rotation_angle(i) = rotation_angle(i) + rot_disp_incr(i)
                
                sum_abs_disp = sum_abs_disp + abs(x_disp_incr(i)) + abs(z_disp_incr(i))
            end do
            
            ! 静止判定
            if (num_particles > 0) then
                avg_abs_disp = sum_abs_disp / real(num_particles, 8) / 2.0d0
                if (avg_abs_disp < (time_step * time_step * GRAVITY_ACCEL * 1.0d-1)) then
                    judge_static = 1
                else
                    judge_static = 0
                end if
            else
                judge_static = 1
            end if
            
        else if (phase == 2) then
            ! フェーズ2: 速度更新のみ（新しい位置での力を使用）
            do i = 1, num_particles
                ! 速度更新: v(t+3Δt/2) = v(t+Δt/2) + a(t+Δt) * Δt
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * dt
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - grav) * dt
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * dt
            end do
            judge_static = 0
        end if
    end subroutine nposit_leapfrog_sub

    !> 粒子iと粒子/壁jとの間の実際の接触力（法線方向およびせん断方向）を計算するサブルーチン
    subroutine actf_sub(p_i, p_j, contact_slot_idx_for_pi, angle_sin, angle_cos, initial_overlap, en_coeff, et_coeff)
        ! p_i: 主となる粒子のインデックス
        ! p_j: 他方の粒子インデックス (<= num_particles の場合) または壁ID (> num_particles の場合)
        ! contact_slot_idx_for_pi: p_i の接触配列における p_j のスロット
        ! angle_sin, angle_cos: p_i から p_j の中心/接触点への法線ベクトル成分
        ! initial_overlap: 幾何学的な重なり量、接触していれば正
        ! en_coeff, et_coeff: 反発係数 (normal, tangential)
        use simulation_constants_mod, only: ni_max
        use simulation_parameters_mod, only: time_step, reference_overlap
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none

        integer, intent(in) :: p_i, p_j, contact_slot_idx_for_pi
        real(8), intent(in) :: angle_sin, angle_cos, initial_overlap
        real(8), intent(in) :: en_coeff, et_coeff ! 反発係数 (normal, tangential)

        real(8) :: ri_val, rj_val, effective_mass
        real(8) :: kn_normal_stiffness, ks_shear_stiffness ! 法線・せん断バネ定数 Kn, Ks
        real(8) :: damping_coeff_normal, damping_coeff_shear ! 法線・せん断粘性係数 ηn, ηs
        real(8) :: rel_disp_normal_incr, rel_disp_shear_incr ! 法線・せん断方向の相対変位増分 Δun, Δus
        real(8) :: damping_force_normal, damping_force_shear ! 法線・せん断方向の粘性抵抗力 dn, ds
        real(8) :: total_normal_force, total_shear_force     ! 全法線力 fn, 全せん断力 fs
        real(8) :: friction_coeff_current      ! 現在の摩擦係数 μ
        real(8) :: critical_time_step_check    ! 安定性チェック用の時間刻み (ddt)

        real(8) :: x_disp_incr_pi, z_disp_incr_pi, rot_disp_incr_pi ! 粒子iの変位増分
        real(8) :: x_disp_incr_pj, z_disp_incr_pj, rot_disp_incr_pj ! 粒子j(または壁=0)の変位増分
        real(8) :: mass_pi, mass_pj                                 ! 粒子i,jの質量
        real(8) :: r_eff ! 宣言をここに移動
        logical :: is_tracked_pair
        real(8) :: vx_rel, vz_rel, v_rel_n
        real(8) :: tau, alpha, wd, delta_theory
        real(8) :: s1, s2, denom, A_c, B_c
        real(8) :: vx_rel_current, vz_rel_current, v_rel_n_current, v_rel_theory

        ri_val = radius(p_i)
        x_disp_incr_pi = x_disp_incr(p_i)
        z_disp_incr_pi = z_disp_incr(p_i)
        rot_disp_incr_pi = rot_disp_incr(p_i)
        mass_pi = mass(p_i)

        if (p_j <= num_particles) then ! 粒子間
            rj_val = radius(p_j)
            x_disp_incr_pj = x_disp_incr(p_j)
            z_disp_incr_pj = z_disp_incr(p_j)
            rot_disp_incr_pj = rot_disp_incr(p_j)
            mass_pj = mass(p_j)
            if (mass_pi + mass_pj > 1.0d-20) then
                 effective_mass = mass_pi * mass_pj / (mass_pi + mass_pj) ! 等価質量 m_eff = m1*m2/(m1+m2)
            else
                 effective_mass = mass_pi ! フォールバック
            end if
            friction_coeff_current = friction_coeff_particle
            r_eff = ri_val * rj_val / (ri_val + rj_val) ! 等価半径 Reff = ri*rj/(ri+rj)
            
            ! Hertz接触理論に基づく法線剛性(δ_refに基づく) (粒子間)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle / (1.0d0 - poisson_ratio_particle**2) * &
                                 sqrt(reference_overlap)
            ! write(*,*) '粒子間'
            ! write(*,*) 'kn_normal_stiffness=', kn_normal_stiffness

            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        else ! 粒子-壁接触
            rj_val = 0.0d0 ! 壁の半径は実質無限大、または変位にrjは使用しない
            x_disp_incr_pj = 0.0d0 ! 壁は動かないと仮定
            z_disp_incr_pj = 0.0d0
            rot_disp_incr_pj = 0.0d0
            effective_mass = mass_pi   ! 等価質量 m_eff = m1
            friction_coeff_current = friction_coeff_wall
            r_eff = ri_val ! 粒子-壁接触では粒子半径を使用
            
            ! Hertz接触理論に基づく法線剛性(δ_refに基づく) (粒子-壁)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle * young_modulus_wall / &
                                 ((1.0d0-poisson_ratio_particle**2)*young_modulus_wall + &
                                  (1.0d0-poisson_ratio_wall**2)*young_modulus_particle) * &
                                 sqrt(reference_overlap)
            ! write(*,*) '粒子-壁'
            ! write(*,*) 'kn_normal_stiffness=', kn_normal_stiffness
            
            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        end if
        ks_shear_stiffness = kn_normal_stiffness * shear_to_normal_stiffness_ratio
        ! write(*,*) 'ks_shear_stiffness=', ks_shear_stiffness

        ! 完全弾性衝突 (e=1) のための粘性係数設定
        if (en_coeff >= 0.99d0) then
            ! 完全弾性衝突の場合、粘性をほぼゼロに設定
            damping_coeff_normal = 0.0d0
            damping_coeff_shear = 0.0d0
        else
            ! 粘性係数 (反発係数に基づいて計算)
            if (effective_mass > 0.0d0 .and. kn_normal_stiffness > 0.0d0 .and. en_coeff > 1.0d-6) then
                damping_coeff_normal = -2.0d0 * log(en_coeff) * sqrt(effective_mass * kn_normal_stiffness / (log(en_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_normal = 0.0d0
            end if
            if (effective_mass > 0.0d0 .and. ks_shear_stiffness > 0.0d0 .and. et_coeff > 1.0d-6) then
                damping_coeff_shear = -2.0d0 * log(et_coeff) * sqrt(effective_mass * ks_shear_stiffness / (log(et_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_shear = 0.0d0
            end if
        end if

        ! 安定性基準のチェック (元のddt、レイリー時間刻みに関連)
        if (kn_normal_stiffness > 1.0d-12) then
           critical_time_step_check = 0.1d0 * sqrt(effective_mass / kn_normal_stiffness)
           if (critical_time_step_check < time_step .and. critical_time_step_check > 1.0d-12) then 
                write(*,*) '警告: 安定性基準違反 - 推奨時間刻み:', critical_time_step_check, '現在:', time_step
                write(*,*) '  Kn=', kn_normal_stiffness, ' M_eff=', effective_mass, ' 粒子:', p_i, p_j
           end if
        end if

        ! 相対変位増分 (現時間ステップ time_step における)
        ! angle成分は粒子p_iからp_jへの法線ベクトルを定義
        rel_disp_normal_incr = (x_disp_incr_pi - x_disp_incr_pj) * angle_cos + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_sin
        rel_disp_shear_incr = -(x_disp_incr_pi - x_disp_incr_pj) * angle_sin + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_cos + &
                               (ri_val * rot_disp_incr_pi + rj_val * rot_disp_incr_pj)

        ! 弾性力成分の更新 (式3.7, 3.10)
        if (abs(normal_force_contact(p_i, contact_slot_idx_for_pi)) < 1.0d-8) then ! 新規接触
            ! 新規接触の場合、弾性力を重なり量に基づいて初期化
            normal_force_contact(p_i, contact_slot_idx_for_pi) = kn_normal_stiffness * initial_overlap ! 弾性項抜き出す部分
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0 ! せん断力は初期化時にゼロ
        else
            ! 既存接触の場合、増分で更新
            normal_force_contact(p_i, contact_slot_idx_for_pi) = normal_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                 kn_normal_stiffness * rel_disp_normal_incr 
            shear_force_contact(p_i, contact_slot_idx_for_pi) = shear_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                ks_shear_stiffness * rel_disp_shear_incr
        end if
        
        ! 粘性抵抗力成分の計算 (式3.6, 3.9)
        if (time_step > 1.0d-20) then
             damping_force_normal = damping_coeff_normal * rel_disp_normal_incr / time_step 
             damping_force_shear = damping_coeff_shear * rel_disp_shear_incr / time_step
        else
             damping_force_normal = 0.0d0
             damping_force_shear  = 0.0d0
        end if

        ! 引張力のチェック (粒子が引き離される場合) - 付着力は考慮しない
        if (normal_force_contact(p_i, contact_slot_idx_for_pi) < 0.0d0) then
            normal_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            damping_force_normal = 0.0d0
            damping_force_shear = 0.0d0
            contact_partner_idx(p_i, contact_slot_idx_for_pi) = 0 ! 引張なら接触を切る
            return ! 引き離される場合は力なし
        end if

        ! クーロンの摩擦法則を適用 (式3.11)
        if (abs(shear_force_contact(p_i, contact_slot_idx_for_pi)) > &
            friction_coeff_current * normal_force_contact(p_i, contact_slot_idx_for_pi)) then
            shear_force_contact(p_i, contact_slot_idx_for_pi) = friction_coeff_current * &
                normal_force_contact(p_i, contact_slot_idx_for_pi) * &
                sign(1.0d0, shear_force_contact(p_i, contact_slot_idx_for_pi))
            damping_force_shear = 0.0d0 ! 滑りが発生している場合はせん断粘性なし
        end if

        ! 粘性を含む合計の力 (式3.8, 3.12)
        total_normal_force = normal_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_normal 
        total_shear_force = shear_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_shear

        ! 粒子p_iに力を適用 (式3.13)
        ! 法線力は中心を結ぶ線に沿って作用 (angle_cos, angle_sin で定義される iからjへの方向)
        ! せん断力はそれに垂直。
        x_force_sum(p_i) = x_force_sum(p_i) - total_normal_force * angle_cos + total_shear_force * angle_sin
        z_force_sum(p_i) = z_force_sum(p_i) - total_normal_force * angle_sin - total_shear_force * angle_cos    
        moment_sum(p_i) = moment_sum(p_i) - ri_val * total_shear_force

        ! 粒子p_jに反作用力を適用 (相手が粒子の場合)
        if (p_j <= num_particles .and. contact_slot_idx_for_pi <= 10) then ! 元の jk < 10 は粒子間接触のチェック
            x_force_sum(p_j) = x_force_sum(p_j) + total_normal_force * angle_cos - total_shear_force * angle_sin
            z_force_sum(p_j) = z_force_sum(p_j) + total_normal_force * angle_sin + total_shear_force * angle_cos
            moment_sum(p_j) = moment_sum(p_j) - rj_val * total_shear_force ! p_jに対するせん断力は同じ大きさ、逆向きの回転効果
            
            ! デバッグ出力 (検証モードのみ)
            ! if (validation_mode .and. abs(total_normal_force) > 1.0d-6) then
            !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5)') 'Force p',p_i,' -> p',p_j,': Fn=',total_normal_force,', Fs=',total_shear_force
            !     write(*,'(A,ES12.5,A,ES12.5)') '  overlap=',initial_overlap,', dist=',sqrt((x_coord(p_i)-x_coord(p_j))**2 + (z_coord(p_i)-z_coord(p_j))**2)
            ! end if
        end if
        
        ! 現在のオーバーラップを記録（次ステップで接触開始検出に使用）
        previous_overlap(p_i, contact_slot_idx_for_pi) = initial_overlap
    end subroutine actf_sub

    !> グラフィック用データを出力するサブルーチン
    subroutine gfout_sub(iter_step, time_val, rmax_val)
        use simulation_constants_mod, only: nj_max, GRAVITY_ACCEL
        use simulation_parameters_mod, only: container_width, container_height, time_step
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: iter_step    ! 現在のイテレーションステップ
        real(8), intent(in) :: time_val, rmax_val ! 現在時刻、最大粒子半径
        integer :: i,j
        real(8) :: dt, grav
        real(8), allocatable :: vx_out(:), vz_out(:), rotation_vel_out(:)

        if (iter_step == 1) then
            open(unit=10, file='data/graph11.d', status='replace', action='write')
            open(unit=11, file='data/graph21.d', status='replace', action='write')
        end if

        write(10,*) num_particles, time_val, container_width, container_height, rmax_val
        if (num_particles > 0) then
            ! 出力用補正速度の準備（蛙飛び法のみ 0.5*dt*加速度で補正）
            dt = time_step
            grav = GRAVITY_ACCEL

            allocate(vx_out(num_particles), vz_out(num_particles), rotation_vel_out(num_particles))
            ! 出力用補正速度の計算
            do i = 1, num_particles
                vx_out(i) = x_vel(i) - 0.5d0 * dt * (x_force_sum(i) / mass(i))
                vz_out(i) = z_vel(i) - 0.5d0 * dt * (z_force_sum(i) / mass(i) - grav)
                rotation_vel_out(i) = rotation_vel(i) - 0.5d0 * dt * (moment_sum(i) / moment_inertia(i))
            end do

            write(10,'(1000(ES12.5,1X,ES12.5,1X,ES12.5,2X))') (sngl(x_coord(i)), sngl(z_coord(i)), sngl(radius(i)), i=1,num_particles)
            write(10,'(1000(ES12.5,1X,ES12.5,1X,ES12.5,2X))') (sngl(vx_out(i)), sngl(vz_out(i)), sngl(rotation_vel_out(i)), i=1,num_particles)
            write(10,'(1000(ES12.5,2X))') (sngl(rotation_angle(i)), i=1,num_particles)
        end if
        
        ! 接触力の出力 (オプション、graph21.dより)
        write(11,*) 'Time: ', time_val 
        if (num_particles > 0) then
            do i = 1, num_particles
                 write(11,'(A,I5,A,I5)') 'Particle: ', i, ' NumContacts: ', count(contact_partner_idx(i,1:nj_max) > 0)
                 write(11,'(2X,A,13(ES10.3,1X))') 'ShearF: ', (shear_force_contact(i,j), j=1,nj_max)
                 write(11,'(2X,A,13(ES10.3,1X))') 'NormalF:', (normal_force_contact(i,j), j=1,nj_max)
                 write(11,'(2X,A,13(I5,2X))')    'Partner:', (contact_partner_idx(i,j), j=1,nj_max)
            end do
        end if
    end subroutine gfout_sub

    !> バックアップデータを出力するサブルーチン
    subroutine bfout_sub
        use simulation_constants_mod
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod
        implicit none
        integer :: i, j
        real(8) :: rmax_dummy_val ! 元のbfoutはrmaxを必要とするが、メインの呼び出しからは渡されない。
                                  ! リストアに不可欠でないか、粒子半径から取得すると仮定。
        
        if (num_particles > 0) then
           rmax_dummy_val = maxval(radius(1:num_particles))
        else
           rmax_dummy_val = 0.0d0
        end if

        open(unit=13, file='data/backl.d', status='replace', action='write')

        write(13,*) num_particles, cells_x_dir, cells_z_dir, particle_gen_layers
        write(13,*) rmax_dummy_val, 0.0d0, container_width, container_height, cell_size, time_step ! current_timeではなく初期t=0を保存すると仮定
        write(13,*) particle_density, friction_coeff_particle, friction_coeff_wall, GRAVITY_ACCEL, PI_VAL
        write(13,*) young_modulus_particle, young_modulus_wall, poisson_ratio_particle, poisson_ratio_wall, shear_to_normal_stiffness_ratio
        
        if (num_particles > 0) then
            write(13,*) (mass(i), moment_inertia(i), i=1,num_particles)
            write(13,*) (x_coord(i), z_coord(i), radius(i), charge(i), i=1,num_particles)
            write(13,*) (x_disp_incr(i), z_disp_incr(i), rot_disp_incr(i), i=1,num_particles) ! u,v,f (dpm)
            write(13,*) (x_vel(i), z_vel(i), rotation_vel(i), i=1,num_particles)            ! u0,v0,f0
            do i = 1, num_particles
                write(13,*) (shear_force_contact(i,j), normal_force_contact(i,j), j=1,nj_max)
                write(13,*) (contact_partner_idx(i,j), j=1,nj_max)
            end do
        end if
        close(13)
    end subroutine bfout_sub

    !> 擬似乱数を生成するサブルーチン
    subroutine custom_random(seed_io, random_val_out)
        implicit none
        integer, intent(inout) :: seed_io          ! ジェネレータの現在の状態を保持
        real(8), intent(out)   :: random_val_out   ! [0,1) の一様乱数

        ! 元のコードの定数とロジックを可能な限り再現
        seed_io = seed_io * 65539 
        if (seed_io < 0) then
             seed_io = (seed_io + 2147483647) + 1 
        end if
        random_val_out = dble(seed_io) * 0.4656613d-9 ! 元の正規化定数 (1.0 / 2147483648.0)

    end subroutine custom_random

end program two_dimensional_pem