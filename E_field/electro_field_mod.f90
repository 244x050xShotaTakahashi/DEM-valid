! *****************************************************
! 背景電場モジュール
! ・背景電場データの読み込み
! ・双線形補間で電場を求める
! *****************************************************
module electro_field_mod
  implicit none
  integer, parameter :: nx = 101
  integer, parameter :: ny = 101
  real(8), allocatable :: xgrid(:), ygrid(:)
  real(8), allocatable :: Ex_grid(:,:), Ey_grid(:,:)

contains
  subroutine load_field()
    character(len=*), parameter :: filename = "./E_field/data/electro_field_2D.dat"
    integer :: i, j

    open(10, file=filename, status='old', action='read')
    read(10, *) nx, ny

    allocate(xgrid(nx), ygrid(ny))
    allocate(Ex_grid(ny,nx), Ey_grid(ny,nx))

    do i = 1, ny
       do j = 1, nx
          read(10, *) xgrid(j), ygrid(i), Ex_grid(i,j), Ey_grid(i,j)
       end do
    end do
    close(10)
  end subroutine load_field


  ! -----------------------------
  ! 双線形補間で電場を求める
  ! -----------------------------
  subroutine interpolate_field(xp, yp, Ex_p, Ey_p)
    real(8), intent(in) :: xp, yp
    real(8), intent(out) :: Ex_p, Ey_p
    integer :: i, j
    real(8) :: wx, wy

    ! --- x のセルを特定 ---
    if (xp <= xgrid(1)) then
       j = 1
    else if (xp >= xgrid(nx)) then
       j = nx - 1
    else
       do j = 1, nx-1
          if (xgrid(j) <= xp .and. xp < xgrid(j+1)) exit
       end do
    end if

    ! --- y のセルを特定 ---
    if (yp <= ygrid(1)) then
       i = 1
    else if (yp >= ygrid(ny)) then
       i = ny - 1
    else
       do i = 1, ny-1
          if (ygrid(i) <= yp .and. yp < ygrid(i+1)) exit
       end do
    end if

    ! --- セル内の重み ---
    wx = (xp - xgrid(j)) / (xgrid(j+1) - xgrid(j))
    wy = (yp - ygrid(i)) / (ygrid(i+1) - ygrid(i))

    ! --- 双線形補間 ---
    Ex_p = (1-wx)*(1-wy)*Ex_grid(i,j) + wx*(1-wy)*Ex_grid(i,j+1) + &
           (1-wx)*wy*Ex_grid(i+1,j) + wx*wy*Ex_grid(i+1,j+1)

    Ey_p = (1-wx)*(1-wy)*Ey_grid(i,j) + wx*(1-wy)*Ey_grid(i,j+1) + &
           (1-wx)*wy*Ey_grid(i+1,j) + wx*wy*Ey_grid(i+1,j+1)

end module electro_field_mod
