module coulomb_force_mod
  implicit none
  real(8), parameter :: eps0 = 8.854187817d-12
contains

! *****************************************************
! 粒子 i に働く静電力を計算
! ・背景電場 E(x,z)
! ・粒子間クーロン力
! ・画像電荷力（誘電体壁）
! *****************************************************
subroutine compute_coulomb_force(i, np, x, z, q, Ex_p, Ez_p, eps_d, wall_y, Fx, Fz)
  implicit none
  integer, intent(in) :: i, np
  real(8), intent(in) :: x(np), z(np), q(np)
  real(8), intent(in) :: Ex_p, Ez_p     ! 背景電場
  real(8), intent(in) :: eps_d          ! 壁の誘電率（アクリルなど）
  real(8), intent(in) :: wall_z         ! 壁の位置（例：z=0 の下壁）
  real(8), intent(out):: Fx, Fz

  integer :: k
  real(8) :: rx, rz, r2, r3
  real(8) :: Fc_x, Fc_z
  real(8) :: F_img, d, n_z
  real(8) :: coeff

  ! 初期化
  Fx = 0.d0
  Fz = 0.d0

  ! --------------------------------------
  ! (1) 背景電場からのクーロン力
  ! --------------------------------------
  Fx = Fx + q(i) * Ex_p
  Fz = Fz + q(i) * Ez_p


  ! --------------------------------------
  ! (2) 粒子間クーロン相互作用
  ! --------------------------------------
  do k = 1, np
     if (k == i) cycle

     rx = x(i) - x(k)
     rz = z(i) - z(k)
     r2 = rx*rx + rz*rz
     r3 = r2 * sqrt(r2)

     coeff = (q(i)*q(k)) / (4.d0 * 3.14159265358979d0 * eps0 * r3)

     Fx = Fx + coeff * rx
     Fz = Fz + coeff * rz
  end do


  ! --------------------------------------
  ! (3) 電気影像力
  ! --------------------------------------
  ! 壁面の法線方向（ここでは z+ 方向）
  n_z = 1.d0

  ! 粒子中心から壁までの距離
  d = z(i) - wall_z
  if (d > 0.d0) then
     F_img = - ( q(i)*q(i) / (4.d0*3.14159265358979d0*eps0) ) * &
             ( (eps_d - eps0)/(eps_d + eps0) ) * (1.d0/(2.d0*d))**2

     Fz = Fz + F_img * n_z
  end if

end subroutine compute_coulomb_force

end module coulomb_force_mod
