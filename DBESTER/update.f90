subroutine update

use defvar
!$ use omp_lib

implicit none

	real*8 :: vminusx,vminusy,vminusz,ex,ey,ez,bx,by,bz,alpha,beta,ganma
	integer :: realx,realy,realz,realx1,realy1,realz1
        integer :: istep, m

	do istep=1,nstep

		do m=1,np(1) 
		
                  !linear interpolation for E and B
                  realx = floor(pbuf(m)%x)
                  realy = floor(pbuf(m)%y)
                  realz = floor(pbuf(m)%z)
                  realx1 = realx + 1
                  realy1 = realy + 1
                  realz1 = realz + 1
                  alpha = pbuf(m)%x-realx
                  beta  = pbuf(m)%y-realy
                  ganma = pbuf(m)%z-realz

                  ex = pex(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pex(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pex(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pex(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
                      +pex(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pex(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pex(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pex(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)
                  ey = pey(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pey(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pey(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pey(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
                      +pey(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pey(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pey(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pey(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)
                  ez = pez(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pez(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pez(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pez(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
                      +pez(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pez(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pez(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pez(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)
                  bx = pbx(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pbx(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pbx(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pbx(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
                      +pbx(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pbx(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pbx(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pbx(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)
                  by = pby(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pby(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pby(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pby(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pby(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pby(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pby(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)
                  bz = pbz(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
                      +pbz(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
                      +pbz(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
                      +pbz(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
                      +pbz(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
                      +pbz(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
                      +pbz(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
                      +pbz(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)

      print*, "step:",istep
      print*, "pbuf", &
     & pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%vx,pbuf(m)%vy,pbuf(m)%vz
      print*, "realxyz01", &
     & realx,realy,realz,realx1,realy1,realz1
      print*, "abg", &
     alpha,beta,ganma
      print*, "ebxyz", &
     ex,ey,ez,bx,by,bz
      print*, "qm(1)", &
     qm(1)
      print*, "b to x", qm(1)*(pbuf(m)%vy*bz-pbuf(m)%vz*by)

		!velocity
			pbuf(m)%vx = pbuf(m)%vx + qm(1)*ex
			pbuf(m)%vy = pbuf(m)%vy + qm(1)*ey
			pbuf(m)%vz = pbuf(m)%vz + qm(1)*ez

			vminusx = pbuf(m)%vx + qm(1)*(pbuf(m)%vy*bz-pbuf(m)%vz*by)
			vminusy = pbuf(m)%vy + qm(1)*(pbuf(m)%vz*bx-pbuf(m)%vx*bz)
			vminusz = pbuf(m)%vz + qm(1)*(pbuf(m)%vx*by-pbuf(m)%vy*bx)

                        pbuf(m)%vx = pbuf(m)%vx + 2*qm(1)/(1+qm(1)*qm(1)*(bx*bx+by*by+bz*bz))*(vminusy*bz-vminusz*by)
                        pbuf(m)%vy = pbuf(m)%vy + 2*qm(1)/(1+qm(1)*qm(1)*(bx*bx+by*by+bz*bz))*(vminusz*bx-vminusx*bz)
                        pbuf(m)%vz = pbuf(m)%vz + 2*qm(1)/(1+qm(1)*qm(1)*(bx*bx+by*by+bz*bz))*(vminusx*by-vminusy*bx)

			pbuf(m)%vx = pbuf(m)%vx + qm(1)*ex
			pbuf(m)%vy = pbuf(m)%vy + qm(1)*ey
			pbuf(m)%vz = pbuf(m)%vz + qm(1)*ez

		!position
			pbuf(m)%x =pbuf(m)%x+pbuf(m)%vx*2
			pbuf(m)%y =pbuf(m)%y+pbuf(m)%vy*2
			pbuf(m)%z =pbuf(m)%z+pbuf(m)%vz*2

      if (pbuf(m)%x<0) then
         pbuf(m)%x = dims(1)+pbuf(m)%x
      end if
      if (pbuf(m)%y<0) then
         pbuf(m)%y = dims(2)+pbuf(m)%y
      end if
      if (pbuf(m)%z<0) then
         pbuf(m)%z = dims(3)+pbuf(m)%z
      end if
      if (pbuf(m)%x>dims(1)) then
         pbuf(m)%x = pbuf(m)%x-dims(1)
      end if
      if (pbuf(m)%y>dims(2)) then
         pbuf(m)%y = pbuf(m)%y-dims(2)
      end if
      if (pbuf(m)%z>dims(3)) then
         pbuf(m)%z = pbuf(m)%z-dims(3)
      end if

			pdatax(m,istep)=pbuf(m)%x
			pdatay(m,istep)=pbuf(m)%y
			pdataz(m,istep)=pbuf(m)%z
			pdatavx(m,istep)=pbuf(m)%vx
			pdatavy(m,istep)=pbuf(m)%vy
			pdatavz(m,istep)=pbuf(m)%vz

		end do
 
		do m=np(1)+1,np(1)+np(2) 
			!ion

			!linear interpolation for E and B
			realx = floor(pbuf(m)%x)
			realy = floor(pbuf(m)%y)
			realz = floor(pbuf(m)%z)
                        realx1 = realx + 1
                        realy1 = realy + 1
                        realz1 = realz + 1
			alpha = pbuf(m)%x-realx
			beta  = pbuf(m)%y-realy
			ganma = pbuf(m)%z-realz

			ex = pex(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pex(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pex(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pex(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pex(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pex(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pex(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pex(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)

			ey = pey(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pey(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pey(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pey(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pey(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pey(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pey(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pey(realx1,realy1,realz1)*(alpha)*(beta)*(ganma) 

			ez = pez(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pez(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pez(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pez(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pez(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pez(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pez(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pez(realx1,realy1,realz1)*(alpha)*(beta)*(ganma) 

			bx = pbx(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pbx(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pbx(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pbx(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pbx(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pbx(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pbx(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pbx(realx1,realy1,realz1)*(alpha)*(beta)*(ganma) 

			by = pby(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pby(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pby(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pby(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pby(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pby(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pby(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pby(realx1,realy1,realz1)*(alpha)*(beta)*(ganma) 

			bz = pbz(realx,realy,realz)*(1-alpha)*(1-beta)*(1-ganma) &
				+pbz(realx1,realy,realz)*(alpha)*(1-beta)*(1-ganma) &
				+pbz(realx,realy1,realz)*(1-alpha)*(beta)*(1-ganma) &
				+pbz(realx1,realy1,realz)*(alpha)*(beta)*(1-ganma) &
				+pbz(realx,realy,realz1)*(1-alpha)*(1-beta)*(ganma) &
				+pbz(realx1,realy,realz1)*(alpha)*(1-beta)*(ganma) &
				+pbz(realx,realy1,realz1)*(1-alpha)*(beta)*(ganma) &
				+pbz(realx1,realy1,realz1)*(alpha)*(beta)*(ganma)

		!velocity
			pbuf(m)%vx = pbuf(m)%vx + qm(2)*ex
			pbuf(m)%vy = pbuf(m)%vy + qm(2)*ey
			pbuf(m)%vz = pbuf(m)%vz + qm(2)*ez

                        vminusx = pbuf(m)%vx + qm(2)*(pbuf(m)%vy*bz-pbuf(m)%vz*by)
                        vminusy = pbuf(m)%vy + qm(2)*(pbuf(m)%vz*bx-pbuf(m)%vx*bz)
                        vminusz = pbuf(m)%vz + qm(2)*(pbuf(m)%vx*by-pbuf(m)%vy*bx)

                        pbuf(m)%vx = pbuf(m)%vx + 2*qm(2)/(1+qm(2)*qm(2)*(bx*bx+by*by+bz*bz))*(vminusy*bz-vminusz*by)
                        pbuf(m)%vy = pbuf(m)%vy + 2*qm(2)/(1+qm(2)*qm(2)*(bx*bx+by*by+bz*bz))*(vminusz*bx-vminusx*bz)
                        pbuf(m)%vz = pbuf(m)%vz + 2*qm(2)/(1+qm(2)*qm(2)*(bx*bx+by*by+bz*bz))*(vminusx*by-vminusy*bx)

			pbuf(m)%vx = pbuf(m)%vx + qm(2)*ex
			pbuf(m)%vy = pbuf(m)%vy + qm(2)*ey
			pbuf(m)%vz = pbuf(m)%vz + qm(2)*ez

		!position		
			pbuf(m)%x =pbuf(m)%x+pbuf(m)%vx*2
			pbuf(m)%y =pbuf(m)%y+pbuf(m)%vy*2
			pbuf(m)%z =pbuf(m)%z+pbuf(m)%vz*2

      if (pbuf(m)%x<0) then
         pbuf(m)%x = dims(1)+pbuf(m)%x
      end if
      if (pbuf(m)%y<0) then
         pbuf(m)%y = dims(2)+pbuf(m)%y
      end if
      if (pbuf(m)%z<0) then
         pbuf(m)%z = dims(3)+pbuf(m)%z
      end if
      if (pbuf(m)%x>dims(1)) then
         pbuf(m)%x = pbuf(m)%x-dims(1)
      end if
      if (pbuf(m)%y>dims(2)) then
         pbuf(m)%y = pbuf(m)%y-dims(2)
      end if
      if (pbuf(m)%z>dims(3)) then
         pbuf(m)%z = pbuf(m)%z-dims(3)
      end if

			pdatax(m,istep)=pbuf(m)%x
			pdatay(m,istep)=pbuf(m)%y
			pdataz(m,istep)=pbuf(m)%z
			pdatavx(m,istep)=pbuf(m)%vx
			pdatavy(m,istep)=pbuf(m)%vy
			pdatavz(m,istep)=pbuf(m)%vz

		end do

	end do

	end subroutine update
