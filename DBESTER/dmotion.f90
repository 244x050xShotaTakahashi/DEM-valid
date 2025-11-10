program dmotion

implicit none

type dust
   real(kind=8) :: x,y,z,vx,vy,vz
   integer(kind=4) :: nid
end type dust

integer(kind=4) :: i,j,k,i1,j1,k1,m,n,sind,eind,ns,ne,N,stride
integer(kind=4) :: ACC_R, ACC_W, ACC_RW, ACC_C
integer(kind=4) :: EX,EY,EZ
integer(kind=4) :: idt1, idt2, idt3, idt4, idt5, idt6, ids1, ids2, ids3, ist, ist1, ist2,iind
integer(kind=8) :: totdim(3),dustp(3)
integer(kind=8) :: lbhole,ubhole,xholec,yholec,rhole
character(len=40) :: fnorgx, fnorgy, fnorgz, fnorgx2, fnorgy2, fnorgz2
character(len=16) :: sdsnamex, sdsnamey, sdsnamez, sdsnamex2, sdsnamey2, sdsnamez2
character(len=20) :: str(16)
real(kind=8),allocatable :: eb(:,:,:,:),wrk(:,:,:,:)
real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2,xlocal,ylocal,zlocal
real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8,dt,xl,yl,zl,xu,yu,zu,eex,eey,eez
real(kind=8) :: vxt,vyt,vzt,q,qm,r,g,conv,pi
type(dust) :: pbuf(4)

parameter(ACC_R=1, ACC_W=2, ACC_RW=3, ACC_C=4)
parameter(EX=1,EY=2,EZ=3)
parameter(pi=3.14159265358979324)

call getarg(1,str(1))
call getarg(2,str(2))
call getarg(3,str(3))
call getarg(4,str(4))
call getarg(5,str(5))
call getarg(6,str(6))
call getarg(7,str(7))
call getarg(8,str(8))
call getarg(9,str(9))
call getarg(10,str(10))
call getarg(11,str(11))
call getarg(12,str(12))
call getarg(13,str(13))
call getarg(14,str(14))
call getarg(15,str(15))
call getarg(16,str(16))
call getarg(17,str(17))

read(str(1),*) totdim(1)
read(str(2),*) totdim(2)
read(str(3),*) totdim(3)
read(str(4),*) dustp(1)
read(str(5),*) dustp(2)
read(str(6),*) dustp(3)
read(str(7),*) r
read(str(8),*) lbhole
read(str(9),*) ubhole
read(str(10),*) xholec
read(str(11),*) yholec
read(str(12),*) rhole
read(str(13),*) stride
read(str(14),*) N
read(str(15),*) conv
read(str(16),*) sind
read(str(17),*) eind

allocate(eb(3,-1:totdim(1),-1:totdim(2),-1:totdim(3)))
allocate(wrk(3,0:totdim(1),0:totdim(2),0:totdim(3)))

q = 1.602176565e-19
dt = 1.0E-2
ns = 1; ne = 1;
xl = 0; yl = 0; zl = 0;
xu = totdim(1); yu = totdim(2); zu = totdim(3);
qm = q / ((4/3)*pi*r**3*2.5) * dt
g = 1.6 * dt

!------------------- initialize dust data
do m=ns,ne
   pbuf(m)%x = dustp(1); pbuf(m)%y = dustp(2); pbuf(m)%z = dustp(3);
   pbuf(m)%vz = 1.0
end do

!------------------- import E-data
!write(fnorgx,'(a,i3.3,a,i3.3,a)') 'ex00000000_id:', sind, '-', eind ,'.h5'
!write(fnorgy,'(a,i3.3,a,i3.3,a)') 'ey00000000_id:', sind, '-', eind ,'.h5'
!write(fnorgz,'(a,i3.3,a,i3.3,a)') 'ez00000000_id:', sind, '-', eind ,'.h5'
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'ex', sind, '-', eind ,'a.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'ey', sind, '-', eind ,'a.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'ez', sind, '-', eind ,'a.h5'

!write(sdsnamex,'(a,i4.4)') fnorgx(1:2),eind
!write(sdsnamey,'(a,i4.4)') fnorgy(1:2),eind
!write(sdsnamez,'(a,i4.4)') fnorgz(1:2),eind
write(sdsnamex,'(a,i4.4)') fnorgx(1:2),sind
write(sdsnamey,'(a,i4.4)') fnorgy(1:2),sind
write(sdsnamez,'(a,i4.4)') fnorgz(1:2),sind

call hdfinit()

call hdfopen(fnorgx,ids1,ACC_R)
call hdfopen(fnorgy,ids2,ACC_R)
call hdfopen(fnorgz,ids3,ACC_R)

call read3d(ids1,sdsnamex,totdim+1,eb(EX,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids2,sdsnamey,totdim+1,eb(EY,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids3,sdsnamez,totdim+1,eb(EZ,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)

call hdfclose(ids1,ist)
call hdfclose(ids2,ist)
call hdfclose(ids3,ist)

!-------------------- S > R
eb(:,:,:,:) = eb(:,:,:,:)/conv

!-------------------- relocation of pe-fields
do k=0,zu-zl
   do j=0,yu-yl
      do i=0,xu-xl
         wrk(EX,i,j,k) = 0.5d0*(eb(EX,i,j,k) + eb(EX,i-1,j,k))
         wrk(EY,i,j,k) = 0.5d0*(eb(EY,i,j,k) + eb(EY,i,j-1,k))
         wrk(EZ,i,j,k) = 0.5d0*(eb(EZ,i,j,k) + eb(EZ,i,j,k-1))
      end do
   end do
end do

!-------------------- open output file
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'dpx00000000_id:', 0, '-', N/stride ,'.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'dpy00000000_id:', 0, '-', N/stride ,'.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'dpz00000000_id:', 0, '-', N/stride ,'.h5'

write(fnorgx2,'(a,i3.3,a,i3.3,a)') 'dvx00000000_id:', 0, '-', N/stride ,'.h5'
write(fnorgy2,'(a,i3.3,a,i3.3,a)') 'dvy00000000_id:', 0, '-', N/stride ,'.h5'
write(fnorgz2,'(a,i3.3,a,i3.3,a)') 'dvz00000000_id:', 0, '-', N/stride ,'.h5'

call hdfopen(fnorgx,idt1,ACC_C)
call hdfopen(fnorgy,idt2,ACC_C)
call hdfopen(fnorgz,idt3,ACC_C)

call hdfopen(fnorgx2,idt4,ACC_C)
call hdfopen(fnorgy2,idt5,ACC_C)
call hdfopen(fnorgz2,idt6,ACC_C)

open(11,file='output.dat', status='replace')

!-------------------- update x,y,z,vx,vy,vz
do n = 0,N
   do m=ns,ne
      if(pbuf(m)%nid.eq.-1) cycle
      
      xlocal = pbuf(m)%x
      ylocal = pbuf(m)%y
      zlocal = pbuf(m)%z 
      
      i = int(xlocal)
      j = int(ylocal)
      k = int(zlocal)
      
      i1 = i + 1
      j1 = j + 1
      k1 = k + 1
   
      x1 = xlocal - i
      y1 = ylocal - j
      z1 = zlocal - k
   
      xy1 = x1*y1
      xz1 = x1*z1
      yz1 = y1*z1
      z2 = 1.0d0 - z1
      xz2 = x1*z2
      yz2 = y1*z2
      v3 = xy1*z1
      v2 = xz1 - v3
      v4 = yz1 - v3
      v1 = z1 - xz1 - v4
      v7 = xy1*z2
      v6 = xz2 - v7
      v8 = yz2 - v7
      v5 = z2 - xz2 - v8
      
      !----------- e-fields interpolation
      eex = wrk(EX,i ,j ,k1)*v1 + wrk(EX,i1,j ,k1)*v2 &
     &    + wrk(EX,i1,j1,k1)*v3 + wrk(EX,i ,j1,k1)*v4 &
     &    + wrk(EX,i ,j ,k )*v5 + wrk(EX,i1,j ,k )*v6 &
     &    + wrk(EX,i1,j1,k )*v7 + wrk(EX,i ,j1,k )*v8
      eey = wrk(EY,i ,j ,k1)*v1 + wrk(EY,i1,j ,k1)*v2 &
     &    + wrk(EY,i1,j1,k1)*v3 + wrk(EY,i ,j1,k1)*v4 &
     &    + wrk(EY,i ,j ,k )*v5 + wrk(EY,i1,j ,k )*v6 &
     &    + wrk(EY,i1,j1,k )*v7 + wrk(EY,i ,j1,k )*v8
      eez = wrk(EZ,i ,j ,k1)*v1 + wrk(EZ,i1,j ,k1)*v2 &
     &    + wrk(EZ,i1,j1,k1)*v3 + wrk(EZ,i ,j1,k1)*v4 &
     &    + wrk(EZ,i ,j ,k )*v5 + wrk(EZ,i1,j ,k )*v6 &
     &    + wrk(EZ,i1,j1,k )*v7 + wrk(EZ,i ,j1,k )*v8

      !----------- update particle velocities
      pbuf(m)%vx = pbuf(m)%vx + eex*qm
      pbuf(m)%vy = pbuf(m)%vy + eey*qm
      pbuf(m)%vz = pbuf(m)%vz + eez*qm - g
 
      !----------- update particle positions
      pbuf(m)%x = pbuf(m)%x + pbuf(m)%vx*dt
      pbuf(m)%y = pbuf(m)%y + pbuf(m)%vy*dt
      pbuf(m)%z = pbuf(m)%z + pbuf(m)%vz*dt
      
      !----------- check particle positions
      if(pbuf(m)%x<xl.or.xu<pbuf(m)%x)     pbuf(m)%nid = -1
      if(pbuf(m)%y<yl.or.yu<pbuf(m)%y)     pbuf(m)%nid = -1
      if(pbuf(m)%z<lbhole.or.zu<pbuf(m)%z) pbuf(m)%nid = -1
      if(pbuf(m)%z<ubhole.and.rhole<((pbuf(m)%x-xholec)**2 &
     & +(pbuf(m)%y-yholec)**2)**0.5) pbuf(m)%nid = -1
      
   end do
   
   !----------- output particle positionos
   if(mod(n,stride)==0) then
      write(sdsnamex,'(a,i4.4)') fnorgx(1:3),n/stride
      write(sdsnamey,'(a,i4.4)') fnorgy(1:3),n/stride
      write(sdsnamez,'(a,i4.4)') fnorgz(1:3),n/stride
         
      write(sdsnamex2,'(a,i4.4)') fnorgx2(1:3),n/stride
      write(sdsnamey2,'(a,i4.4)') fnorgy2(1:3),n/stride
      write(sdsnamez2,'(a,i4.4)') fnorgz2(1:3),n/stride
         
      call wrt1d(idt1,sdsnamex,ne,pbuf(ns:ne)%x,ist1,ist2)
      call wrt1d(idt2,sdsnamey,ne,pbuf(ns:ne)%y,ist1,ist2)
      call wrt1d(idt3,sdsnamez,ne,pbuf(ns:ne)%z,ist1,ist2)
     
      call wrt1d(idt4,sdsnamex2,ne,pbuf(ns:ne)%vx,ist1,ist2)
      call wrt1d(idt5,sdsnamey2,ne,pbuf(ns:ne)%vy,ist1,ist2)
      call wrt1d(idt6,sdsnamez2,ne,pbuf(ns:ne)%vz,ist1,ist2)

      write(11,'(i5,3f12.6)') n, pbuf(1)%x, pbuf(1)%y, pbuf(1)%z
   end if

   write(6,*) n

end do

call hdfclose(idt1,ist)
call hdfclose(idt2,ist)
call hdfclose(idt3,ist) 

call hdfclose(idt4,ist)
call hdfclose(idt5,ist)
call hdfclose(idt6,ist)

call hdffinalize()

close(11)

end program dmotion
