program dmotion2

implicit none

type dust
   real(kind=8) :: x,y,z,vx,vy,vz,q
   integer(kind=4) :: nid
end type dust

integer(kind=4) :: i,j,k,i1,j1,k1,m,n,sind,eind,ns,ne,N,stride,dnum,c
integer(kind=4) :: ACC_R, ACC_W, ACC_RW, ACC_C
integer(kind=4) :: X,Y,EX,EY,EZ,J1X,J1Y,J1Z,J2X,J2Y,J2Z,J3X,J3Y,J3Z
integer(kind=4) :: idt1, idt2, idt3, idt4, idt5, idt6, idt7, ids1, ids2, ids3, ist, ist1, ist2,iind
integer(kind=8) :: totdim(3),dustp(3)
integer(kind=8) :: lbhole,ubhole,xholec,yholec,rhole,xholec2,yholec2
character(len=40) :: fnorgx, fnorgy, fnorgz, fnorgx2, fnorgy2, fnorgz2, fnorg, fvtk
character(len=16) :: sdsnamex, sdsnamey, sdsnamez, sdsnamex2, sdsnamey2, sdsnamez2, sdsname
character(len=20) :: str(30)
real(kind=8),allocatable :: ej(:,:,:,:),wrk(:,:,:,:),random(:,:)
real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2,xlocal,ylocal,zlocal
real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8,dt,xl,yl,zl,xu,yu,zu
real(kind=8) :: eex,eey,eez,jj1,jj1x,jj1y,jj1z,jj2,jj2x,jj2y,jj2z,jj3,jj3x,jj3y,jj3z,jph,jjph,jplus,jminus
real(kind=8) :: vxt,vyt,vzt,q,qm,r,v0,g,conve,convj,pi,dm,dv,ds,dq,den,theta,distance,distance2,e,a,bar
type(dust),allocatable :: pbuf(:)

parameter(ACC_R=1, ACC_W=2, ACC_RW=3, ACC_C=4)
parameter(X=1,Y=2)
parameter(EX=1,EY=2,EZ=3,J1X=4,J1Y=5,J1Z=6,J2X=7,J2Y=8,J2Z=9,J3X=10,J3Y=11,J3Z=12)
parameter(pi=3.14159265358979324)
parameter(e=1.602176565E-19)
parameter(g=1.6)
parameter(den=2.5E+3)
parameter(jph=4.5E-6)

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
call getarg(18,str(18))
call getarg(19,str(19))

read(str(1),*) totdim(1)
read(str(2),*) totdim(2)
read(str(3),*) totdim(3)
read(str(4),*) r
read(str(5),*) v0
read(str(6),*) lbhole
read(str(7),*) ubhole
read(str(8),*) xholec
read(str(9),*) yholec
read(str(10),*) rhole
read(str(11),*) theta
read(str(12),*) dnum
read(str(13),*) stride
read(str(14),*) N
read(str(15),*) conve
read(str(16),*) convj
read(str(17),*) sind
read(str(18),*) eind
read(str(19),*) dt


allocate(ej(12,-1:totdim(1),-1:totdim(2),-1:totdim(3)))
allocate(wrk(12,0:totdim(1),0:totdim(2),0:totdim(3)))
allocate(random(dnum,2))
allocate(pbuf(dnum))

!------------------- initialize field data
xl = 0; yl = 0; zl = 0;
xu = totdim(1); yu = totdim(2); zu = totdim(3);

!------------------- initialize dust data
ds = 4*pi*r**2
dv = (4/3)*pi*r**3
dm = dv*den
jjph = ds*0.5*jph
ns = 1; ne = dnum;
xholec2 = xholec + tan(theta*PI/180.0)*(ubhole-lbhole)
yholec2 = yholec

do m=ns,ne
   call system_clock(count=c)
   call random_seed(put=(/c/))
   call random_number(random)
   pbuf(m)%x = random(m,X)*totdim(1)
   pbuf(m)%y = random(m,Y)*totdim(2)
   pbuf(m)%vz = v0
   distance  = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
   distance2 = ((pbuf(m)%x-xholec2)**2+(pbuf(m)%y-yholec2)**2)**0.5
   ! upper hole
   if(rhole < distance) then
      pbuf(m)%z = ubhole
      pbuf(m)%q = e
   ! lower hole
   else
      pbuf(m)%z = lbhole
      ! shade
      if(rhole < distance2) then
         pbuf(m)%q = -e
      ! sunshine
      else
         pbuf(m)%q = e
      end if
   end if
end do

!----------- output initialized data
open(11,file='output.dat', status='replace')
      !--------- output VTK
      write(fvtk,'(a,i4.4,a)') 'dust00000000_id:', iind, '.vtk'
      open(20,file=fvtk,status='replace')

      !--------- VTK format
      write(20,"('# vtk DataFile Version 3.0')")
      write(20,'(a)') fvtk
      write(20,"('ASCII ')")
      write(20,"('DATASET UNSTRUCTURED_GRID')")
      write(20,"('POINTS ',i4,' float')") dnum
      do m=ns,ne
         write(20,"(3(f10.4,1x))") pbuf(m)%x, pbuf(m)%y, pbuf(m)%z
      end do
      write(20,"('CELLS ',i4,'  ',i5)") dnum, dnum*2
      do m=ns,ne
         write(20,"('1 ',i)") m-1
      end do
      write(20,"('CELL_TYPES ',i4)") dnum
      do m=ns,ne
         write(20,"('1')")
      end do
   
      write(20,"('POINT_DATA ',i4)") dnum

      !---------- SCALARS
      write(20,"('SCALARS charge float')")
      write(20,"('LOOKUP_TABLE default')")
      do m=ns,ne
         write(20,*) pbuf(m)%q 
      end do
      
      ! !---------- VECTORS
      ! write(20,"('VECTORS velocity float')")
      ! do m=ns,ne
      !    write(20,'(3f)') pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz
      ! end do
       
      close(20)
      iind = iind + 1

call hdfinit()

!------------------- import E-data
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'ex', sind, '-', eind ,'a.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'ey', sind, '-', eind ,'a.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'ez', sind, '-', eind ,'a.h5'

write(sdsnamex,'(a,i4.4)') fnorgx(1:2),sind
write(sdsnamey,'(a,i4.4)') fnorgy(1:2),sind
write(sdsnamez,'(a,i4.4)') fnorgz(1:2),sind

call hdfopen(fnorgx,ids1,ACC_R)
call hdfopen(fnorgy,ids2,ACC_R)
call hdfopen(fnorgz,ids3,ACC_R)

call read3d(ids1,sdsnamex,totdim+1,ej(EX,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids2,sdsnamey,totdim+1,ej(EY,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids3,sdsnamez,totdim+1,ej(EZ,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)

call hdfclose(ids1,ist)
call hdfclose(ids2,ist)
call hdfclose(ids3,ist)

!------------------- import J1-data
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'j1x', sind, '-', eind ,'a.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'j1y', sind, '-', eind ,'a.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'j1z', sind, '-', eind ,'a.h5'

write(sdsnamex,'(a,i4.4)') fnorgx(1:3),sind
write(sdsnamey,'(a,i4.4)') fnorgy(1:3),sind
write(sdsnamez,'(a,i4.4)') fnorgz(1:3),sind

call hdfopen(fnorgx,ids1,ACC_R)
call hdfopen(fnorgy,ids2,ACC_R)
call hdfopen(fnorgz,ids3,ACC_R)

call read3d(ids1,sdsnamex,totdim+1,ej(J1X,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids2,sdsnamey,totdim+1,ej(J1Y,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids3,sdsnamez,totdim+1,ej(J1Z,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)

call hdfclose(ids1,ist)
call hdfclose(ids2,ist)
call hdfclose(ids3,ist)

!------------------- import J2-data
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'j2x', sind, '-', eind ,'a.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'j2y', sind, '-', eind ,'a.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'j2z', sind, '-', eind ,'a.h5'

write(sdsnamex,'(a,i4.4)') fnorgx(1:3),sind
write(sdsnamey,'(a,i4.4)') fnorgy(1:3),sind
write(sdsnamez,'(a,i4.4)') fnorgz(1:3),sind

call hdfopen(fnorgx,ids1,ACC_R)
call hdfopen(fnorgy,ids2,ACC_R)
call hdfopen(fnorgz,ids3,ACC_R)

call read3d(ids1,sdsnamex,totdim+1,ej(J2X,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids2,sdsnamey,totdim+1,ej(J2Y,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids3,sdsnamez,totdim+1,ej(J2Z,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)

call hdfclose(ids1,ist)
call hdfclose(ids2,ist)
call hdfclose(ids3,ist)

!------------------- import J3-data
write(fnorgx,'(a,i3.3,a,i3.3,a)') 'j3x', sind, '-', eind ,'a.h5'
write(fnorgy,'(a,i3.3,a,i3.3,a)') 'j3y', sind, '-', eind ,'a.h5'
write(fnorgz,'(a,i3.3,a,i3.3,a)') 'j3z', sind, '-', eind ,'a.h5'

write(sdsnamex,'(a,i4.4)') fnorgx(1:3),sind
write(sdsnamey,'(a,i4.4)') fnorgy(1:3),sind
write(sdsnamez,'(a,i4.4)') fnorgz(1:3),sind

call hdfopen(fnorgx,ids1,ACC_R)
call hdfopen(fnorgy,ids2,ACC_R)
call hdfopen(fnorgz,ids3,ACC_R)

call read3d(ids1,sdsnamex,totdim+1,ej(J3X,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids2,sdsnamey,totdim+1,ej(J3Y,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
call read3d(ids3,sdsnamez,totdim+1,ej(J3Z,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)

call hdfclose(ids1,ist)
call hdfclose(ids2,ist)
call hdfclose(ids3,ist)

call hdffinalize()

!-------------------- S > R
ej(EX:EZ,:,:,:) = ej(EX:EZ,:,:,:)/conve
ej(J1X:J3Z,:,:,:) = ej(J1X:J3Z,:,:,:)/convj

!-------------------- relocation of pej-fields
do k=0,zu-zl
   do j=0,yu-yl
      do i=0,xu-xl
         wrk(EX,i,j,k) = 0.5d0*(ej(EX,i,j,k) + ej(EX,i-1,j,k))
         wrk(EY,i,j,k) = 0.5d0*(ej(EY,i,j,k) + ej(EY,i,j-1,k))
         wrk(EZ,i,j,k) = 0.5d0*(ej(EZ,i,j,k) + ej(EZ,i,j,k-1))
         wrk(J1X,i,j,k) = 0.5d0*(ej(J1X,i,j,k) + ej(J1X,i-1,j,k))
         wrk(J1Y,i,j,k) = 0.5d0*(ej(J1Y,i,j,k) + ej(J1Y,i,j-1,k))
         wrk(J1Z,i,j,k) = 0.5d0*(ej(J1Z,i,j,k) + ej(J1Z,i,j,k-1))
         wrk(J2X,i,j,k) = 0.5d0*(ej(J2X,i,j,k) + ej(J2X,i-1,j,k))
         wrk(J2Y,i,j,k) = 0.5d0*(ej(J2Y,i,j,k) + ej(J2Y,i,j-1,k))
         wrk(J2Z,i,j,k) = 0.5d0*(ej(J2Z,i,j,k) + ej(J2Z,i,j,k-1))
         wrk(J3X,i,j,k) = 0.5d0*(ej(J3X,i,j,k) + ej(J3X,i-1,j,k))
         wrk(J3Y,i,j,k) = 0.5d0*(ej(J3Y,i,j,k) + ej(J3Y,i,j-1,k))
         wrk(J3Z,i,j,k) = 0.5d0*(ej(J3Z,i,j,k) + ej(J3Z,i,j,k-1))
      end do
   end do
end do

!-------------------- update x,y,z,vx,vy,vz
do n = 1,N
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
      
       if(k==lbhole) then
          v5=0; v6=0; v7=0; v8=0;
       end if

      eez = wrk(EZ,i ,j ,k1)*v1 + wrk(EZ,i1,j ,k1)*v2 &
     &    + wrk(EZ,i1,j1,k1)*v3 + wrk(EZ,i ,j1,k1)*v4 &
     &    + wrk(EZ,i ,j ,k )*v5 + wrk(EZ,i1,j ,k )*v6 &
     &    + wrk(EZ,i1,j1,k )*v7 + wrk(EZ,i ,j1,k )*v8

      !----------- j1-fields interpolation
      jj1x = wrk(J1X,i ,j ,k1)*v1 + wrk(J1X,i1,j ,k1)*v2 &
     &    + wrk(J1X,i1,j1,k1)*v3 + wrk(J1X,i ,j1,k1)*v4 &
     &    + wrk(J1X,i ,j ,k )*v5 + wrk(J1X,i1,j ,k )*v6 &
     &    + wrk(J1X,i1,j1,k )*v7 + wrk(J1X,i ,j1,k )*v8
      jj1y = wrk(J1Y,i ,j ,k1)*v1 + wrk(J1Y,i1,j ,k1)*v2 &
     &    + wrk(J1Y,i1,j1,k1)*v3 + wrk(J1Y,i ,j1,k1)*v4 &
     &    + wrk(J1Y,i ,j ,k )*v5 + wrk(J1Y,i1,j ,k )*v6 &
     &    + wrk(J1Y,i1,j1,k )*v7 + wrk(J1Y,i ,j1,k )*v8
      jj1z = wrk(J1Z,i ,j ,k1)*v1 + wrk(J1Z,i1,j ,k1)*v2 &
     &    + wrk(J1Z,i1,j1,k1)*v3 + wrk(J1Z,i ,j1,k1)*v4 &
     &    + wrk(J1Z,i ,j ,k )*v5 + wrk(J1Z,i1,j ,k )*v6 &
     &    + wrk(J1Z,i1,j1,k )*v7 + wrk(J1Z,i ,j1,k )*v8

      !----------- j2-fields interpolation
      jj2x = wrk(J2X,i ,j ,k1)*v1 + wrk(J2X,i1,j ,k1)*v2 &
     &    + wrk(J2X,i1,j1,k1)*v3 + wrk(J2X,i ,j1,k1)*v4 &
     &    + wrk(J2X,i ,j ,k )*v5 + wrk(J2X,i1,j ,k )*v6 &
     &    + wrk(J2X,i1,j1,k )*v7 + wrk(J2X,i ,j1,k )*v8
      jj2y = wrk(J2Y,i ,j ,k1)*v1 + wrk(J2Y,i1,j ,k1)*v2 &
     &    + wrk(J2Y,i1,j1,k1)*v3 + wrk(J2Y,i ,j1,k1)*v4 &
     &    + wrk(J2Y,i ,j ,k )*v5 + wrk(J2Y,i1,j ,k )*v6 &
     &    + wrk(J2Y,i1,j1,k )*v7 + wrk(J2Y,i ,j1,k )*v8
      jj2z = wrk(J2Z,i ,j ,k1)*v1 + wrk(J2Z,i1,j ,k1)*v2 &
     &    + wrk(J2Z,i1,j1,k1)*v3 + wrk(J2Z,i ,j1,k1)*v4 &
     &    + wrk(J2Z,i ,j ,k )*v5 + wrk(J2Z,i1,j ,k )*v6 &
     &    + wrk(J2Z,i1,j1,k )*v7 + wrk(J2Z,i ,j1,k )*v8

      !----------- j3-fields interpolation
      jj3x = wrk(J3X,i ,j ,k1)*v1 + wrk(J3X,i1,j ,k1)*v2 &
     &    + wrk(J3X,i1,j1,k1)*v3 + wrk(J3X,i ,j1,k1)*v4 &
     &    + wrk(J3X,i ,j ,k )*v5 + wrk(J3X,i1,j ,k )*v6 &
     &    + wrk(J3X,i1,j1,k )*v7 + wrk(J3X,i ,j1,k )*v8
      jj3y = wrk(J3Y,i ,j ,k1)*v1 + wrk(J3Y,i1,j ,k1)*v2 &
     &    + wrk(J3Y,i1,j1,k1)*v3 + wrk(J3Y,i ,j1,k1)*v4 &
     &    + wrk(J3Y,i ,j ,k )*v5 + wrk(J3Y,i1,j ,k )*v6 &
     &    + wrk(J3Y,i1,j1,k )*v7 + wrk(J3Y,i ,j1,k )*v8
      jj3z = wrk(J3Z,i ,j ,k1)*v1 + wrk(J3Z,i1,j ,k1)*v2 &
     &    + wrk(J3Z,i1,j1,k1)*v3 + wrk(J3Z,i ,j1,k1)*v4 &
     &    + wrk(J3Z,i ,j ,k )*v5 + wrk(J3Z,i1,j ,k )*v6 &
     &    + wrk(J3Z,i1,j1,k )*v7 + wrk(J3Z,i ,j1,k )*v8

      !----------- caluculate electric current
      jj1x = jj1x * ds * 0.5
      jj1y = jj1y * ds * 0.5
      jj1z = jj1z * ds * 0.5
      jj2x = jj2x * ds * 0.5
      jj2y = jj2y * ds * 0.5
      jj2z = jj2z * ds * 0.5
      jj3x = jj3x * ds * 0.5
      jj3y = jj3y * ds * 0.5
      jj3z = jj3z * ds * 0.5

      !----------- caluculate J magnitude into dust
      jj1 = (jj1x**2+jj1y**2+jj1z**2)**0.5
      jj2 = (jj2x**2+jj2y**2+jj2z**2)**0.5
      jj3 = (jj3x**2+jj3y**2+jj3z**2)**0.5

      q = (jj2 - jj1 + jj3)*dt
      write(11,*) q

      call system_clock(count=c)
      call random_seed(put=(/c/))
      call random_number(a)
      if(q < 0) then
         a = a*(-e)
         if(q<a) pbuf(m)%q = pbuf(m)%q - e
      end if

      if(q <= 0) then
         a = a*e
         if(a<q) pbuf(m)%q = pbuf(m)%q + e
      end if

      !write(11,*) n, pbuf(1)%x, pbuf(1)%y, pbuf(1)%z, pbuf(1)%q

      !----------- update dust velocities
      pbuf(m)%vx = pbuf(m)%vx + (eex*pbuf(m)%q/dm)*dt
      pbuf(m)%vy = pbuf(m)%vy + (eey*pbuf(m)%q/dm)*dt
      pbuf(m)%vz = pbuf(m)%vz + (eez*pbuf(m)%q/dm)*dt - g*dt
 
      !----------- update dust positions
      pbuf(m)%x = pbuf(m)%x + pbuf(m)%vx*dt
      pbuf(m)%y = pbuf(m)%y + pbuf(m)%vy*dt
      pbuf(m)%z = pbuf(m)%z + pbuf(m)%vz*dt
      
      !----------- check dust positions
      if(pbuf(m)%x<xl.or.xu<pbuf(m)%x)     pbuf(m)%nid = -1
      if(pbuf(m)%y<yl.or.yu<pbuf(m)%y)     pbuf(m)%nid = -1
      if(pbuf(m)%z<lbhole.or.zu<pbuf(m)%z) pbuf(m)%nid = -1
      if(pbuf(m)%z<ubhole) then
         distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
         if(rhole < distance) then
            pbuf(m)%nid = -1
         end if
      end if

      !----------- reinitialize dust data
      if(pbuf(m)%nid == -1) then
         call system_clock(count=c)
         call random_seed(put=(/c/))
         call random_number(random)
         pbuf(m)%x = random(m,X)*totdim(1)
         pbuf(m)%y = random(m,Y)*totdim(2)
         pbuf(m)%vz = v0
         distance  = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
         distance2 = ((pbuf(m)%x-xholec2)**2+(pbuf(m)%y-yholec2)**2)**0.5
         ! upper hole
         if(rhole < distance) then
            pbuf(m)%z = ubhole
            pbuf(m)%q = e
            ! lower hole
         else
            pbuf(m)%z = lbhole
            ! shade
            if(rhole < distance2) then
               pbuf(m)%q = -e
               ! sunshine
            else
               pbuf(m)%q = e
            end if
         end if
         pbuf(m)%nid = 1
      end if

   end do

   !----------- output dust data
   if(mod(n,stride)==0) then
     
      write(fvtk,'(a,i4.4,a)') 'dust00000000_id:', iind, '.vtk'
      open(20,file=fvtk,status='replace')

      !--------- VTK format
      write(20,"('# vtk DataFile Version 3.0')")
      write(20,'(a)') fvtk
      write(20,"('ASCII ')")
      write(20,"('DATASET UNSTRUCTURED_GRID')")
      write(20,"('POINTS ',i4,' float')") dnum
      do m=ns,ne
         write(20,"(3(f10.4,1x))") pbuf(m)%x, pbuf(m)%y, pbuf(m)%z
      end do
      write(20,"('CELLS ',i4,'  ',i5)") dnum, dnum*2
      do m=ns,ne
         write(20,"('1 ',i)") m-1
      end do
      write(20,"('CELL_TYPES ',i4)") dnum
      do m=ns,ne
         write(20,"('1')")
      end do
   
      write(20,"('POINT_DATA ',i4)") dnum

      !---------- SCALARS
      write(20,"('SCALARS charge float')")
      write(20,"('LOOKUP_TABLE default')")
      do m=ns,ne
         write(20,*) pbuf(m)%q 
      end do 
      
      ! !---------- VECTORS
      ! write(20,"('VECTORS velocity float')")
      ! do m=ns,ne
      !    write(20,*) pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz
      ! end do
     
        
      close(20)
      iind = iind + 1
   end if
   
   write(6,*) n
end do
close(11)
end program dmotion2
