program dmotion14

  implicit none

  type dust
     real(kind=8) :: x,y,z,vx,vy,vz
     integer(kind=4) :: q,nid
  end type dust

  integer(kind=4) :: i,j,k,i1,j1,k1,m,n,sind,eind,ns,dnum,cc,edind,nstep
  integer(kind=4) :: stvtk,stdat,stnden,lavtk,ladat,nx,ny,nz,jmode,oml,dat,vtk,nden,reflect,fix,periodic,xyfix
  integer(kind=4) :: ACC_R, ACC_W, ACC_RW, ACC_C
  integer(kind=4) :: X,Y,EX,EY,EZ,J1X,J1Y,J1Z,J2X,J2Y,J2Z,J3X,J3Y,J3Z,VDE,VDI,VDP
  integer(kind=4) :: ids1, ids2, ids3, ist, ist1, ist2,iind,iind2,idt
  integer(kind=8) :: totdim(3),xl,yl,zl,xu,yu,zu
  integer(kind=8) :: lbhole,ubhole,xholec,yholec,rhole,xholec2,yholec2
  character(len=40) :: fnorgx, fnorgy, fnorgz, fvtk
  character(len=16) :: sdsnamex, sdsnamey, sdsnamez, sdsname
  real(kind=8),allocatable :: ej(:,:,:,:),wrk(:,:,:,:),random(:,:),random2(:),domain(:,:,:)
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2,xlocal,ylocal,zlocal,time0,time1,omp_get_wtime
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8,dt
  real(kind=8) :: eex,eey,eez,jj1,jj1x,jj1y,jj1z,jj2,jj2x,jj2y,jj2z,jj3,jj3x,jj3y,jj3z,jph,jjph,Te,Ti,Tph,Vd,Cd,e0,bk
  real(kind=8) :: q,r,v0,g,conve,convj,pi,dm,dv,ds,den,theta,theta2,theta3,distance,distance2,distanceb,e
  real(kind=8) :: intx,inty,intx1,inty1,intx2,inty2,xc,yc,xd,yd,dss,apl,bpl,t,s,kk,d,cpx,cpy,cqx,cqy,a,b,c,xb,yb,zb,ll
  character(len=40) :: fnd
  type(dust),allocatable :: pbuf(:)

  parameter(ACC_R=1, ACC_W=2, ACC_RW=3, ACC_C=4)
  parameter(X=1,Y=2)
  parameter(EX=1,EY=2,EZ=3,J1X=4,J1Y=5,J1Z=6,J2X=7,J2Y=8,J2Z=9,J3X=10,J3Y=11,J3Z=12,VDE=13,VDI=14,VDP=15)
  parameter(pi=3.14159265358979324)
  parameter(e=1.602176565E-19)
  parameter(e0=8.854187817E-12 )
  parameter(bk=1.38064852E-23)

  namelist /jobcon/ nstep,dat,vtk,nden,oml,jmode,reflect,fix,periodic,xyfix
  namelist /inphdf/ sind,eind
  namelist /digcon/ stvtk,lavtk,stdat,ladat,stnden
  namelist /conv/   conve,convj
  namelist /tmgrid/ dt,nx,ny,nz
  namelist /plasma/ Te,Ti,Tph,jph
  namelist /grain/  r,den,v0,dnum
  namelist /ptcond/ xholec,yholec,rhole,lbhole,ubhole,theta,g

  open(10,file='dust.inp')
  read(10,jobcon)
  read(10,inphdf)
  read(10,digcon)
  read(10,conv)
  read(10,tmgrid)
  read(10,plasma)
  read(10,grain)
  read(10,ptcond)
  close(10)

  ns = 1; n = 1
  iind  = 0; iind2 = 0;
  totdim(1) = nx; totdim(2) = ny; totdim(3) = nz;

  !------------------- initialize field data
  xl = 0; yl = 0; zl = 0;
  xu = totdim(1); yu = totdim(2); zu = totdim(3);

  !------------------- initialize dust data
  ds = pi*r**2
  dv = (4/3)*pi*r**3
  dm = dv*den
  jjph = jph*ds

  xholec2 = xholec + tan(theta*pi/180.0)*(ubhole-lbhole)
  yholec2 = yholec

  Te  = Te  * 11604 ! [eV] to [K]
  Ti  = Ti  * 11604 ! [eV] to [K]
  Tph = Tph * 11604 ! [eV] to [K]

  allocate(ej(15,-1:totdim(1),-1:totdim(2),-1:totdim(3)))
  allocate(wrk(15,0:totdim(1),0:totdim(2),0:totdim(3)))
  allocate(random(dnum,2))
  allocate(random2(dnum))
  allocate(pbuf(dnum))
  allocate(domain(0:totdim(1),0:totdim(2),0:totdim(3)))

  call hdfinit()

  time0=omp_get_wtime()
  do m=ns,dnum
     call system_clock(count=cc)
     call random_seed(put=(/cc/))
     call random_number(random)
     pbuf(m)%x = random(m,X)*totdim(1)
     pbuf(m)%y = random(m,Y)*totdim(2)
     pbuf(m)%vz = v0
     distance  = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
     distance2 = ((pbuf(m)%x-xholec2)**2+(pbuf(m)%y-yholec2)**2)**0.5
     ! upper hole
     if(rhole < distance) then
        pbuf(m)%z = ubhole
        pbuf(m)%q = 1
        ! lower hole
     else
        pbuf(m)%z = lbhole
        ! shade
        if(rhole < distance2) then
           pbuf(m)%q = -1
           ! sunshine
        else
           pbuf(m)%q = 1
        end if
     end if
  end do


  !----------- output initialized data
  if(dat==1)then
     do m=ns,dnum
        write(fvtk,'(a,i4.4,a)') 'dust:', m, '.dat'
        open(11,file=fvtk, status='replace')
        if(ladat==1) then
           write(11,*) 0,pbuf(m)%z-ubhole,pbuf(m)%vz,pbuf(m)%q
        end if
        close(11)
     end do
  end if

  if(vtk==1.and.lavtk==1) then
     !--------- output VTK
     write(fvtk,'(a,i4.4,a)') 'dust00000000_id:', iind, '.vtk'
     open(20,file=fvtk,status='replace')

     !--------- VTK format
     write(20,"('# vtk DataFile Version 3.0')")
     write(20,'(a)') fvtk
     write(20,"('ASCII ')")
     write(20,"('DATASET UNSTRUCTURED_GRID')")
     write(20,"('POINTS ',i7,' float')") dnum
     do m=ns,dnum
        write(20,"(3(f10.4,1x))") pbuf(m)%x, pbuf(m)%y, pbuf(m)%z
     end do
     write(20,"('CELLS ',i7,'  ',i7)") dnum, dnum*2
     do m=ns,dnum
        write(20,"('1 ',i)") m-1
     end do
     write(20,"('CELL_TYPES ',i7)") dnum
     do m=ns,dnum
        write(20,"('1')")
     end do

     write(20,"('POINT_DATA ',i7)") dnum

     !---------- SCALARS
     write(20,"('SCALARS charge float')")
     write(20,"('LOOKUP_TABLE default')")
     do m=ns,dnum
        write(20,*) pbuf(m)%q
     end do

     ! !---------- VECTORS
     ! write(20,"('VECTORS velocity float')")
     ! do m=ns,dnum
     !    write(20,'(3f)') pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz
     ! end do

     close(20)

     iind = iind + 1
  end if

  if(nden==1) then
     edind = nstep/stnden
     write(fnd,'(a,i4.4,a,i4.4,a)') 'fnd_id:', 0,'-', edind ,'.h5' ! (input-jx)
     call hdfopen(fnd,idt,ACC_C) ! Open HDF5 file for merged data (output)

     do m=ns,dnum
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

        domain(i ,j ,k1) = domain(i ,j ,k1) + v1
        domain(i1,j ,k1) = domain(i1,j ,k1) + v2
        domain(i1,j1,k1) = domain(i1,j1,k1) + v3
        domain(i ,j1,k1) = domain(i ,j1,k1) + v4
        domain(i ,j ,k ) = domain(i ,j ,k ) + v5
        domain(i1,j ,k ) = domain(i1,j ,k ) + v6
        domain(i1,j1,k ) = domain(i1,j1,k ) + v7
        domain(i ,j1,k ) = domain(i ,j1,k ) + v8      
     end do

     write(sdsname,'(i4.4)') iind2 ! Define dataset name
     call wrt3d(idt,sdsname,totdim+1,domain,ist1,ist2) ! Write whole domain data
     iind2 = iind2 + 1
  end if

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

  if(jmode.eq.1) then     
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
  end if

  if(jmode.eq.2) then
     write(fnorgx,'(a,i3.3,a,i3.3,a)') 'vde', sind, '-', eind ,'a.h5'
     write(fnorgy,'(a,i3.3,a,i3.3,a)') 'vdi', sind, '-', eind ,'a.h5'
     write(fnorgz,'(a,i3.3,a,i3.3,a)') 'vdp', sind, '-', eind ,'a.h5'

     write(sdsnamex,'(a,i4.4)') fnorgx(1:3),sind
     write(sdsnamey,'(a,i4.4)') fnorgy(1:3),sind
     write(sdsnamez,'(a,i4.4)') fnorgz(1:3),sind

     call hdfopen(fnorgx,ids1,ACC_R)
     call hdfopen(fnorgy,ids2,ACC_R)
     call hdfopen(fnorgz,ids3,ACC_R)

     call read3d(ids1,sdsnamex,totdim+1,ej(VDE,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
     call read3d(ids2,sdsnamey,totdim+1,ej(VDI,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
     call read3d(ids3,sdsnamez,totdim+1,ej(VDP,0:totdim(1),0:totdim(2),0:totdim(3)),ist1,ist2)
     ej(VDE:VDP,:,:,:) = abs(ej(VDE:VDP,:,:,:))

     call hdfclose(ids1,ist)
     call hdfclose(ids2,ist)
     call hdfclose(ids3,ist)
  end if

  !-------------------- S > R
  ej(EX:EZ,:,:,:) = ej(EX:EZ,:,:,:)/conve
  ej(J1X:VDP,:,:,:) = ej(J1X:VDP,:,:,:)/convj

  !-------------------- relocation of pej-fields
  do k=0,zu-zl
     do j=0,yu-yl
        do i=0,xu-xl
           wrk(EX,i,j,k) = 0.5d0*(ej(EX,i,j,k) + ej(EX,i-1,j,k))
           wrk(EY,i,j,k) = 0.5d0*(ej(EY,i,j,k) + ej(EY,i,j-1,k))
           wrk(EZ,i,j,k) = 0.5d0*(ej(EZ,i,j,k) + ej(EZ,i,j,k-1))
           if(jmode.eq.1) then
              wrk(J1X,i,j,k) = 0.5d0*(ej(J1X,i,j,k) + ej(J1X,i-1,j,k))
              wrk(J1Y,i,j,k) = 0.5d0*(ej(J1Y,i,j,k) + ej(J1Y,i,j-1,k))
              wrk(J1Z,i,j,k) = 0.5d0*(ej(J1Z,i,j,k) + ej(J1Z,i,j,k-1))
              wrk(J2X,i,j,k) = 0.5d0*(ej(J2X,i,j,k) + ej(J2X,i-1,j,k))
              wrk(J2Y,i,j,k) = 0.5d0*(ej(J2Y,i,j,k) + ej(J2Y,i,j-1,k))
              wrk(J2Z,i,j,k) = 0.5d0*(ej(J2Z,i,j,k) + ej(J2Z,i,j,k-1))
              wrk(J3X,i,j,k) = 0.5d0*(ej(J3X,i,j,k) + ej(J3X,i-1,j,k))
              wrk(J3Y,i,j,k) = 0.5d0*(ej(J3Y,i,j,k) + ej(J3Y,i,j-1,k))
              wrk(J3Z,i,j,k) = 0.5d0*(ej(J3Z,i,j,k) + ej(J3Z,i,j,k-1))
           end if
           if(jmode.eq.2) then
              wrk(VDE,i,j,k) = 0.5d0*(ej(VDE,i,j,k) + ej(VDE,i-1,j,k))
              wrk(VDI,i,j,k) = 0.5d0*(ej(VDI,i,j,k) + ej(VDI,i,j-1,k))
              wrk(VDP,i,j,k) = 0.5d0*(ej(VDP,i,j,k) + ej(VDP,i,j,k-1))
           end if
        end do
     end do
  end do

  do k=0,zu-zl
     do j=0,yu-yl
        do i=0,xu-xl
           distance  = ((i-xholec)**2+(j-yholec)**2)**0.5
           if(rhole<distance.and.k==ubhole) then
              wrk(EZ,i,j,k) = 2*wrk(EZ,i,j,k+1) - wrk(EZ,i,j,k+2)
           else if(distance<rhole.and.k==lbhole) then
              wrk(EZ,i,j,k) = 2*wrk(EZ,i,j,k+1) - wrk(EZ,i,j,k+2)
           end if
        end do
     end do
  end do


  !-------------------- update x,y,z,vx,vy,vz
  do n = 1,nstep

     call system_clock(count=cc)
     call random_seed(put=(/cc/))
     call random_number(random2)
     
     !$omp parallel do default(firstprivate) shared(ns,dnum,pbuf,xholec,yholec,theta,lbhole,ubhole,wrk,ds,r,rhole,Te,Ti,Tph,jjph,dt,oml,xyfix,g,dm,reflect,periodic,fix,xl,xu,yl,yu,zl,zu)
     do m=ns,dnum
        if(pbuf(m)%nid.eq.-1) cycle

        xholec2 = xholec + tan(theta*pi/180.0)*(ubhole-pbuf(m)%z)
        yholec2 = yholec
        distance  = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
        distance2 = ((pbuf(m)%x-xholec2)**2+(pbuf(m)%y-yholec2)**2)**0.5

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

        !if(k==lbhole) then
        !   v5=0; v6=0; v7=0; v8=0;
        !end if

        eez = wrk(EZ,i ,j ,k1)*v1 + wrk(EZ,i1,j ,k1)*v2 &
             &    + wrk(EZ,i1,j1,k1)*v3 + wrk(EZ,i ,j1,k1)*v4 &
             &    + wrk(EZ,i ,j ,k )*v5 + wrk(EZ,i1,j ,k )*v6 &
             &    + wrk(EZ,i1,j1,k )*v7 + wrk(EZ,i ,j1,k )*v8

        if(jmode.eq.1) then
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

           !----------- caluculate J magnitude into dust
           jj1 = (jj1x**2+jj1y**2+jj1z**2)**0.5*ds
           jj2 = (jj2x**2+jj2y**2+jj2z**2)**0.5*ds
           jj3 = (jj3x**2+jj3y**2+jj3z**2)**0.5*ds
        end if

        if(jmode.eq.2) then
           jj1 = wrk(VDE,i ,j ,k1)*v1 + wrk(VDE,i1,j ,k1)*v2 &
                &    + wrk(VDE,i1,j1,k1)*v3 + wrk(VDE,i ,j1,k1)*v4 &
                &    + wrk(VDE,i ,j ,k )*v5 + wrk(VDE,i1,j ,k )*v6 &
                &    + wrk(VDE,i1,j1,k )*v7 + wrk(VDE,i ,j1,k )*v8
           jj2 = wrk(VDI,i ,j ,k1)*v1 + wrk(VDI,i1,j ,k1)*v2 &
                &    + wrk(VDI,i1,j1,k1)*v3 + wrk(VDI,i ,j1,k1)*v4 &
                &    + wrk(VDI,i ,j ,k )*v5 + wrk(VDI,i1,j ,k )*v6 &
                &    + wrk(VDI,i1,j1,k )*v7 + wrk(VDI,i ,j1,k )*v8
           jj3 = wrk(VDP,i ,j ,k1)*v1 + wrk(VDP,i1,j ,k1)*v2 &
                &    + wrk(VDP,i1,j1,k1)*v3 + wrk(VDP,i ,j1,k1)*v4 &
                &    + wrk(VDP,i ,j ,k )*v5 + wrk(VDP,i1,j ,k )*v6 &
                &    + wrk(VDP,i1,j1,k )*v7 + wrk(VDP,i ,j1,k )*v8

           jj1 = jj1*4*ds
           jj2 = jj2*4*ds
           jj3 = jj3*4*ds
        end if

        !----------- caluculate dust potential
        Cd = 4*pi*e0*r
        Vd = (pbuf(m)%q*e)/Cd
        !write(11,*) 'Vd',Vd

        if(oml.eq.1) then
           !----------- caluculate charge
           !upper hole
           if(rhole < distance) then
              if(0<Vd) then
                 q = (jjph*exp(-(e*Vd)/(bk*Tph))*(1+(e*Vd)/(bk*Tph)) + jj2*exp(-(e*Vd)/(bk*Ti)) - jj1*(1+(e*Vd)/(bk*Te)) - jj3*(1+(e*Vd)/(bk*Tph)))*dt
                 !q = (jjph + jj2- jj1 - jj3)*dt
              else
                 q = (jjph + jj2*(1-(e*Vd)/(bk*Te)) - jj1*exp((e*Vd)/(bk*Ti)) - jj3*exp((e*Vd)/(bk*Tph)))*dt
              end if
              ! lower hole  
           else
              ! shade
              if(rhole < distance2) then
                 if(0<Vd) then
                    q = (jj2*exp(-(e*Vd)/(bk*Ti)) - jj1*(1+(e*Vd)/(bk*Te)) - jj3*(1+(e*Vd)/(bk*Tph)))*dt
                 else
                    q = (jj2*(1-(e*Vd)/(bk*Te)) - jj1*exp((e*Vd)/(bk*Ti)) - jj3*exp((e*Vd)/(bk*Tph)))*dt
                 end if
                 ! sunshine
              else 
                 if(0<Vd) then
                    q = (jjph*exp(-(e*Vd)/(bk*Tph))*(1+(e*Vd)/(bk*Tph)) + jj2*exp(-(e*Vd)/(bk*Ti)) - jj1*(1+(e*Vd)/(bk*Te)) - jj3*(1+(e*Vd)/(bk*Tph)))*dt
                 else
                    q = (jjph + jj2*(1-(e*Vd)/(bk*Te)) - jj1*exp((e*Vd)/(bk*Ti)) - jj3*exp((e*Vd)/(bk*Tph)))*dt
                 end if
              end if
           end if
        end if

        if(oml.eq.2) then
           !upper hole
           if(rhole < distance) then
              q = (jjph + jj2 - jj1 - jj3)*dt
              ! lower hole
           else
              ! shade
              if(rhole < distance2) then
                 q = (jj2 - jj1 - jj3)*dt
                 ! sunshine
              else
                 q = (jjph + jj2 - jj1 - jj3)*dt
              end if
           end if
        end if

        if(q < 0) then
           random2(m) = random2(m)*(-e)
           if(q<random2(m)) pbuf(m)%q = pbuf(m)%q - 1
        else  
           random2(m) = random2(m)*e
           !write(11,*) random2(m),q
           if(random2(m)<q) pbuf(m)%q = pbuf(m)%q + 1
        end if

        !----------- update dust velocities
        if(xyfix.eq.1) then
           pbuf(m)%vz = pbuf(m)%vz + (eez*pbuf(m)%q*e/dm)*dt - g*dt
        else
           pbuf(m)%vx = pbuf(m)%vx + (eex*pbuf(m)%q*e/dm)*dt
           pbuf(m)%vy = pbuf(m)%vy + (eey*pbuf(m)%q*e/dm)*dt
           pbuf(m)%vz = pbuf(m)%vz + (eez*pbuf(m)%q*e/dm)*dt - g*dt
        end if

        xb = pbuf(m)%x
        yb = pbuf(m)%y
        zb = pbuf(m)%z

        !----------- update dust positions
        pbuf(m)%x = pbuf(m)%x + pbuf(m)%vx*dt
        pbuf(m)%y = pbuf(m)%y + pbuf(m)%vy*dt
        pbuf(m)%z = pbuf(m)%z + pbuf(m)%vz*dt
        distanceb = ((xb-xholec)**2+(yb-yholec)**2)**0.5
        distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5

        if(reflect.eq.1) then
           !----------- reflect dusts
           if(lbhole<pbuf(m)%z.and.pbuf(m)%z<ubhole.and.rhole<distance.and.distanceb<rhole) then
              
              !----------- calculate a point of interaction
              a = pbuf(m)%y - yb
              b = xb - pbuf(m)%x
              c = pbuf(m)%x*yb - xb*pbuf(m)%y
              
              ll = a*a+b*b
              kk = a*xholec+b*yholec+c
              d = ll*rhole*rhole-kk*kk
              
              dss = sqrt(d)
              apl = a/ll
              bpl = b/ll
              xc = xholec-apl*kk
              yc = yholec-bpl*kk
              xd = bpl*dss
              yd = apl*dss
              
              intx1 = xc-xd
              inty1 = yc+yd
              intx2 = xc+xd
              inty2 = yc-yd

              if((pbuf(m)%x-intx1)**2+(pbuf(m)%y-inty1)**2 < (pbuf(m)%x-intx2)**2+(pbuf(m)%y-inty2)**2) then
                 intx = intx1
                 inty = inty1
              else
                 intx = intx2
                 inty = inty2
              end if

              !----------- calculate an angle
              cpx = xholec - intx
              cpy = yholec - inty
              cqx = xb - intx
              cqy = yb - inty
              s = cpx*cqy - cpy*cqx
              t = cpx*cqx + cpy*cqy
              theta2 = atan2(s,t)

              if(s>0.0) then
                 theta3 = (pi/2 - theta2)*2
              else
                 theta3 = -(pi/2 - theta2)*2
              end if

              !----------- update dust point
              pbuf(m)%x = cos(theta3)*pbuf(m)%x-sin(theta3)*pbuf(m)%y+intx-intx*cos(theta3)+inty*sin(theta3)
              pbuf(m)%y = sin(theta3)*pbuf(m)%x+cos(theta3)*pbuf(m)%y+inty-intx*sin(theta3)-inty*cos(theta3)
              pbuf(m)%vx = cos(theta3)*pbuf(m)%vx-sin(theta3)*pbuf(m)%vy
              pbuf(m)%vy = sin(theta3)*pbuf(m)%vx+cos(theta3)*pbuf(m)%vy
              distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
           end if
        end if

        if(periodic.eq.1) then
           !----------- periodic boundary
           if(ubhole<=pbuf(m)%z) then
              if(pbuf(m)%x<xl) pbuf(m)%x = pbuf(m)%x + xu
              if(xu<pbuf(m)%x) pbuf(m)%x = pbuf(m)%x - xu
              if(pbuf(m)%y<yl) pbuf(m)%y = pbuf(m)%y + yu
              if(yu<pbuf(m)%y) pbuf(m)%y = pbuf(m)%y - yu
           end if
        end if

        if(fix.eq.1) then
           !----------- fix dusts
           if(pbuf(m)%z<lbhole.and.distance<rhole) then
              pbuf(m)%x = (lbhole*(xb-pbuf(m)%x)-(xb*pbuf(m)%z-pbuf(m)%x*zb))/(zb-pbuf(m)%z)
              pbuf(m)%y = (lbhole*(yb-pbuf(m)%y)-(yb*pbuf(m)%z-pbuf(m)%y*zb))/(zb-pbuf(m)%z)
              pbuf(m)%z = lbhole
              distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
           else if(pbuf(m)%z<ubhole.and.rhole<distance) then
              pbuf(m)%x = (ubhole*(xb-pbuf(m)%x)-(xb*pbuf(m)%z-pbuf(m)%x*zb))/(zb-pbuf(m)%z)
              pbuf(m)%y = (ubhole*(yb-pbuf(m)%y)-(yb*pbuf(m)%z-pbuf(m)%y*zb))/(zb-pbuf(m)%z)
              pbuf(m)%z = ubhole
              distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
           end if
        end if

        !----------- check dust positions
        if(pbuf(m)%x<xl.or.xu<pbuf(m)%x.or.pbuf(m)%y<yl.or.yu<pbuf(m)%y.or.zu<pbuf(m)%z.or.pbuf(m)%z<lbhole.or.(pbuf(m)%z<ubhole.and.rhole<distance)) pbuf(m)%nid = -1
        
     end do
     !$omp end parallel do
     
     if(dat.ne.1) then
        do m=ns,dnum
           if(pbuf(m)%nid == -1) then
              if(m == dnum) then
                 dnum = dnum - 1
              else
                 pbuf(m)%x=pbuf(dnum)%x
                 pbuf(m)%y=pbuf(dnum)%y
                 pbuf(m)%z=pbuf(dnum)%z
                 pbuf(m)%vx=pbuf(dnum)%vx
                 pbuf(m)%vy=pbuf(dnum)%vy
                 pbuf(m)%vz=pbuf(dnum)%vz
                 pbuf(m)%q=pbuf(dnum)%q
                 pbuf(m)%nid = pbuf(dnum)%nid
                 dnum = dnum - 1
              end if
           end if
        end do
     end if
    
     !----------- output dust data
     if(dnum==0) ns = 0

     if(dat==1) then
        if(mod(n,stdat)==0.and.n>=ladat) then
           do m=ns,dnum
              if(pbuf(m)%nid.eq.-1) cycle
              write(fvtk,'(a,i4.4,a)') 'dust:', m, '.dat'
              open(11,file=fvtk, position='append')
              write(11,*) (n-ladat)*dt,pbuf(m)%z-ubhole,pbuf(m)%vz,pbuf(m)%q
              close(11)
           end do
        end if
     end if

     if(vtk.eq.1) then
        if((mod(n,stvtk)==0.and.n>=lavtk).or.dnum==0) then
           write(fvtk,'(a,i4.4,a)') 'dust00000000_id:', iind, '.vtk'
           open(20,file=fvtk,status='replace')

           !--------- VTK format
           write(20,"('# vtk DataFile Version 3.0')")
           write(20,'(a)') fvtk
           write(20,"('ASCII ')")
           write(20,"('DATASET UNSTRUCTURED_GRID')")
           write(20,"('POINTS ',i7,' float')") dnum
           do m=ns,dnum
              write(20,"(3(f10.4,1x))") pbuf(m)%x, pbuf(m)%y, pbuf(m)%z
           end do
           write(20,"('CELLS ',i7,'  ',i7)") dnum, dnum*2
           do m=ns,dnum
              write(20,"('1 ',i)") m-1
           end do
           write(20,"('CELL_TYPES ',i7)") dnum
           do m=ns,dnum
              write(20,"('1')")
           end do

           write(20,"('POINT_DATA ',i7)") dnum

           !---------- SCALARS
           write(20,"('SCALARS charge float')")
           write(20,"('LOOKUP_TABLE default')")
           do m=ns,dnum
              write(20,*) pbuf(m)%q
           end do

           ! !---------- VECTORS
           ! write(20,"('VECTORS velocity float')")
           ! do m=ns,dnum
           !    write(20,*) pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz
           ! end do

           close(20)

           iind = iind + 1
        end if
     end if

     if(nden.eq.1) then
        if(mod(n,stnden)==0.or.dnum==0) then
           domain(:,:,:)=0

           do m=ns,dnum
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

              domain(i ,j ,k1) = domain(i ,j ,k1) + v1
              domain(i1,j ,k1) = domain(i1,j ,k1) + v2
              domain(i1,j1,k1) = domain(i1,j1,k1) + v3
              domain(i ,j1,k1) = domain(i ,j1,k1) + v4
              domain(i ,j ,k ) = domain(i ,j ,k ) + v5
              domain(i1,j ,k ) = domain(i1,j ,k ) + v6
              domain(i1,j1,k ) = domain(i1,j1,k ) + v7
              domain(i ,j1,k ) = domain(i ,j1,k ) + v8      
           end do

           write(sdsname,'(i4.4)') iind2 ! Define dataset name
           call wrt3d(idt,sdsname,totdim+1,domain,ist1,ist2) ! Write whole domain data
           iind2 = iind2 + 1
        end if
     end if

     write(6,*) n

     if(dnum==0) exit

  end do
  
  call hdfclose(idt,ist)
  call hdffinalize()

  time1 = omp_get_wtime()

  write(6,*) time1 - time0

  write(6,*) 'Finish'
end program dmotion14
