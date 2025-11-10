program dmotion_simple

  implicit none

  type dust
     real(kind=8) :: x,y,z,vx,vy,vz
     integer(kind=4) :: q,nid
  end type dust

  integer(kind=4) :: i,j,k,i1,j1,k1,m,n,ns,dnum,cc,iind
  integer(kind=4) :: stvtk,stdat,lavtk,ladat,nx,ny,nz,oml,dat,vtk,reflect,fix,periodic,xyfix
  integer(kind=8) :: totdim(3),xl,yl,zl,xu,yu,zu
  integer(kind=8) :: lbhole,ubhole,xholec,yholec,rhole,xholec2,yholec2
  character(len=40) :: fvtk
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2,xlocal,ylocal,zlocal,time0,time1
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8,dt
  real(kind=8) :: eex,eey,eez,Te,Ti,Tph,Vd,Cd,e0,bk
  real(kind=8) :: q,r,v0,g,pi,dm,dv,ds,den,theta,distance,e
  real(kind=8) :: xb,yb,zb
  type(dust),allocatable :: pbuf(:)

  parameter(pi=3.14159265358979324)
  parameter(e=1.602176565E-19)
  parameter(e0=8.854187817E-12 )
  parameter(bk=1.38064852E-23)

  ! Simple parameter settings (no input file)
  nx = 160; ny = 160; nz = 800
  ns = 1
  iind = 0
  totdim(1) = nx; totdim(2) = ny; totdim(3) = nz

  !------------------- initialize field data
  xl = 0; yl = 0; zl = 0
  xu = totdim(1); yu = totdim(2); zu = totdim(3)

  !------------------- dust parameters
  r = 2.0E-8         ! dust radius 
  den = 2.5E+3       ! dust density
  v0 = 0             ! initial velocity Z direction
  dnum = 100         ! number of dusts
  dt = 1.0E-4        ! time step
  g = 1.6            ! gravity

  ! Geometry parameters
  xholec = 80; yholec = 80; rhole = 25
  lbhole = 15; ubhole = 60; theta = 30

  !------------------- initialize dust data
  ds = pi*r**2
  dv = (4/3)*pi*r**3
  dm = dv*den
  xholec2 = xholec + tan(theta*pi/180.0)*(ubhole-lbhole)
  yholec2 = yholec

  allocate(pbuf(dnum))

  ! Initialize dust positions
  do m=ns,dnum
     call random_seed()
     call random_number(xlocal)
     call random_number(ylocal)
     pbuf(m)%x = xlocal*totdim(1)
     pbuf(m)%y = ylocal*totdim(2)
     pbuf(m)%vz = v0
     distance = ((pbuf(m)%x-xholec)**2+(pbuf(m)%y-yholec)**2)**0.5
     ! upper hole
     if(rhole < distance) then
        pbuf(m)%z = ubhole
        pbuf(m)%q = 1
     else
        pbuf(m)%z = lbhole
        pbuf(m)%q = 1
     end if
     pbuf(m)%nid = 1
  end do

  ! Simple simulation loop (no electric field, just gravity)
  do n = 1, 10000
     
     do m=ns,dnum
        if(pbuf(m)%nid.eq.-1) cycle

        ! Simple gravity only
        pbuf(m)%vz = pbuf(m)%vz - g*dt
        
        ! Update positions
        pbuf(m)%x = pbuf(m)%x + pbuf(m)%vx*dt
        pbuf(m)%y = pbuf(m)%y + pbuf(m)%vy*dt
        pbuf(m)%z = pbuf(m)%z + pbuf(m)%vz*dt

        ! Check boundaries
        if(pbuf(m)%x<xl.or.xu<pbuf(m)%x.or.pbuf(m)%y<yl.or.yu<pbuf(m)%y) then
           pbuf(m)%nid = -1
        end if
        if(pbuf(m)%z<lbhole.or.zu<pbuf(m)%z) then
           pbuf(m)%nid = -1
        end if
        
     end do

     ! Remove dead particles
     do m=ns,dnum
        if(pbuf(m)%nid == -1) then
           if(m == dnum) then
              dnum = dnum - 1
           else
              pbuf(m) = pbuf(dnum)
              dnum = dnum - 1
           end if
        end if
     end do

     ! Output VTK file every 1000 steps
     if(mod(n,1000)==0) then
        write(fvtk,'(a,i4.4,a)') 'dust_simple_', n/1000, '.vtk'
        open(20,file=fvtk,status='replace')

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
        write(20,"('SCALARS charge float')")
        write(20,"('LOOKUP_TABLE default')")
        do m=ns,dnum
           write(20,*) pbuf(m)%q
        end do

        close(20)
        write(6,*) 'Output: ', fvtk, ' particles: ', dnum
     end if

     if(dnum==0) exit
     
  end do

  write(6,*) 'Simulation finished'
end program dmotion_simple
