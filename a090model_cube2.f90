program model_column
  implicit none
  ! geometric parameters
  real :: h0,r0,h1,r1,l,a !20210406UENO
  ! SPH parameters
  integer :: np
  real :: dp,hsml
  real,allocatable :: pos(:,:)
  integer,allocatable :: fix(:,:),mat(:)
  ! material parameters
  real :: yng,poisn,rho,c,phi,visc
  ! simulation parameters
  real :: dt
  integer :: nstep
  ! program parameter
  integer :: i,j,k,m
  integer :: maxnp
  real :: x,y,z,z1,r

  a=0.90
  r0= 0.05   !9.7*0.01
  r1=0.06  !20210406UENO
  h0=0.1
  h1=0.15 !20210406UENO
  l=0.1

  dp=0.005
  hsml=dp*1.35

  rho=2.6e3
  yng=1.0e6
  c=2.0e2
  phi=15.0
  poisn=0.30
  visc=5.0

  dt=2.5e-5
  nstep=60000

  open(1,file='space.d')
  write(1,'(2e12.4)') -dp,l*10+dp
  write(1,'(2e12.4)') -dp,l*10+dp
  write(1,'(2e12.4)') -dp,l+dp
  close(1)

  open(2,file='material.d')
  write(2,'(i2)') 2
  write(2,'(i2,7e12.4)') 1,rho,yng,poisn,c,phi,visc,0.0
  write(2,'(i2,7e12.4)') 2,rho,yng,poisn,c,phi,visc,0.0
  close(2)

  open(3,file='control.d')
  write(3,'(a)') 'c200phi15_'
  write(3,'(2i3)') 0,0
  write(3,'(i10)') nstep
  write(3,'(e12.4)') dt
  write(3,'(e12.4)') 0.005
  write(3,'(i3,2e12.4)') 1,0.1,0.1 ! 0:no damp. 1:art. visc, 2:Rayleigh damp. 
  close(3) 

  maxnp=int(l*10/dp)*int(l*10/dp)*int(l/dp)
  allocate(pos(3,maxnp))
  allocate(fix(3,maxnp))
  allocate(mat(maxnp))
  np=0
  do i=1,int(l*10/dp+1)  ! kisoziban
    !x=-dp+i*dp
    x=float(i-1)*dp
    do j=1,int(l*10/dp+1)
      !y=-dp+j*dp
      y=float(j-1)*dp
      do k=1,int(h0/dp+1)
        !z=-dp*2+k*dp
        z=-dp*2+k*dp
        if ( z < dp*0.1 ) then
          np=np+1
          pos(1,np)=x
          pos(2,np)=y
          pos(3,np)=z
          fix(1,np)=1
          fix(2,np)=1
          fix(3,np)=1
          mat(np)=1
          endif
        end do
       end do
      end do  
       
     do i=1,int(l*10/dp+1)   !entyuu
    !x=-dp+i*dp
    x=float(i-1)*dp
    do j=1,int(l*10/dp+1)
      !y=-dp+j*dp
      y=float(j-1)*dp
      do k=1,int(h0/dp+1)
        !z=-dp*2+k*dp
        z=-dp*2+k*dp
        if ( z > dp*0.1 ) then
          r = (x-(l*5.0))**2 + (y-(l*5.0))**2       
          if ( r0**2 < r.AND.r < r1**2 ) then
                 np=np+1
                 pos(1,np)=x
                 pos(2,np)=y
                 pos(3,np)=z
                 fix(1,np)=1
                 fix(2,np)=1
                 fix(3,np)=1
                 mat(np)=1
           end if
         end if
        end do 
      end do
    end do

  do i=1,int(l*10/dp+1)   !kyousitai
    x=float(i-1)*dp
    do j=1,int(l*10/dp+1)
      y=float(j-1)*dp
      do k=1,int(h0/dp+1)
        z=-dp*2+k*dp
        if ( z > dp*0.1 ) then
          r = (x-(l*5.0))**2 + (y-(l*5.0))**2       
          if ( r <r0**2 ) then
                 np=np+1
                 pos(1,np)=x
                 pos(2,np)=y
                 pos(3,np)=z
                 fix(1,np)=0
                 fix(2,np)=0
                 fix(3,np)=0
                 mat(np)=2
           end if
         end if
      end do
    end do
   end do

     do i=1,int(l*10/dp+1)   !saikaban
    x=float(i-1)*dp
    do j=1,int(l*10/dp+1)
      y=float(j-1)*dp
      do k=1,int(h0/dp+1)
        z=-dp*2+k*dp
       do m=1,int(h1/dp+1)
        z1=-dp*2+k*dp
        if ( z > dp*0.1 ) then
          r = (x-(l*5.0))**2 + (y-(l*5.0))**2       
          if (  z1 > z.AND. r <r0**2 ) then
                 np=np+1
                 pos(1,np)=x
                 pos(2,np)=y
                 pos(3,np)=z1
                 fix(1,np)=1
                 fix(2,np)=1
                 fix(3,np)=1
                 mat(np)=1
           end if
         end if
        end do 
      end do
    end do
   end do

  open(4,file='particle.d')
  write(4,'(i8)') np
  write(4,'(e12.4)') dp
  write(4,'(e12.4)') hsml
  do i=1,np
    write(4,'(i8,3e12.4,3i3,i3)') i,pos(1,i),pos(2,i),pos(3,i),&
      & fix(1,i),fix(2,i),fix(3,i),mat(i)
  end do
  close(4)
end program model_column
