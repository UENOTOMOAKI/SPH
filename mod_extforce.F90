module mod_extforce
  !$use omp_lib
  use mod_utils
  use mod_const, only : grav
  use mod_vtk
  integer,parameter :: maxlen=100
contains

  subroutine add_eqmotion(np,xeqm,yeqm,zeqm,acc)
    type(vec3d),intent(inout) :: acc(:)
    integer,intent(in) :: np
    real(4),intent(in) :: xeqm,yeqm,zeqm
    integer :: ip
    !$omp do
    do ip=1,np
       acc(ip)%x=acc(ip)%x+dble(xeqm)
       acc(ip)%y=acc(ip)%y+dble(yeqm)
       acc(ip)%z=acc(ip)%z+dble(zeqm)
    end do
    !$omp end do
  end subroutine add_eqmotion
  
  subroutine add_gravity(np,acc)
    type(vec3d),intent(inout) :: acc(:)
    integer,intent(in) :: np
    integer :: i
    !$omp do
    do i=1,np
       acc(i)%z=acc(i)%z-grav
    end do
    !$omp end do
  end subroutine add_gravity

  subroutine readeqm(filename,factor,nstep,dt,eqm)
    character(*),intent(in) :: filename
    integer,intent(in) :: nstep
    real(8),intent(in) :: dt
    real(4),intent(in) :: factor
    real(4),allocatable,intent(out) :: eqm(:)
    real(4),allocatable :: eqm0(:)
    integer :: iunit
    integer :: nd
    real(4) :: tinc,dt1
    integer :: i
    integer,save :: count = 1
    character(2) :: s_buffer
    character(len=maxlen) :: ofname
    if(.not.allocated(eqm)) allocate(eqm(nstep))
    iunit = getUnit()
    dt1=real(dt)
    open(iunit,file=trim(filename))
    read(iunit,*)
    read(iunit,*) nd,tinc
    if(.not.allocated(eqm0)) allocate(eqm0(nd))
    do i=1,nd
       read(iunit,*) eqm0(i)
    end do
    close(iunit)
    call interp(eqm,eqm0,factor,nstep,dt1,nd,tinc)
    write(s_buffer,'(i2.2)') count
    iunit=getUnit()
    ofname='eqm'//trim(s_buffer)//'.d'
    open(iunit,file=ofname)
    do i=1,nd
       write(iunit,*) tinc*real(i),eqm0(i)
       if ( tinc*real(i) >= real(nstep*dt) ) exit
    end do
    close(iunit)
    count=count+1
  end subroutine readeqm

  subroutine zeroeqm(nstep,eqm)
    integer,intent(in) :: nstep
    real(4),allocatable,intent(out) :: eqm(:)
    if(.not.allocated(eqm)) allocate(eqm(nstep))
    eqm=0.0
  end subroutine zeroeqm

  SUBROUTINE INTERP(X,RA,AMPL,NSTEP,DELTA,ND,DT)
    implicit none
    real(4) :: X(*),RA(*),ampl,delta,dt
    real(4) ::  time,ds,ts,te,de
    integer nstep,nd
    integer i,n,n1
    RA(ND + 1) = RA(ND)
    DO I = 1,NSTEP
       TIME = DFLOAT(I) * DELTA
       N = TIME / DT
       N1 = N + 1
       IF(N > 0) GO TO 250
       DS = 0.d0
       TS = 0.d0
       GO TO 270
250    CONTINUE
       DS = AMPL * RA(N)
       TS = DFLOAT(N) * DT
270    CONTINUE
       DE = AMPL * RA(N1)
       TE = TS + DT
       X(I) = DS + (DE - DS) * (TIME - TS) /(TE - TS)
    END DO
  end SUBROUTINE INTERP

    subroutine readeqm_kobe(nstep,dt,xeqm,yeqm,zeqm)
    integer,intent(in) :: nstep
    real(8),intent(in) :: dt
    real(4) :: factor = 0.01
    real(4),allocatable,intent(out) :: xeqm(:),yeqm(:),zeqm(:)
    real(4),allocatable :: xeqm0(:),yeqm0(:),zeqm0(:)
    integer :: iunit
    integer :: nd
    real(4) :: tinc,dt1
    integer :: i
    if(.not.allocated(xeqm)) allocate(xeqm(nstep))
    if(.not.allocated(yeqm)) allocate(yeqm(nstep))
    if(.not.allocated(zeqm)) allocate(zeqm(nstep))
    
    iunit = getUnit()
    dt1=real(dt)
    open(iunit,file='Kobe_acc.eqm')
    read(iunit,*) nd,tinc
    if(.not.allocated(xeqm0)) allocate(xeqm0(nd))
    if(.not.allocated(yeqm0)) allocate(yeqm0(nd))
    if(.not.allocated(zeqm0)) allocate(zeqm0(nd))
    do i=1,nd
       read(iunit,*) xeqm0(i),yeqm0(i),zeqm0(i)
    end do
    close(iunit)

    call interp(xeqm,xeqm0,factor,nstep,dt1,nd,tinc)
    if ( dim == 2 ) then
       yeqm0(:) = 0.0
       yeqm(:) = 0.0d0
    end if
    call interp(yeqm,yeqm0,factor,nstep,dt1,nd,tinc)
    call interp(zeqm,zeqm0,factor,nstep,dt1,nd,tinc)

    iunit=getUnit()
    open(iunit,file='eqm.d')
    do i=1,nd
       write(iunit,*) tinc*real(i),xeqm0(i),yeqm0(i),zeqm0(i)
       if ( tinc*real(i) >= real(nstep*dt) ) exit
    end do
    close(iunit)
  end subroutine readeqm_kobe
end module mod_extforce
