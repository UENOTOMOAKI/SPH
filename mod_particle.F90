module mod_particle
  !$use omp_lib
  use mod_const,only : max_pair
  use mod_simpara
  use mod_type, only : vec3d, tensor
  use mod_utils,only : getUnit
  use mod_vtk
  use mod_box
  use mod_extforce
  use mod_time_integral
  use mod_density
  use mod_xsph
  use mod_kernel

  implicit none

  type matprop
     real(8) yng,poisn,c1,phi1,c2,phi2,psi,vsound,fluid_viscosity,beta_coef
  end type matprop

  type(vec3d),allocatable :: pos(:)
  type(vec3d),allocatable :: vel(:)
  type(vec3d),allocatable :: acc(:)
  type(vec3d),allocatable :: acc_pre(:)
  type(vec3d),allocatable :: dsp(:)
  type(vec3d),allocatable :: intf0(:)
  type(vec3d),allocatable :: intf1(:)
  type(matprop),allocatable :: mat_prop(:)
  real(8),allocatable :: stress(:,:),dev_stress(:,:),solid_stress(:,:),fluid_stress(:,:)
  real,allocatable :: fluidized_time(:)
  real(8),allocatable :: pressure(:)
  logical,allocatable :: if_solid(:)
  real(8),allocatable :: adpstrain(:)
  integer,allocatable :: xfix(:),yfix(:),zfix(:)
  integer,allocatable :: id(:)
  integer,allocatable :: mat(:)
  real(8),allocatable :: rho(:),drho(:),ams(:),rho0(:)
  integer,allocatable :: istate(:)
  integer  :: nbox(3)
  integer,allocatable :: np_in_box(:)
  integer,allocatable :: box(:,:)

  type(vec3d),allocatable :: avel(:)  ! for XPSH

  type(tensor),allocatable :: dvdx(:)
  real(8),allocatable :: deps(:,:)

  integer,allocatable :: npair(:)
  integer,allocatable :: pair_list(:,:)
  real(8),allocatable :: pair_w(:,:)
  real(8),allocatable :: pair_dwdx(:,:,:)
  real(8),allocatable :: wsum(:)

contains

   real(8) function calc_beta(ip,current_time,coef)
      real(8), intent(in) :: current_time, coef
      integer, intent(in) :: ip
      real(8) :: beta
      beta = dexp(-coef*(current_time-fluidized_time(ip)))
   end function

  subroutine run_solid()
    integer :: np
    real(8) :: alpha,beta
    integer :: anamode  ! 0=new, 1=continued
    integer :: eqmode ! input eq. motion: 0=no,1=yes
    real(8) :: dp
    real(8) :: xsp(2),ysp(2),zsp(2)
    real(8) :: hsml
    character(len=100) :: ofhead
    real(4),allocatable :: xeqm(:),yeqm(:),zeqm(:)
    integer :: nyfix
    integer :: i,j,k,it,it0,d,ip
    real(4) :: t0,t1,t2
    logical :: lflt = .FALSE.
    integer :: damping = 0  ! 0 = no damping, 1 = artifical viscosity, 2 = Rayleigh damping
    integer :: ounit

    call input_control_params(anamode,eqmode,ofhead,damping,alpha,beta)
    if ( damping == 1 ) then
       write(*,*) '# Artifical damping is applied.'
    else if ( damping == 2 ) then
       write(*,*) '# Rayleigh damping is applied.'
    else
       write(*,*) '# No damping is applied.'
    end if

#ifdef NO_UPDATE_POS
    write(*,*) '# Position is not updated.'
#endif

#ifdef NO_GRAVITY
    write(*,*) '# No gravity load.'
#endif

    write(*,*) '# PAIR LIST Algorithm is used.'

    call input_space(xsp,ysp,zsp)
    call input_particles(np,dp,hsml)
    call input_materials(np)

    ! check whether 2d or 3d
    nyfix=count(mask=yfix==1)
    write(*,*) 'nyfix=',nyfix,np
    if ( nyfix == np ) then
       dim = 2
       write(*,*) 'Two dimensional analysis.'
    end if

    call set_mass(np,dp)
    call check_courant_condition(np,dt,dp)
    call output_space_vtk(xsp,ysp,zsp)
    call init_box(xsp,ysp,zsp,hsml,nbox,np_in_box,box)
    call vtk_out_box(hsml,xsp,ysp,zsp,nbox)

    if ( eqmode == 1 ) then
       !call readeqm('Kobe_NS.d',0.01,nstep,dt,xeqm)
       !call readeqm('Kobe_EW.d',0.01,nstep,dt,yeqm)
       !call readeqm('Kobe_UD.d',0.01,nstep,dt,zeqm)
       call readeqm('ricker.eqm',1.00,nstep,dt,xeqm)
       call readeqm('zero.eqm',1.00,nstep,dt,yeqm)
       call readeqm('zero.eqm',1.00,nstep,dt,zeqm)
    else if ( eqmode == 2 ) then
       call readeqm_kobe(nstep,dt,xeqm,yeqm,zeqm)
    else if ( eqmode == 3 ) then
       call readeqm('random.eqm',1.00,nstep,dt,xeqm)
       call zeroeqm(nstep,yeqm)
       call zeroeqm(nstep,zeqm)
    else if ( eqmode == 4 ) then
       call readeqm('sine600gal.eqm',1.00,nstep,dt,xeqm)
       call zeroeqm(nstep,yeqm)
       call zeroeqm(nstep,zeqm)
    end if

    do ip=1,np
       acc(ip)%x=0.0d0
       acc(ip)%y=0.0d0
       acc(ip)%z=0.0d0
       vel(ip)%x=0.0d0
       vel(ip)%y=0.0d0
       vel(ip)%z=0.0d0
       dsp(ip)%x=0.0d0
       dsp(ip)%y=0.0d0
       dsp(ip)%z=0.0d0
       avel(ip)%x=0.0d0
       avel(ip)%y=0.0d0
       avel(ip)%z=0.0d0
       if_solid(ip) = .TRUE.
       adpstrain(ip)=0.0d0
       istate(ip)=0
    end do
    dev_stress(:,:)=0.0d0
    it0=0
    ostep=0

    if ( anamode /= 0 ) then
       call read_binary_data(np,it0,ostep)
    end if

    ounit=getUnit()
    open(ounit,file='time_step.bi',form='unformatted')
    write(ounit) nstep,np

    call vtkout(np,0,ofhead)

    !$omp parallel default(none) &
    !$omp & shared(it0,nstep,np,acc_pre,acc,xsp,ysp,zsp,hsml,pos,&
    !$omp & lflt,dt,xeqm,yeqm,zeqm,eqmode,nout,ostep,ofhead,&
    !$omp & nbox,np_in_box,box,drho,ams,rho,vel,xfix,yfix,zfix,avel,dsp,&
    !$omp & npair,pair_list,pair_w,pair_dwdx,damping,alpha,beta,wsum,ounit,if_solid)
    do it=it0+1,it0+nstep

       !$omp do
       do i=1,np
          acc_pre(i)%x=acc(i)%x
          acc_pre(i)%y=acc(i)%y
          acc_pre(i)%z=acc(i)%z
       end do
      !$omp end do

       call set_box(np,xsp,ysp,zsp,hsml,pos,nbox,np_in_box,box)

       call make_interaction(np,hsml,xsp,ysp,zsp)

!       !$omp single
!       lflt = .FALSE.
!       if ( mod(it,5).eq.0 ) lflt = .TRUE.
!       !$omp end single

       call con_density2(np,ams,pos,vel,drho,npair,pair_list,pair_dwdx,rho)
       !call sum_density2(np,hsml,ams,npair,pair_list,pair_w,wsum,rho)
       !call particle_density2(np,hsml,ams,npair,pair_list,pair_w,rho)
       !call kernel_correction(np,ams,rho,pos,npair,pair_list,pair_dwdx)

       call calc_pressure(np)

       call intforce(np,dt,hsml,xsp,ysp,zsp,lflt,it)

       !if ( it*dt < 3.0d0 ) if_solid(:)=.TRUE.

#ifdef NO_GRAVITY
       ! nothing to do
#else
       call add_gravity(np,acc)
#endif

       if ( eqmode > 0 ) then
          call add_eqmotion(np,xeqm(it-it0),yeqm(it-it0),zeqm(it-it0),acc)
       end if

       if ( damping == 1 ) then
          call artvisc2(np,hsml,alpha,beta,0.01d0)
       else if ( damping == 2 ) then
          call rayleigh_damping(np,dt,alpha,beta)
       end if

#ifdef XY_PLANE
       call xy_plane(0.0d0,np,ams,dt)
#endif
       !call update_vel_hht(np,dt,0.70d0,acc_pre,acc,vel,xfix,yfix,zfix)
       !call xsph(np,hsml,xsp,ysp,zsp,ams,rho,avel,vel,pos,nbox,np_in_box,box)
       !call update_pos_hht(np,dt,xsp,ysp,zsp,0.70d0,acc_pre,acc,vel,pos,xfix,yfix,zfix,dsp,avel)
       call xsph2(np,hsml,ams,rho,avel,vel,pos,npair,pair_list,pair_w)
       call update_velocity_velret(np,dt,acc_pre,acc,vel,pos,dsp,xfix,yfix,zfix)

      call update_density(np,dt,rho,drho)

       !$omp do
         do ip=1,np
            if ( pos(ip)%x < xsp(1) ) pos(ip)%x = xsp(1)
            if ( pos(ip)%x > xsp(2) ) pos(ip)%x = xsp(2)
            if ( pos(ip)%y < ysp(1) ) pos(ip)%y = ysp(1)
            if ( pos(ip)%y > ysp(2) ) pos(ip)%y = ysp(2)
            if ( pos(ip)%z < zsp(1) ) pos(ip)%z = zsp(1)
            if ( pos(ip)%z > zsp(2) ) pos(ip)%z = zsp(2)
         end do
       !$omp end do

       !$omp single
       if ( mod(it,nout).eq.0 ) then
          write(*,'(a,i8)') ' STEP:',it
          write(*,*) 'max. of npair =',maxval(npair)
          ostep=ostep+1
          call vtkout(np,ostep,ofhead)
          call output_binary_data_time_step(ounit,it*dt)
       end if
       !$omp end single
    end do
    !$omp end parallel
    call output_binary_data(np,nstep,ostep)
    close(ounit)

  end subroutine run_solid

#ifdef XY_PLANE
  subroutine xy_plane(z,np,ams,dt)
    real(8),intent(in) :: z
    integer,intent(in) :: np
    real(8),intent(in) :: ams(:),dt
    integer :: ip
    real(8) :: er=0.0d0 ! repulsion coefficient
    real(8) :: fx,fy,fz,cd,r
    real(8) :: mu=0.6d0  ! friction coefficient
    cd = 1.0d0/(dt*0.25d0)
    !$omp do private(fx,fy,fz,r)
    do ip=1,np
       if ( pos(ip)%z < z ) then
          fz = -ams(ip)*((1.0+er)*vel(ip)%z-cd*(z-pos(ip)%z))
          acc(ip)%z = acc(ip)%z + fz
          fx = -cd*vel(ip)%x
          fy = -cd*vel(ip)%y
          r = dsqrt(fx*fx+fy*fy)/dabs(mu*fz)
          if ( r > 1.0d0 ) then
             fx = fx/r
             fy = fy/r
          end if
          acc(ip)%x = acc(ip)%x + fx
          acc(ip)%y = acc(ip)%y + fy
       end if
    end do
    !$omp end do
  end subroutine xy_plane
#endif

  subroutine alloc_particles(np)
    integer :: np
    allocate(pos(np))
    allocate(vel(np))
    allocate(acc(np))
    allocate(acc_pre(np))
    allocate(intf0(np))
    allocate(intf1(np))
    allocate(dsp(np))
    allocate(xfix(np))
    allocate(yfix(np))
    allocate(zfix(np))
    allocate(id(np))
    allocate(mat(np))
    allocate(rho(np))
    allocate(drho(np))
    allocate(rho0(np))
    allocate(ams(np))
    allocate(adpstrain(np))
    allocate(istate(np))
    allocate(pressure(np))
    allocate(if_solid(np))
    allocate(mat_prop(np))
    allocate(stress(6,np))
    allocate(dev_stress(6,np))
    allocate(solid_stress(6,np))
    allocate(fluid_stress(6,np))
    allocate(fluidized_time(np))
    allocate(dvdx(np))
    allocate(deps(6,np))
    allocate(avel(np))
    allocate(npair(np),pair_list(max_pair,np))
    allocate(pair_w(max_pair,np),pair_dwdx(3,max_pair,np))
    allocate(wsum(np))
  end subroutine alloc_particles

  subroutine input_space(xsp,ysp,zsp)
    real(8),intent(out) :: xsp(2),ysp(2),zsp(2)
    integer :: sunit
    sunit=getUnit()
    open(sunit,file='space.d')
    read(sunit,*) xsp(1),xsp(2)
    read(sunit,*) ysp(1),ysp(2)
    read(sunit,*) zsp(1),zsp(2)
    close(sunit)
  end subroutine input_space

  subroutine input_control_params(anamode,eqmode,ofhead,damping,alpha,beta)
    integer,intent(out) :: anamode,eqmode
    real(8),intent(out) :: alpha,beta
    integer,intent(out) :: damping ! 1 = Artificial viscosity, 2 = Rayleigh damping
    character(*) :: ofhead
    integer :: cunit
    real(8) :: dtout
    cunit=getUnit()
    open(cunit,file='control.d',form='formatted')
    read(cunit,*) ofhead
    write(*,*) ofhead
    read(cunit,*) anamode,eqmode
    read(cunit,*) nstep
    read(cunit,*) dt
    read(cunit,*) dtout
    read(cunit,*) damping,alpha,beta
    close(cunit)
    nout=int(dtout/dt)
    if ( nout.eq.0 ) stop
  end subroutine input_control_params

  subroutine input_particles(np,dp,hsml)
    integer,intent(out) :: np
    real(8),intent(out) :: dp
    real(8),intent(out) :: hsml
    integer :: i
    integer :: iunit
    iunit=getUnit()
    open(iunit,file='particle.d',form='formatted',status='old')
    read(iunit,*) np
    read(iunit,*) dp
    read(iunit,*) hsml

    call alloc_particles(np)

    do i=1,np
       read(iunit,*) id(i),pos(i)%x,pos(i)%y,pos(i)%z,xfix(i),yfix(i),zfix(i),&
         mat(i)
    end do
    close(iunit)

    !call output_model_vtk(np,id,pos,mat,ams)

  end subroutine input_particles

  subroutine input_materials(np)
    integer,intent(in) :: np
    real(8),allocatable :: r(:),y(:),p(:),cc1(:),p1(:),visc(:),beta_coef(:)
    integer :: nmat
    integer :: m
    integer :: i,j
    integer :: munit
    munit=getUnit()
    open(munit,file='material.d')
    read(munit,*) nmat
    allocate(r(nmat),y(nmat),p(nmat),cc1(nmat),p1(nmat),visc(nmat),beta_coef(nmat))
    do j=1,nmat
       read(munit,*) i,r(i),y(i),p(i),cc1(i),p1(i),visc(i),beta_coef(i)
    end do
    do i=1,np
       m=mat(i)
       rho(i)=r(m)
       rho0(i)=r(m)
       mat_prop(i)%yng=y(m)
       mat_prop(i)%poisn=p(m)
       mat_prop(i)%c1=cc1(m)
       mat_prop(i)%phi1=p1(m)
       mat_prop(i)%fluid_viscosity=visc(m)
       mat_prop(i)%beta_coef=beta_coef(m)
    end do
    deallocate(r,y,p,cc1,p1,visc)
    close(munit)
  end subroutine input_materials

  subroutine check_courant_condition(np,dt,dp)
    integer,intent(in) :: np
    real(8),intent(in) :: dt,dp
    integer :: ip
    real(8) :: dt0,dt1,sound_speed,g,k
    dt0 = dt*1000.0d0
    do ip=1,np
       g=0.5d0*mat_prop(ip)%yng/(1.0d0+mat_prop(ip)%poisn)
       k=mat_prop(ip)%yng/(3.0d0*(1.0d0-2.0*mat_prop(ip)%poisn))
       sound_speed=dsqrt(4.0d0*g/(3.0d0*rho(ip))+k/rho(ip))
       dt1 = 0.2d0*dp/sound_speed
       if ( dt0 > dt1 ) dt0=dt1
    end do
    if ( dt0 < dt ) then
       write(*,*) '*** Time increment dt must be smaller than ',dt0
       write(*,*) '*** However dt is set to be ',dt
       stop
    end if
  end subroutine check_courant_condition

  subroutine set_mass(np,dp)
    integer,intent(in) :: np
    !    real(8),intent(in) :: rho(max_np),dp
    real(8),intent(in) :: dp
    integer :: i
    if ( dim == 2 ) then
       do i=1,np
          ams(i)=dp*dp*rho(i)
       end do
    else if ( dim == 3 ) then
       do i=1,np
          ams(i)=dp*dp*dp*rho(i)
       end do
    end if
  end subroutine set_mass

  subroutine read_binary_data(np,istep,ostep)
    implicit none
    integer,intent(in) :: np
    integer,intent(out) :: istep,ostep
    integer :: i,j,ounit,np0
    ounit=getUnit()
    open(ounit,file='data.bi',form='UNFORMATTED')
    read(ounit) np0
    read(ounit) istep
    read(ounit) ostep
    read(ounit) pos !(pos(i)%x,i=1,np),(pos(i)%y,i=1,np),(pos(i)%z,i=1,np)
    read(ounit) vel !(xvel(i),i=1,np),(yvel(i),i=1,np),(zvel(i),i=1,np)
    read(ounit) acc !(xacc(i),i=1,np),(yacc(i),i=1,np),(zacc(i),i=1,np)
    read(ounit) stress !((stress(j,i),j=1,6),i=1,np)
    read(ounit) adpstrain !(adpstrain(i),i=1,np)
    read(ounit) id !(id(i),i=1,np)
    read(ounit) istate !(istate(i),i=1,np)
    read(ounit) if_solid
    read(ounit) fluidized_time
    close(ounit)
    if ( np /= np0 ) stop 'data.bi: Number of particle is different!'
  end subroutine read_binary_data

  subroutine vtkout(np,ostep,ofhead)
    integer,parameter :: maxlen=100
    integer,intent(in) :: np,ostep
    character(*),intent(in) :: ofhead
    character(len=maxlen) s_buffer
    character(len=maxlen) vtkf
    character*1 head
    real*8 cc,rr
    integer :: j
    integer, allocatable :: dmy(:)
    integer :: ip
    integer :: vunit
    real(8) :: var(np)
    write(s_buffer,'(i10.10)') ostep
    vtkf=trim(ofhead)//trim(s_buffer)//'.vtk'
    vunit=getUnit()
    call vtk_unstructured_header(vunit,vtkf)
    call vtk_unstructured_points(vunit,np,pos)
    call vtk_unstructured_cells(vunit,np)
    call vtk_data_point(vunit,np)
    call vtk_var_scalar_int(vunit,np,'id',id)
    call vtk_var_scalar_int(vunit,np,'mat',mat)
    do ip=1,np
      if ( if_solid(ip) ) then
        istate(ip) = 0
      else
        istate(ip) = 1
      end if
    end do
    call vtk_var_scalar_int(vunit,np,'istate',istate)
    call vtk_var_scalar(vunit,np,'adpstrain',adpstrain)
    call vtk_var_scalar(vunit,np,'pressure',pressure)
    do ip=1,np
       var(ip)=stress(1,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_xx',var)
    do ip=1,np
       var(ip)=stress(2,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_yy',var)
    do ip=1,np
       var(ip)=stress(3,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_zz',var)
    do ip=1,np
       var(ip)=stress(4,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_xy',var)
    do ip=1,np
       var(ip)=stress(5,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_yz',var)
    do ip=1,np
       var(ip)=stress(6,ip)
    end do
    call vtk_var_scalar(vunit,np,'sig_zx',var)
    call vtk_var_scalar(vunit,np,'density',rho)
    call vtk_var_scalar(vunit,np,'wsum',wsum)
    call vtk_var_vector(vunit,np,'displacement',dsp)
    call vtk_var_vector(vunit,np,'velocity',vel)
    call vtk_var_vector(vunit,np,'acceleration',acc)
    call vtk_var_vector(vunit,np,'averaged_velocity',avel)
    !call vtk_var_tensor(vunit,np,'stress',stress)
    close(vunit)
    !deallocate(var)
  end subroutine vtkout


  subroutine output_binary_data(np,istep,ostep)
    integer,intent(in) :: np,istep,ostep
    integer :: i,j,ounit
    ounit=getUnit()
    open(ounit,file='data.bi',form='UNFORMATTED')
    write(ounit) np
    write(ounit) istep
    write(ounit) ostep
    write(ounit) pos !(pos(i)%x,i=1,np),(pos(i)%y,i=1,np),(pos(i)%z,i=1,np)
    write(ounit) vel !(xvel(i),i=1,np),(yvel(i),i=1,np),(zvel(i),i=1,np)
    write(ounit) acc !(xacc(i),i=1,np),(yacc(i),i=1,np),(zacc(i),i=1,np)
    write(ounit) stress !((stress(j,i),j=1,6),i=1,np)
    write(ounit) adpstrain !(adpstrain(i),i=1,np)
    write(ounit) id !(id(i),i=1,np)
    write(ounit) istate !(istate(i),i=1,np)
    write(ounit) if_solid
    write(ounit) fluidized_time
    close(ounit)
  end subroutine output_binary_data

  subroutine output_binary_data_time_step(ounit,t)
    integer,intent(in) :: ounit
    real(8),intent(in) :: t
    integer :: i,j
    write(ounit) t
    write(ounit) pos !(pos(i)%x,i=1,np),(pos(i)%y,i=1,np),(pos(i)%z,i=1,np)
    write(ounit) vel !(xvel(i),i=1,np),(yvel(i),i=1,np),(zvel(i),i=1,np)
    write(ounit) acc !(xacc(i),i=1,np),(yacc(i),i=1,np),(zacc(i),i=1,np)
    write(ounit) dsp
    write(ounit) stress !((stress(j,i),j=1,6),i=1,np)
    write(ounit) adpstrain !(adpstrain(i),i=1,np)
    !write(ounit) id !(id(i),i=1,np)
    write(ounit) istate !(istate(i),i=1,np)
    write(ounit) if_solid
    write(ounit) fluidized_time
  end subroutine output_binary_data_time_step

  subroutine intforce(np,dt,hsml,xsp,ysp,zsp,lflt,it)
    implicit none
    integer,intent(in) :: np,it
    real(8),intent(in) :: dt,hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    logical,intent(in) :: lflt
    integer :: i,j

    !$omp do
    do i=1,np
       intf0(i)%x=intf1(i)%x
       intf0(i)%y=intf1(i)%y
       intf0(i)%z=intf1(i)%z
    end do
    !$omp end do

    call calc_dvdx2(np)

    call calc_deps(np,dt)

    !if ( lflt ) call shepard_filter_tensor(np,hsml,xsp,ysp,zsp,deps)

    !call calc_stress_two_phase_model(np,dt*it)
    call calc_stress_bingham_model(np,dt*it)

    !    if ( lflt ) call shepard_filter_tensor(np,hsml,xsp,ysp,zsp,stress)

    call calc_acc2(np)

    !$omp do
    do i=1,np
       intf1(i)%x=acc(i)%x
       intf1(i)%y=acc(i)%y
       intf1(i)%z=acc(i)%z
    end do
    !$omp end do

  end subroutine intforce

  subroutine calc_acc2(np)
    implicit none
    integer,intent(in) :: np
    integer :: ip,jp,j
    real(8) :: r,r2
    real(8) :: dx(3),w,dwdx(3),wsum,ri,rj,amsj,ri2_inv,rj2_inv
    integer :: nbx,nby,nbz,indx1

    !$omp do
    do ip=1,np
       acc(ip)%x=0.0d0
       acc(ip)%y=0.0d0
       acc(ip)%z=0.0d0
    end do
    !$omp end do

    !$omp do private(ip,jp,j,ri,rj,amsj,dwdx,ri2_inv,rj2_inv)
    do ip=1,np
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          dwdx(1)=pair_dwdx(1,j,ip)
          dwdx(2)=pair_dwdx(2,j,ip)
          dwdx(3)=pair_dwdx(3,j,ip)
          ri2_inv=1.0d0/(rho(ip)*rho(ip))
          rj2_inv=1.0d0/(rho(jp)*rho(jp))
          amsj=ams(jp)
          acc(ip)%x=acc(ip)%x &
               +amsj*((stress(1,ip)*ri2_inv+stress(1,jp)*rj2_inv)*dwdx(1) &
               +(stress(4,ip)*ri2_inv+stress(4,jp)*rj2_inv)*dwdx(2) &
               +(stress(6,ip)*ri2_inv+stress(6,jp)*rj2_inv)*dwdx(3))
          acc(ip)%y=acc(ip)%y &
               +amsj*((stress(2,ip)*ri2_inv+stress(2,jp)*rj2_inv)*dwdx(2) &
               +(stress(4,ip)*ri2_inv+stress(4,jp)*rj2_inv)*dwdx(1) &
               +(stress(5,ip)*ri2_inv+stress(5,jp)*rj2_inv)*dwdx(3))
          acc(ip)%z=acc(ip)%z &
               +amsj*((stress(3,ip)*ri2_inv+stress(3,jp)*rj2_inv)*dwdx(3) &
               +(stress(5,ip)*ri2_inv+stress(5,jp)*rj2_inv)*dwdx(2) &
               +(stress(6,ip)*ri2_inv+stress(6,jp)*rj2_inv)*dwdx(1))
       end do
    end do
    !$omp end do
  end subroutine calc_acc2

  subroutine calc_deps(np,dt)
    integer,intent(in) :: np
    real(8) :: dt
    integer :: ip

    !$omp do
    do ip=1,np
       deps(1,ip)=dvdx(ip)%xx
       deps(2,ip)=dvdx(ip)%yy
       deps(3,ip)=dvdx(ip)%zz
       deps(4,ip)=(dvdx(ip)%xy+dvdx(ip)%yx)*0.5d0
       deps(5,ip)=(dvdx(ip)%yz+dvdx(ip)%zy)*0.5d0
       deps(6,ip)=(dvdx(ip)%xz+dvdx(ip)%zx)*0.5d0
    end do
    !$omp end do

!    !$omp master
!    if ( dim == 2 ) write(*,*) 'deps_yy=',maxval(deps(3,:))
!    !$omp end master
  end subroutine calc_deps

  subroutine calc_sm_and_j2(np,stress,sm,dsj2)
    integer,intent(in) :: np
    real(8),intent(in) :: stress(:,:)
    real(8),intent(out) :: sm(:),dsj2(:)
    real(8) :: s1,s2,s3,tx,ty,tz,sx,sy,sz,aj2
    integer :: ip
    do ip=1,np
       s1=stress(1,ip)
       s2=stress(2,ip)
       s3=stress(3,ip)
       tx=stress(4,ip)
       ty=stress(5,ip)
       tz=stress(6,ip)
       sm(ip)=(s1+s2+s3)/3.0d0
       sx=s1-sm(ip)
       sy=s2-sm(ip)
       sz=s3-sm(ip)
       aj2=0.5d0*(sx*sx+sy*sy+sz*sz)+tx*tx+ty*ty+tz*tz
       dsj2(ip)=dsqrt(aj2)
    end do
  end subroutine calc_sm_and_j2

  ! subroutine calc_primary_stress(np,stress,s1,s2,s3)
  !   integer,intent(in) :: np
  !   real(8),intent(in) :: stress(:,:)
  !   real(8),intent(out) :: s1(np),s2(np),s3(np)
  !   integer :: ip,info
  !   real(8) :: sm,mat(3,3),w(3),work(9)
  !   !$omp do private(sm,mat,w,work,info)
  !   do ip=1,np
  !      sm = stress(1,ip) + stress(2,ip) + stress(3,ip)
  !      mat(1,1)=stress(1,ip)
  !      mat(1,2)=stress(4,ip)
  !      mat(1,3)=stress(6,ip)
  !      mat(2,2)=stress(2,ip)
  !      mat(2,3)=stress(5,ip)
  !      mat(3,3)=stress(3,ip)
  !      call dsyev('N','U',3,mat,3,w,work,9,info)
  !      if ( info /= 0 ) print *, 'info=',info
  !      s1(ip)=w(1)
  !      s2(ip)=w(2)
  !      s3(ip)=w(3)
  !   end do
  !   !$omp end do
  ! end subroutine

  subroutine calc_pressure(np)
    !compression positive
    integer, intent(in) :: np
    integer :: ip
    real(8) :: p1,e1,k1,c0
    !$omp do private(p1,e1,k1,c0)
    do ip=1,np
      p1=mat_prop(ip)%poisn
      e1=mat_prop(ip)%yng
      k1=e1/(3.0d0*(1.0d0-2.0*p1))
      c0=k1/rho0(ip)
      !pressure(ip)=c0*(rho(ip)-rho0(ip))
      pressure(ip) = c0*rho0(ip)/7.0d0*((rho(ip)/rho0(ip))**7-1.0d0)
    end do
    !$omp end do
  end subroutine calc_pressure

  subroutine calc_stress_two_phase_model(np,current_time)
    ! compression positive
    integer,intent(in) :: np
    real(8),intent(in) :: current_time
    integer :: ip,j
    real(8) :: c,fai,psai,p1,e1,k1
    real(8) :: alp1,alp2,g,s1,s2,s3
    real(8) :: tx,ty,tz,sm,sx,sy,sz
    real(8) :: aj2,dsj2,ai1,kk,fa
    real(8) :: dsig(6)
    real(8) :: depsrr,se,dlambda,rn
    real(8) :: omega(6)
    real(8) :: dpstrain(6),deps2
    real(8) :: pi180
    real(8) :: visc,visc0,tau_y,gamma
    real(8) :: beta,depsm,dev_deps(6)

    !$omp do private(pi180,p1,e1,k1,c,fai,psai,alp1,alp2,g,&
    !$omp & s1,s2,s3,tx,ty,tz,sm,sx,sy,sz,aj2,dsj2,ai1,kk,fa,depsrr,se,dsig,j,&
    !$omp & dlambda,omega,rn,visc,beta,visc0,tau_y,gamma,depsm,dev_deps)
    do ip=1,np
      pi180=pi/180.0d0
      p1=mat_prop(ip)%poisn
      e1=mat_prop(ip)%yng
      k1=e1/(3.0d0*(1.0d0-2.0*p1))
      c=mat_prop(ip)%c1
      fai=mat_prop(ip)%phi1*pi180
      alp1=2.0d0*dsin(fai)/(dsqrt(3.0d0)*(3.0d0-dsin(fai))) 
      if (if_solid(ip)) then
         g=0.5d0*e1/(1.0d0+p1)
         depsrr=deps(1,ip)+deps(2,ip)+deps(3,ip)
         dsig(1)=2.0d0*g*(deps(1,ip)-depsrr/3.0d0)
         dsig(2)=2.0d0*g*(deps(2,ip)-depsrr/3.0d0)
         dsig(3)=2.0d0*g*(deps(3,ip)-depsrr/3.0d0)
         dsig(4)=2.0d0*g*deps(4,ip)
         dsig(5)=2.0d0*g*deps(5,ip)
         dsig(6)=2.0d0*g*deps(6,ip)
         omega(4)=0.5d0*(dvdx(ip)%xy-dvdx(ip)%yx)
         omega(5)=0.5d0*(dvdx(ip)%yz-dvdx(ip)%zy)
         omega(6)=0.5d0*(dvdx(ip)%zx-dvdx(ip)%xz)
         ! Jaumann Stress Rate
         dsig(1) = dsig(1) + 2.0d0*omega(4)*dev_stress(4,ip) - 2.0d0*omega(6)*dev_stress(6,ip)
         dsig(2) = dsig(2) + 2.0d0*omega(5)*dev_stress(5,ip) - 2.0d0*omega(4)*dev_stress(4,ip)
         dsig(3) = dsig(3) + 2.0d0*omega(6)*dev_stress(6,ip) - 2.0d0*omega(5)*dev_stress(5,ip)
         dsig(4) = dsig(4) - omega(4)*(dev_stress(1,ip)-dev_stress(2,ip)) - omega(6)*dev_stress(5,ip) + omega(5)*dev_stress(6,ip)
         dsig(5) = dsig(5) - omega(5)*(dev_stress(2,ip)-dev_stress(3,ip)) - omega(4)*dev_stress(6,ip) + omega(6)*dev_stress(4,ip)
         dsig(6) = dsig(6) - omega(6)*(dev_stress(1,ip)-dev_stress(3,ip)) - omega(5)*dev_stress(4,ip) + omega(4)*dev_stress(5,ip)
         dev_stress(1,ip)=dev_stress(1,ip)+dsig(1)*dt
         dev_stress(2,ip)=dev_stress(2,ip)+dsig(2)*dt
         dev_stress(3,ip)=dev_stress(3,ip)+dsig(3)*dt
         dev_stress(4,ip)=dev_stress(4,ip)+dsig(4)*dt
         dev_stress(5,ip)=dev_stress(5,ip)+dsig(5)*dt
         dev_stress(6,ip)=dev_stress(6,ip)+dsig(6)*dt
         solid_stress(1,ip) = -pressure(ip) + dev_stress(1,ip)
         solid_stress(2,ip) = -pressure(ip) + dev_stress(2,ip)
         solid_stress(3,ip) = -pressure(ip) + dev_stress(3,ip)
         solid_stress(4,ip) =  dev_stress(4,ip)
         solid_stress(5,ip) =  dev_stress(5,ip)
         solid_stress(6,ip) =  dev_stress(6,ip)
         sm=(solid_stress(1,ip)+solid_stress(2,ip)+solid_stress(3,ip))/3.0d0
         sx=solid_stress(1,ip)-sm
         sy=solid_stress(2,ip)-sm
         sz=solid_stress(3,ip)-sm
         tx=solid_stress(4,ip)
         ty=solid_stress(5,ip)
         tz=solid_stress(6,ip)
         aj2=0.5d0*(sx*sx+sy*sy+sz*sz)+tx*tx+ty*ty+tz*tz
         dsj2=dsqrt(aj2)
         ai1=3.0d0*sm
         kk=6.0d0*c/(dsqrt(3.0d0)*(3.0d0-dsin(fai)))
         alp2=dsin(psai)
         fa=dsj2-alp1*ai1-kk
         if ( fa.gt.0.0d0 ) then
            !if ( current_time > 1.0 ) then
               if_solid(ip) = .FALSE.
               fluidized_time(ip) = current_time
            !endif
         end if
         do j=1,6
            stress(j,ip) = solid_stress(j,ip)
         end do
      else 
         visc0 = mat_prop(ip)%fluid_viscosity
         if ( pressure(ip) < 0.0d0 ) pressure(ip)=0.0d0
         tau_y = dabs(c + pressure(ip)*tan(fai))
         depsm = deps(1,ip) + deps(2,ip) + deps(3,ip)
         do j=1,6
            dev_deps(j) = deps(j,ip) - depsm/3.0d0
         end do
         gamma = dsqrt(2.0d0/3.0d0*(dev_deps(1)**2+dev_deps(2)**2+dev_deps(3)**2 + &
         &  2.0*(dev_deps(4)**2+dev_deps(5)**2+dev_deps(6)**2)))
         if ( gamma < 1.0d-3 ) gamma = 1.0d-3
         visc = visc0 + tau_y/gamma
         fluid_stress(1,ip) = -pressure(ip) + 4.0/3.0*visc*deps(1,ip)
         fluid_stress(2,ip) = -pressure(ip) + 4.0/3.0*visc*deps(2,ip)
         fluid_stress(3,ip) = -pressure(ip) + 4.0/3.0*visc*deps(3,ip)
         fluid_stress(4,ip) = 2.0*visc*deps(4,ip)
         fluid_stress(5,ip) = 2.0*visc*deps(5,ip)
         fluid_stress(6,ip) = 2.0*visc*deps(6,ip)
         !beta = calc_beta(ip,current_time,mat_prop(ip)%beta_coef)
         do j=1,6
            !stress(j,ip) = beta*solid_stress(j,ip) + (1.0d0-beta)*fluid_stress(j,ip)
            stress(j,ip) = fluid_stress(j,ip)
         end do
      endif
    end do
    !$omp end do
  end subroutine calc_stress_two_phase_model

  subroutine calc_stress_bingham_model(np,current_time)
   ! compression positive
   integer,intent(in) :: np
   real(8),intent(in) :: current_time
   integer :: ip,j
   real(8) :: c,fai
   real(8) :: pi180,gamma,tau_y
   real(8) :: visc0,visc

   pi180=pi/180.0d0
   !$omp do private(c,fai,j,visc0,visc,gamma,tau_y)
   do ip=1,np
      c=mat_prop(ip)%c1
      fai=mat_prop(ip)%phi1*pi180
      visc0 = mat_prop(ip)%fluid_viscosity
      !!if ( pressure(ip) < 0.0d0 ) pressure(ip)=0.0d0
      !if ( fai .lt. 0.001d0 ) fai = 0.1d0
      !if ( fai > 0.0d0 .and. pressure(ip) < -c/dtan(fai) ) then
      !   pressure(ip) = -c/dtan(fai)
      !   tau_y = 0.0
      !else if ( fai == 0.0d0 .and. pressure(ip) < 0.0d0 ) then
      !   pressure(ip) = 0.0d0
      !   tau_y = 0.0d0
      !else
      tau_y = c + pressure(ip)*dtan(fai)
      if ( tau_y < 0.0d0 ) tau_y = 0.0d0
      gamma = dsqrt(0.5d0*(deps(1,ip)**2+deps(2,ip)**2+deps(3,ip)**2 + &
      &  2.0d0*(deps(4,ip)**2+deps(5,ip)**2+deps(6,ip)**2)))
      if ( gamma < 1.0d-3 ) gamma = 1.0d-3
      visc = visc0 + tau_y/gamma
      if (visc > 1.0d8) visc=1.0d8
      stress(1,ip) = -pressure(ip) + visc*deps(1,ip)
      stress(2,ip) = -pressure(ip) + visc*deps(2,ip)
      stress(3,ip) = -pressure(ip) + visc*deps(3,ip)
      stress(4,ip) = visc*deps(4,ip)
      stress(5,ip) = visc*deps(5,ip)
      stress(6,ip) = visc*deps(6,ip)
   end do
   !$omp end do
end subroutine calc_stress_bingham_model

subroutine calc_dvdx2(np)
    integer,intent(in) :: np
    integer :: ip,jp,j
    real(8) :: rj_inv

    !$omp do
    do ip=1,np
       dvdx(ip)%xx=0.0d0
       dvdx(ip)%yy=0.0d0
       dvdx(ip)%zz=0.0d0
       dvdx(ip)%xy=0.0d0
       dvdx(ip)%yx=0.0d0
       dvdx(ip)%yz=0.0d0
       dvdx(ip)%zy=0.0d0
       dvdx(ip)%zx=0.0d0
       dvdx(ip)%xz=0.0d0
    end do
    !$omp end do

    !$omp do private(ip,jp,j,rj_inv)
    do ip=1,np
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          rj_inv=1.0d0/rho(jp)
          dvdx(ip)%xx=dvdx(ip)%xx-ams(jp)*(vel(ip)%x-vel(jp)%x)*pair_dwdx(1,j,ip)*rj_inv
          dvdx(ip)%xy=dvdx(ip)%xy-ams(jp)*(vel(ip)%x-vel(jp)%x)*pair_dwdx(2,j,ip)*rj_inv
          dvdx(ip)%xz=dvdx(ip)%xz-ams(jp)*(vel(ip)%x-vel(jp)%x)*pair_dwdx(3,j,ip)*rj_inv
          dvdx(ip)%yx=dvdx(ip)%yx-ams(jp)*(vel(ip)%y-vel(jp)%y)*pair_dwdx(1,j,ip)*rj_inv
          dvdx(ip)%yy=dvdx(ip)%yy-ams(jp)*(vel(ip)%y-vel(jp)%y)*pair_dwdx(2,j,ip)*rj_inv
          dvdx(ip)%yz=dvdx(ip)%yz-ams(jp)*(vel(ip)%y-vel(jp)%y)*pair_dwdx(3,j,ip)*rj_inv
          dvdx(ip)%zx=dvdx(ip)%zx-ams(jp)*(vel(ip)%z-vel(jp)%z)*pair_dwdx(1,j,ip)*rj_inv
          dvdx(ip)%zy=dvdx(ip)%zy-ams(jp)*(vel(ip)%z-vel(jp)%z)*pair_dwdx(2,j,ip)*rj_inv
          dvdx(ip)%zz=dvdx(ip)%zz-ams(jp)*(vel(ip)%z-vel(jp)%z)*pair_dwdx(3,j,ip)*rj_inv
       end do
    end do
    !$omp end do
  end subroutine calc_dvdx2

  subroutine shepard_filter_tensor(np,hsml,xsp,ysp,zsp,tensor1)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(inout) :: tensor1(:,:)
    real(8),allocatable :: tensor0(:,:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    integer :: d
    real(8) :: r,r2
    real(8) :: dx(3),w,dwdx(3),ri,rj,amsj
    real(8) :: dv(3)
    real(8) :: wsum
    integer :: nbx,nby,nbz,indx1
    !$omp single
    if ( .not.allocated(tensor0) ) allocate(tensor0(6,np))
    !$omp end single

    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2

    !$omp do
    do ip=1,np
       do d=1,6
          tensor0(d,ip)=tensor1(d,ip)
       end do
    end do
    !$omp end do

    !$omp do private(ix,iy,iz,r,dx,w,dwdx,d,jx,jy,jz,m,kx,ky,kz,jp,rj,amsj,wsum,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       r=0.0d0
       dx=0.0d0
!       call cubic_spline_kernel(r,dx,hsml,w,dwdx)
       w=cubic_spline_kernel_w(r,hsml)
       do d=1,6
          tensor1(d,ip)=ams(ip)*tensor0(d,ip)*w/rho(ip)
       end do
       wsum=ams(ip)*w/rho(ip)
       do jz=-1,1
          do jy=-1,1
             do jx=-1,1
                kx=ix+jx
                ky=iy+jy
                kz=iz+jz
                indx1=nbx*nby*(kz-1)+nbx*(ky-1)+kx
                if ( np_in_box(indx1) > 0 ) then
                   do m=1,np_in_box(indx1)
                      jp=box(m,indx1)
                      dx(1)=pos(ip)%x-pos(jp)%x
                      dx(2)=pos(ip)%y-pos(jp)%y
                      dx(3)=pos(ip)%z-pos(jp)%z
                      r=dsqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                      if ( r < (kh*hsml) ) then
                         !call cubic_spline_kernel(r,dx,hsml,w,dwdx)
                         w=cubic_spline_kernel_w(r,hsml)
                         rj=rho(jp)
                         amsj=ams(jp)
                         wsum=wsum+amsj*w/rj
                         do d=1,6
                            tensor1(d,ip)=tensor1(d,ip)+tensor0(d,jp)*w*amsj/rj
                         end do
                      end if
                   end do
                end if
             end do
          end do
       end do
       !write(*,*) wsum
       do d=1,6
          tensor1(d,ip)=tensor1(d,ip)/wsum
       end do
    end do
    !$omp end do
  end subroutine shepard_filter_tensor

  subroutine rayleigh_damping(np,dt,alpha,beta)
    integer,intent(in) :: np
    real(8),intent(in) :: alpha,beta,dt
    integer :: ip
    !$omp do
    do ip=1,np
      if ( if_solid(ip) ) then
        acc(ip)%x=acc(ip)%x-alpha*vel(ip)%x+beta*(intf1(ip)%x-intf0(ip)%x)/dt
        acc(ip)%y=acc(ip)%y-alpha*vel(ip)%y+beta*(intf1(ip)%y-intf0(ip)%y)/dt
        acc(ip)%z=acc(ip)%z-alpha*vel(ip)%z+beta*(intf1(ip)%z-intf0(ip)%z)/dt
      else
        acc(ip)%x=acc(ip)%x+beta*(intf1(ip)%x-intf0(ip)%x)/dt
        acc(ip)%y=acc(ip)%y+beta*(intf1(ip)%y-intf0(ip)%y)/dt
        acc(ip)%z=acc(ip)%z+beta*(intf1(ip)%z-intf0(ip)%z)/dt
      end if
    end do
    !$omp end do
  end subroutine rayleigh_damping

  subroutine artvisc(np,hsml,xsp,ysp,zsp,alpha,beta,etq)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: alpha,beta,etq
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r
    real(8) :: dx(3),w,dwdx(3),ri,rj,amsj,ri2_inv,rj2_inv
    real(8) :: vr,muv,mrho,piv
    real(8),parameter :: mc=1000.0d0
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2
    !$omp do private(ip,jp,jx,jy,jz,m,kx,ky,kz,dx,r,ri,rj,w,amsj,dwdx,ri2_inv,rj2_inv,vr,muv,mrho,piv,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       do jz=-1,1
          do jy=-1,1
             do jx=-1,1
                kx=ix+jx
                ky=iy+jy
                kz=iz+jz
                indx1=nbx*nby*(kz-1)+nbx*(ky-1)+kx
                if ( np_in_box(indx1) > 0 ) then
                   do m=1,np_in_box(indx1)
                      jp=box(m,indx1)
                      dx(1)=pos(ip)%x-pos(jp)%x
                      dx(2)=pos(ip)%y-pos(jp)%y
                      dx(3)=pos(ip)%z-pos(jp)%z
                      r=dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
                      if ( r < (kh*hsml) ) then
                         vr=(vel(ip)%x-vel(jp)%x)*dx(1)+(vel(ip)%y-vel(jp)%y)*dx(2)+(vel(ip)%z-vel(jp)%z)*dx(3)
                         if ( vr.lt.0.0d0 ) then
                            muv=hsml*vr/(r+hsml*hsml*etq*etq)
                            mrho=(rho(ip)+rho(jp))*0.5d0
                            piv=(beta*muv-alpha*mc)*muv/mrho
                            !call cubic_spline_kernel(r,dx,hsml,w,dwdx)
                            dwdx=cubic_spline_kernel_dwdx(r,dx,hsml)
                            acc(ip)%x=acc(ip)%x - ams(jp)*piv*dwdx(1)
                            acc(ip)%y=acc(ip)%y - ams(jp)*piv*dwdx(2)
                            acc(ip)%z=acc(ip)%z - ams(jp)*piv*dwdx(3)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do
    !$omp end do
  end subroutine artvisc

  subroutine artvisc2(np,hsml,alpha,beta,etq)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: hsml,alpha,beta,etq
    integer :: ip,jp,j
    real(8) :: r
    real(8) :: dx(3)
    real(8) :: vr,muv,mrho,piv,p1,e1,g,k1
    real(8) :: mc,cip,cjp

    !$omp do private(ip,jp,j,dx,r,vr,muv,mrho,piv)
    do ip=1,np
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          dx(1)=pos(ip)%x-pos(jp)%x
          dx(2)=pos(ip)%y-pos(jp)%y
          dx(3)=pos(ip)%z-pos(jp)%z
          r=dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
          vr=(vel(ip)%x-vel(jp)%x)*dx(1)+(vel(ip)%y-vel(jp)%y)*dx(2)+(vel(ip)%z-vel(jp)%z)*dx(3)
          if ( vr.lt.0.0d0 ) then
             muv=hsml*vr/(r+hsml*hsml*etq*etq)
             mrho=(rho(ip)+rho(jp))*0.5d0
             p1=mat_prop(ip)%poisn
             e1=mat_prop(ip)%yng
             k1=e1/(3.0d0*(1.0d0-2.0*p1))
             g=0.5d0*e1/(1.0d0+p1)
             cip=dsqrt((4.0*g/(3.0d0*rho(ip))+k1/rho(ip)))
             p1=mat_prop(jp)%poisn
             e1=mat_prop(jp)%yng
             k1=e1/(3.0d0*(1.0d0-2.0*p1))
             g=0.5d0*e1/(1.0d0+p1)
             cjp=dsqrt((4.0*g/(3.0d0*rho(jp))+k1/rho(jp)))
             mc=(cip+cjp)*0.5d0
             piv=(beta*muv-alpha*mc)*muv/mrho
             acc(ip)%x=acc(ip)%x - ams(jp)*piv*pair_dwdx(1,j,ip)
             acc(ip)%y=acc(ip)%y - ams(jp)*piv*pair_dwdx(2,j,ip)
             acc(ip)%z=acc(ip)%z - ams(jp)*piv*pair_dwdx(3,j,ip)
          end if
       end do
    end do
    !$omp end do
  end subroutine artvisc2

  subroutine make_interaction(np,hsml,xsp,ysp,zsp)
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8) :: dx(3),r
    integer :: ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m
    integer :: nbx,nby,nbz,indx1
    real(8) :: dwdx(3)
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2

    !$omp do private(ip,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,dwdx,jp,indx1)
    do ip=1,np
       npair(ip)=0
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
        do jz=-1,1
          do jy=-1,1
             do jx=-1,1
                kx=ix+jx
                ky=iy+jy
                kz=iz+jz
                indx1=nbx*nby*(kz-1)+nbx*(ky-1)+kx
                if ( np_in_box(indx1) > 0 ) then
                   do m=1,np_in_box(indx1)
                      jp=box(m,indx1)
                      if ( ip /= jp ) then
                         dx(1)=pos(ip)%x-pos(jp)%x
                         dx(2)=pos(ip)%y-pos(jp)%y
                         dx(3)=pos(ip)%z-pos(jp)%z
                         r=dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
                         if ( r < (kh*hsml) ) then
                            npair(ip)=npair(ip)+1
                            if ( npair(ip) > max_pair ) stop 'too much interaction'
                            pair_list(npair(ip),ip)=jp
                            pair_w(npair(ip),ip)=cubic_spline_kernel_w(r,hsml)
                            dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
                            pair_dwdx(1,npair(ip),ip) = dwdx(1)
                            pair_dwdx(2,npair(ip),ip) = dwdx(2)
                            pair_dwdx(3,npair(ip),ip) = dwdx(3)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do
    !$omp end do

    ! !$omp master
    ! write(*,*) maxval(npair,1)
    ! ip=109
    ! do jp=1,npair(ip)
    !    write(104,*) pos(pair_list(jp,ip))%x,pos(pair_list(jp,ip))%z
    ! end do
    ! stop
    ! !$omp end master

  end subroutine make_interaction

  subroutine verification_check_pairs1(np,hsml,ip)
    integer,intent(in) :: np,ip
    real(8),intent(in) :: hsml
    integer :: jp,j,kpair
    real(8) :: dst
    integer :: iunit

    write(*,*) 'ip=',ip
    write(*,*) 'npair=',npair(ip)
    write(*,*) (pair_list(j,ip),j=1,npair(ip))
    iunit = getUnit()
    open(iunit,file='test_pairs.d')
    do j=1,npair(ip)
       jp=pair_list(j,ip)
       write(iunit,*) pos(jp)%x,pos(jp)%y,pos(jp)%z
    end do
    close(iunit)
    kpair=0
    do jp=1,np
       if ( jp /= ip ) then
          dst = (pos(ip)%x-pos(jp)%x)**2 + (pos(ip)%y-pos(jp)%y)**2 + (pos(ip)%z-pos(jp)%z)**2
          dst = dsqrt(dst)
          if ( dst < kh*hsml ) then
             kpair=kpair+1
             write(*,*) kpair,jp,dst
          end if
       end if
    end do
  end subroutine verification_check_pairs1

  subroutine verification_check_pairs2(np,hsml)
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    integer :: ip,jp,kpair
    real(8) :: dst

    do ip=1,np
       kpair=0
       do jp=1,np
          if ( jp /= ip ) then
             dst = (pos(ip)%x-pos(jp)%x)**2 + (pos(ip)%y-pos(jp)%y)**2 + (pos(ip)%z-pos(jp)%z)**2
             dst = dsqrt(dst)
             if ( dst < kh*hsml ) then
                kpair=kpair+1
             end if
          end if
       end do
       if (npair(ip)/=kpair) then
          write(*,*) 'VERIFICATION TEST FAILURE'
          write(*,*) 'ERROR:',ip
       end if
    end do
  end subroutine verification_check_pairs2

end module mod_particle
