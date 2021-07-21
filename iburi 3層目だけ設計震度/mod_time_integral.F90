module mod_time_integral
  !$use omp_lib
  use mod_type,only : vec3d
contains
  subroutine update_vel_hht(np,dt,alpha,acc_pre,acc,vel,xfix,yfix,zfix)
    type(vec3d),intent(in) :: acc_pre(:),acc(:)
    type(vec3d),intent(inout) :: vel(:)
    integer,intent(in) :: xfix(:),yfix(:),zfix(:)
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    real(8),intent(in) :: alpha
    integer :: i
    real(8) :: gamma,gamma_dt

    !$omp do private(gamma,gamma_dt)
    do i=1,np
       gamma=3.0d0/2.0d0-alpha
       gamma_dt=gamma*dt
       if ( xfix(i).eq.1 ) then
          vel(i)%x=0.0d0
       else
          vel(i)%x=vel(i)%x+(dt-gamma_dt)*acc_pre(i)%x+gamma_dt*acc(i)%x
       end if
       if ( yfix(i).eq.1 ) then
          vel(i)%y=0.0d0
       else
          vel(i)%y=vel(i)%y+(dt-gamma_dt)*acc_pre(i)%y+gamma_dt*acc(i)%y
       end if
       if ( zfix(i).eq.1 ) then
          vel(i)%z=0.0d0
       else
          vel(i)%z=vel(i)%z+(dt-gamma_dt)*acc_pre(i)%z+gamma_dt*acc(i)%z
       end if
    end do
    !$omp end do
  end subroutine update_vel_hht

  subroutine update_velocity_velret(np,dt,acc_pre,acc,vel,pos,dsp,xfix,yfix,zfix)
    type(vec3d),intent(in) :: acc_pre(:),acc(:)
    type(vec3d),intent(inout) :: vel(:),pos(:),dsp(:)
    integer,intent(in) :: xfix(:),yfix(:),zfix(:)
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    real(8) :: dt2
    integer :: i
    dt2=dt*dt*0.5d0
#ifdef NO_UPDATE_POS
    !$omp do
    do i=1,np
       if ( xfix(i) == 1 ) then
          vel(i)%x = 0.0d0
       else
          dsp(i)%x = dsp(i)%x + dt*vel(i)%x + dt2*acc_pre(i)%x
          vel(i)%x = vel(i)%x + 0.5d0*dt*(acc(i)%x+acc_pre(i)%x)
       end if
       if ( yfix(i) == 1 ) then
          vel(i)%y = 0.0d0
       else
          dsp(i)%y = dsp(i)%y + dt*vel(i)%y + dt2*acc_pre(i)%y
          vel(i)%y = vel(i)%y + 0.5d0*dt*(acc(i)%y+acc_pre(i)%y)
       end if
       if ( zfix(i) == 1 ) then
          vel(i)%z = 0.0d0
       else
          dsp(i)%z = dsp(i)%z + dt*vel(i)%z + dt2*acc_pre(i)%z
          vel(i)%z = vel(i)%z + 0.5d0*dt*(acc(i)%z+acc_pre(i)%z)
       end if
    end do
    !$omp end do
#else
    !$omp do
    do i=1,np
       if ( xfix(i) == 1 ) then
          vel(i)%x = 0.0d0
       else
          pos(i)%x = pos(i)%x + dt*vel(i)%x + dt2*acc_pre(i)%x
          dsp(i)%x = dsp(i)%x + dt*vel(i)%x + dt2*acc_pre(i)%x
          vel(i)%x = vel(i)%x + 0.5d0*dt*(acc(i)%x+acc_pre(i)%x)
       end if
       if ( yfix(i) == 1 ) then
          vel(i)%y = 0.0d0
       else
          pos(i)%y = pos(i)%y + dt*vel(i)%y + dt2*acc_pre(i)%y
          dsp(i)%y = dsp(i)%y + dt*vel(i)%y + dt2*acc_pre(i)%y
          vel(i)%y = vel(i)%y + 0.5d0*dt*(acc(i)%y+acc_pre(i)%y)
       end if
       if ( zfix(i) == 1 ) then
          vel(i)%z = 0.0d0
       else
          pos(i)%z = pos(i)%z + dt*vel(i)%z + dt2*acc_pre(i)%z
          dsp(i)%z = dsp(i)%z + dt*vel(i)%z + dt2*acc_pre(i)%z
          vel(i)%z = vel(i)%z + 0.5d0*dt*(acc(i)%z+acc_pre(i)%z)
       end if
    end do
    !$omp end do
#endif
  end subroutine update_velocity_velret
  
  subroutine update_vel(np,dt,acc,vel)
    type(vec3d),intent(in) :: acc(:)
    type(vec3d),intent(inout) :: vel(:)
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    integer :: i
    !$omp do
    do i=1,np
      vel(i)%x=vel(i)%x+acc(i)%x*dt
      vel(i)%y=vel(i)%y+acc(i)%y*dt
      vel(i)%z=vel(i)%z+acc(i)%z*dt
    end do
    !$omp end do
  end subroutine update_vel

  subroutine update_pos_hht(np,dt,xsp,ysp,zsp,alpha,acc_pre,acc,vel,pos,xfix,yfix,zfix,dsp,avel)
    type(vec3d),intent(in) :: acc_pre(:),acc(:),avel(:)
    type(vec3d),intent(inout) :: vel(:),pos(:),dsp(:)
    integer,intent(in) :: xfix(:),yfix(:),zfix(:)
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: alpha
    integer :: i
    real(8) :: beta
    real(8) :: dt2

#ifdef NO_UPDATE_POS
    ! nothing to do
#else
    !$omp do private(beta,dt2)
    do i=1,np
      beta=(2.0d0-alpha)**2/4.0d0
      dt2=dt*dt
      if ( xfix(i).eq.0 ) then
        pos(i)%x=pos(i)%x+vel(i)%x*dt+(0.5d0-beta)*dt2*acc_pre(i)%x+beta*dt2*acc(i)%x &
        &  + avel(i)%x*dt
      end if
      if ( pos(i)%x.lt.xsp(1) ) then
        pos(i)%x=xsp(1)
        vel(i)%x=0.0d0
      endif
      if ( pos(i)%x.gt.xsp(2) ) then
        pos(i)%x=xsp(2)
        vel(i)%x=0.0d0
      endif
      if (yfix(i).eq.0) then
        pos(i)%y=pos(i)%y+vel(i)%y*dt+(0.5d0-beta)*dt2*acc_pre(i)%y+beta*dt2*acc(i)%y &
        &  + avel(i)%y*dt
      end if
      if ( pos(i)%y.lt.ysp(1) ) then
        pos(i)%y=ysp(1)
        vel(i)%y=0.0d0
      endif
      if ( pos(i)%y.gt.ysp(2) ) then
        pos(i)%y=ysp(2)
        vel(i)%y=0.0d0
      endif
      if (zfix(i).eq.0) then
        pos(i)%z=pos(i)%z+vel(i)%z*dt+(0.5d0-beta)*dt2*acc_pre(i)%z+beta*dt2*acc(i)%z &
        &  + avel(i)%z*dt
      end if
      if ( pos(i)%z.lt.zsp(1) ) then
        pos(i)%z=zsp(1)
        vel(i)%z=0.0d0
      endif
      if ( pos(i)%z.gt.zsp(2) ) then
        pos(i)%z=zsp(2)
        vel(i)%z=0.0d0
      endif
      if ( xfix(i).eq.0 ) then
        dsp(i)%x=dsp(i)%x+vel(i)%x*dt+(0.5d0-beta)*dt2*acc_pre(i)%x+beta*dt2*acc(i)%x &
        &  + avel(i)%x*dt
      end if
      if ( yfix(i).eq.0 ) then
        dsp(i)%y=dsp(i)%y+vel(i)%y*dt+(0.5d0-beta)*dt2*acc_pre(i)%y+beta*dt2*acc(i)%y &
        &  + avel(i)%y*dt
      endif
      if ( zfix(i).eq.0 ) then
        dsp(i)%z=dsp(i)%z+vel(i)%z*dt+(0.5d0-beta)*dt2*acc_pre(i)%z+beta*dt2*acc(i)%z &
        &  + avel(i)%z*dt
      endif
    end do
    !$omp end do
#endif
  end subroutine update_pos_hht

  subroutine update_pos(np,dt,xsp,ysp,zsp,pos,vel,xfix,yfix,zfix,dsp,avel)
    type(vec3d),intent(in) :: avel(:)
    type(vec3d),intent(inout) :: pos(:),vel(:),dsp(:)
    integer,intent(in) :: xfix(:),yfix(:),zfix(:)
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    integer :: i

#ifdef NO_UPDATE_POS
    ! nothing to do
#else
    !$omp do
    do i=1,np
      if ( xfix(i).eq.0 ) then
        pos(i)%x=pos(i)%x+vel(i)%x*dt + avel(i)%x*dt
      else
        vel(i)%x=0.0d0
      end if
      if ( pos(i)%x.lt.xsp(1) ) then
        pos(i)%x=xsp(1)
        vel(i)%x=0.0d0
      endif
      if ( pos(i)%x.gt.xsp(2) ) then
        pos(i)%x=xsp(2)
        vel(i)%x=0.0d0
      endif
      if (yfix(i).eq.0) then
        pos(i)%y=pos(i)%y+vel(i)%y*dt + avel(i)%y*dt
      else
        vel(i)%y=0.0d0
      end if
      if ( pos(i)%y.lt.ysp(1) ) then
        pos(i)%y=ysp(1)
        vel(i)%y=0.0d0
      endif
      if ( pos(i)%y.gt.ysp(2) ) then
        pos(i)%y=ysp(2)
        vel(i)%y=0.0d0
      endif
      if (zfix(i).eq.0) then
        pos(i)%z=pos(i)%z+vel(i)%z*dt + avel(i)%z*dt
      else
        vel(i)%z=0.0d0
      end if
      if ( pos(i)%z.lt.zsp(1) ) then
        pos(i)%z=zsp(1)
        vel(i)%z=0.0d0
      endif
      if ( pos(i)%z.gt.zsp(2) ) then
        pos(i)%z=zsp(2)
        vel(i)%z=0.0d0
      endif
      if ( xfix(i).eq.0 ) then
        dsp(i)%x=dsp(i)%x+vel(i)%x*dt + avel(i)%x*dt
      end if
      if ( yfix(i).eq.0 ) then
        dsp(i)%y=dsp(i)%y+vel(i)%y*dt + avel(i)%y*dt
      endif
      if ( zfix(i).eq.0 ) then
        dsp(i)%z=dsp(i)%z+vel(i)%z*dt + avel(i)%z*dt
      endif
    end do
    !$omp end do
#endif
  end subroutine update_pos

end module mod_time_integral
