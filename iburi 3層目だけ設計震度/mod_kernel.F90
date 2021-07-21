module mod_kernel
  use mod_const, only : pi,max_pair
  use mod_type, only : vec3d
  use mod_simpara
  use mod_utils,only : getUnit
  implicit none
  integer,parameter :: kh=2
contains
  subroutine cubic_spline_kernel(r,dx,hsml,w,dwdx)
    !$acc routine seq
    real(8),intent(in) :: r,dx(3),hsml
    real(8),intent(out) :: w,dwdx(3)
    real(8) :: q,factor
    integer :: d
    q = r/hsml 
    w = 0.d0 
    dwdx(:) = 0.0d0

    if ( dim == 2 ) then
       factor = 15.0d0/(7.0d0*pi*hsml*hsml)
    else if ( dim == 3 ) then
       factor = 3.0d0/(2.0d0*pi*hsml*hsml*hsml)
    end if

    if (q.ge.0.0d0.and.q.le.1.0d0) then 
       w = factor*(2.d0/3.d0-q*q + q**3/2.d0) 
       do d = 1,3
          dwdx(d) = factor*(-2.d0+3.d0/2.d0*q)/hsml**2*dx(d) 
       enddo
    else if (q.gt.1.0d0.and.q.le.2.0d0) then 
       w = factor*1.0d0/6.0d0*(2.d0-q)**3
       do d=1,3
          dwdx(d) =-factor*1.d0/6.d0*3.d0*(2.d0-q)**2/hsml*(dx(d)/r)
       enddo
    else 
       w=0.d0 
       do d=1,3 
          dwdx(d) = 0.d0
       enddo
    endif
  end subroutine cubic_spline_kernel

#ifdef TEST_CUBIC_SPLINE_KERNEL
  subroutine test_cubic_spline_kernel()
    real(8),parameter :: dp=1.0d0
    real(8) :: hsml
    integer :: nx,ny,nz
    integer :: ix,iy,iz,iunit
    real(8) :: xpos,ypos,zpos,r,w,dx(3),dwdx(3),wsum,dwsum(3)
    nx=11
    ny=11
    nz=1
    hsml=dp*1.3d0

    iunit=getUnit()
    open(iunit,file='test_cubic_cpline_kernel.d')
    
    ! for xy-plane
    do iz=-int(nz*0.5),int(nz*0.5)
       zpos=dp*(iz)
       do iy=-int(ny*0.5),int(ny*0.5)
          ypos=dp*(iy)
          do ix=-int(nx*0.5),int(nx*0.5)
             xpos=dp*(ix)
             dx(1)=xpos
             dx(2)=ypos
             dx(3)=zpos
             r = dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
             if ( r < hsml*kh ) then
                w = cubic_spline_kernel_w(r,hsml)
                dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
             else
                w=0.0d0
                dwdx(:)=0.0d0
             end if
             write(iunit,*) xpos,ypos,w,dwdx(1),dwdx(2)
          end do
          write(iunit,*)
       end do
    end do
    write(iunit,*)

    ! for yz-plane
    nx=1
    ny=11
    nz=11
    do ix=-int(nx*0.5),int(nx*0.5)
       xpos=dp*(ix)
       do iy=-int(ny*0.5),int(ny*0.5)
          ypos=dp*(iy)
          do iz=-int(nz*0.5),int(nz*0.5)
             zpos=dp*(iz)
             dx(1)=xpos
             dx(2)=ypos
             dx(3)=zpos
             r = dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
             if ( r < hsml*kh ) then
                w = cubic_spline_kernel_w(r,hsml)
                dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
             else
                w=0.0d0
                dwdx(:)=0.0d0
             end if
             write(iunit,*) ypos,zpos,w,dwdx(1),dwdx(2),dwdx(3)
          end do
          write(iunit,*)
       end do
    end do
    write(iunit,*)
    
    ! for zx-plane
    nx=11
    ny=1
    nz=11
    do iy=-int(ny*0.5),int(ny*0.5)
       ypos=dp*(iy)
       do iz=-int(nz*0.5),int(nz*0.5)
          zpos=dp*(iz)
          do ix=-int(nx*0.5),int(nx*0.5)
             xpos=dp*(ix)
             dx(1)=xpos
             dx(2)=ypos
             dx(3)=zpos
             r = dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
             if ( r < hsml*kh ) then
                w = cubic_spline_kernel_w(r,hsml)
                dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
             else
                w=0.0d0
                dwdx(:)=0.0d0
             end if
             write(iunit,*) zpos,xpos,w,dwdx(1),dwdx(2),dwdx(3)
          end do
          write(iunit,*)
       end do
    end do
    close(iunit)

    ! for wsum
    nx=11
    ny=11
    nz=11
    wsum=0.0d0
    dwsum(:)=0.0d0
    do iy=-int(ny*0.5),int(ny*0.5)
       ypos=dp*(iy)
       do iz=-int(nz*0.5),int(nz*0.5)
          zpos=dp*(iz)
          do ix=-int(nx*0.5),int(nx*0.5)
             xpos=dp*(ix)
             dx(1)=xpos
             dx(2)=ypos
             dx(3)=zpos
             r = dsqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
             if ( r < hsml*kh ) then
                w = cubic_spline_kernel_w(r,hsml)
                dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
             else
                w=0.0d0
                dwdx(:)=0.0d0
             end if
             wsum=wsum+w
             dwsum(1)=dwsum(1)*dwdx(1)
             dwsum(2)=dwsum(2)*dwdx(2)
             dwsum(3)=dwsum(3)*dwdx(3)
          end do
       end do
    end do
    write(*,*) 'wsum=',wsum
    write(*,*) 'dwsum=',dwsum
    
    stop "### TEST_CUBIC_SPLINE_KERNEL is done."
  end subroutine test_cubic_spline_kernel
#endif 

  function cubic_spline_kernel_w(r,hsml) result(w)
    implicit none
    real(8),intent(in) :: r,hsml
    real(8) :: w
    real(8) :: factor,q

    q=r/hsml
    
    if ( dim == 2 ) then
       factor = 15.0d0/(7.0d0*pi*hsml*hsml)
    else if ( dim == 3 ) then
       factor = 3.0d0/(2.0d0*pi*hsml*hsml*hsml)
    end if

    if (q.ge.0.0d0.and.q.le.1.0d0) then 
       w = factor*(2.d0/3.d0-q*q + q**3/2.d0) 
    else if (q.gt.1.0d0.and.q.le.2.0d0) then 
       w = factor*1.0d0/6.0d0*(2.d0-q)**3
    else 
       w=0.d0 
    endif
  end function cubic_spline_kernel_w

  function  cubic_spline_kernel_dwdx(r,dx,hsml) result(dwdx)
    real(8),intent(in) :: r,dx(3),hsml
    real(8) :: dwdx(3)
    real(8) :: q,factor
    integer :: d
    q = r/hsml 
    dwdx(:) = 0.0d0
    
    if ( dim == 2 ) then
       factor = 15.0d0/(7.0d0*pi*hsml*hsml)
    else if ( dim == 3 ) then
       factor = 3.0d0/(2.0d0*pi*hsml*hsml*hsml)
    end if

    if (q.ge.0.0d0.and.q.le.1.0d0) then 
       do d=1,3
          dwdx(d) = factor*(-2.d0+3.d0/2.d0*q)/hsml**2*dx(d) 
       enddo
    else if (q.gt.1.0d0.and.q.le.2.0d0) then 
       do d=1,3
          dwdx(d) =-factor*1.d0/6.d0*3.d0*(2.d0-q)**2/hsml*(dx(d)/r)
       enddo
    else 
       do d=1,3
          dwdx(d) = 0.d0
       enddo
    endif
  end function cubic_spline_kernel_dwdx

  subroutine kernel_correction(np,ams,rho,pos,npair,pair_list,pair_dwdx)
    integer,intent(in) :: np,npair(:),pair_list(:,:)
    real(8),intent(in) :: ams(:),rho(:)
    type(vec3d),intent(in) :: pos(:)
    real(8),intent(inout) :: pair_dwdx(:,:,:)
    integer :: ip,jp,j,k
    real(8) :: m(3,3),dx(3),lm(3,3,np),det

    if ( dim == 3 ) then
       !$omp do private(ip,jp,j,dx,m,det,k)
       do ip=1,np
          m(:,:)=0.0d0
          do j=1,npair(ip)
             jp=pair_list(j,ip)
             dx(1) = pos(ip)%x-pos(jp)%x
             dx(2) = pos(ip)%y-pos(jp)%y
             dx(3) = pos(ip)%z-pos(jp)%z
             m(1,1)=m(1,1)+pair_dwdx(1,j,ip)*dx(1)*ams(jp)/rho(jp)
             m(2,1)=m(2,1)+pair_dwdx(2,j,ip)*dx(1)*ams(jp)/rho(jp)
             m(3,1)=m(3,1)+pair_dwdx(3,j,ip)*dx(1)*ams(jp)/rho(jp)
             m(1,2)=m(1,2)+pair_dwdx(1,j,ip)*dx(2)*ams(jp)/rho(jp)
             m(2,2)=m(2,2)+pair_dwdx(2,j,ip)*dx(2)*ams(jp)/rho(jp)
             m(3,2)=m(3,2)+pair_dwdx(3,j,ip)*dx(2)*ams(jp)/rho(jp)
             m(1,3)=m(1,3)+pair_dwdx(1,j,ip)*dx(3)*ams(jp)/rho(jp)
             m(2,3)=m(2,3)+pair_dwdx(2,j,ip)*dx(3)*ams(jp)/rho(jp)
             m(3,3)=m(3,3)+pair_dwdx(3,j,ip)*dx(3)*ams(jp)/rho(jp)
          end do
          det = dabs(m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + m(1,3)*m(2,1)*m(3,2) &
               & - m(1,3)*m(2,2)*m(3,1) - m(1,2)*m(2,1)*m(3,3) - m(1,1)*m(2,3)*m(3,2))
          !       if ( det > 0.0001d0 ) then !if ( det > 0.01d0 ) then
          det = 1.0d0/det
          lm(1,1,ip) =  (m(2,2)*m(3,3)-m(2,3)*m(3,2))*det
          lm(1,2,ip) = -(m(2,1)*m(3,3)-m(2,3)*m(3,1))*det
          lm(1,3,ip) =  (m(2,1)*m(3,2)-m(2,2)*m(3,1))*det
          lm(2,1,ip) = -(m(1,2)*m(3,3)-m(1,3)*m(3,2))*det
          lm(2,2,ip) =  (m(1,1)*m(3,3)-m(1,3)*m(3,1))*det
          lm(2,3,ip) = -(m(1,1)*m(3,2)-m(1,2)*m(3,1))*det
          lm(3,1,ip) =  (m(1,2)*m(2,3)-m(1,3)*m(2,2))*det
          lm(3,2,ip) = -(m(1,1)*m(2,3)-m(1,3)*m(2,1))*det
          lm(3,3,ip) =  (m(1,1)*m(2,2)-m(1,2)*m(2,1))*det
          ! else
          !    do j=1,3
          !       do k=1,3
          !          if ( j == k ) then
          !             lm(j,k,ip) = 1.0d0
          !          else
          !             lm(j,k,ip) = 0.0d0
          !          end if
          !       end do
          !    end do
          ! end if
       end do
       !$omp end do
     
       !$omp do private(ip,j,dx)
       do ip=1,np
          do j=1,npair(ip)
             dx(1) = lm(1,1,ip)*pair_dwdx(1,j,ip) + lm(1,2,ip)*pair_dwdx(2,j,ip) + lm(1,3,ip)*pair_dwdx(3,j,ip)
             dx(2) = lm(2,1,ip)*pair_dwdx(1,j,ip) + lm(2,2,ip)*pair_dwdx(2,j,ip) + lm(2,3,ip)*pair_dwdx(3,j,ip)
             dx(3) = lm(3,1,ip)*pair_dwdx(1,j,ip) + lm(3,2,ip)*pair_dwdx(2,j,ip) + lm(3,3,ip)*pair_dwdx(3,j,ip)
             pair_dwdx(1,j,ip) = dx(1)
             pair_dwdx(2,j,ip) = dx(2)
             pair_dwdx(3,j,ip) = dx(3)
          end do
       end do
       !$omp end do
    else if ( dim == 2 ) then
       !$omp do private(ip,jp,j,dx,m,det,k)
       do ip=1,np
          m(:,:)=0.0d0
          do j=1,npair(ip)
             jp=pair_list(j,ip)
             dx(1) = pos(ip)%x-pos(jp)%x
             dx(3) = pos(ip)%z-pos(jp)%z
             m(1,1)=m(1,1)+pair_dwdx(1,j,ip)*dx(1)*ams(jp)/rho(jp)
             m(3,1)=m(3,1)+pair_dwdx(3,j,ip)*dx(1)*ams(jp)/rho(jp)
             m(1,3)=m(1,3)+pair_dwdx(1,j,ip)*dx(3)*ams(jp)/rho(jp)
             m(3,3)=m(3,3)+pair_dwdx(3,j,ip)*dx(3)*ams(jp)/rho(jp)
          end do
          det = dabs(m(1,1)*m(3,3) - m(1,3)*m(3,1))
          !       if ( det > 0.0001d0 ) then !if ( det > 0.01d0 ) then
          det = 1.0d0/det
          lm(1,1,ip) =  m(3,3)*det
          lm(1,3,ip) =  -m(1,3)*det
          lm(3,1,ip) =  -m(3,1)*det
          lm(3,3,ip) =  m(1,1)*det
       end do
       !$omp end do
     
       !$omp do private(ip,j,dx)
       do ip=1,np
          do j=1,npair(ip)
             dx(1) = lm(1,1,ip)*pair_dwdx(1,j,ip) + lm(1,3,ip)*pair_dwdx(3,j,ip)
             dx(3) = lm(3,1,ip)*pair_dwdx(1,j,ip) + lm(3,3,ip)*pair_dwdx(3,j,ip)
             pair_dwdx(1,j,ip) = dx(1)
             pair_dwdx(3,j,ip) = dx(3)
          end do
       end do
       !$omp end do
    end if
  end subroutine kernel_correction
end module mod_kernel
