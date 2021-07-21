module mod_xsph
  !use mod_particle
  !use mod_box
  use mod_type
  use mod_kernel
  implicit none
contains
  subroutine xsph(np,hsml,xsp,ysp,zsp,ams,rho,avel,vel,pos,nbox,np_in_box,box)
    integer,intent(in) :: np
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:),rho(:)
    type(vec3d),intent(in) :: vel(:),pos(:)
    type(vec3d),intent(inout) :: avel(:)
    integer,intent(in) :: nbox(:), np_in_box(:),box(:,:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r
    real(8) :: dx(3),w,dvx(3),wsum
    real(8),parameter :: eps=0.3d0
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)+2
    nby=nbox(2)+2
    nbz=nbox(3)+2
    !$omp do private(ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,w,dvx,wsum,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       avel(ip)%x=0.0d0
       avel(ip)%y=0.0d0
       avel(ip)%z=0.0d0
       r=0.0d0
       dx=0.0d0
       w=cubic_spline_kernel_w(r,hsml)
       wsum=w*ams(ip)/rho(ip)
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
                            dvx(1)=vel(ip)%x-vel(jp)%x
                            dvx(2)=vel(ip)%y-vel(jp)%y
                            dvx(3)=vel(ip)%z-vel(jp)%z
                            w=cubic_spline_kernel_w(r,hsml)
                            avel(ip)%x = avel(ip)%x - ams(jp)*dvx(1)/(rho(jp))*w
                            avel(ip)%y = avel(ip)%y - ams(jp)*dvx(2)/(rho(jp))*w
                            avel(ip)%z = avel(ip)%z - ams(jp)*dvx(3)/(rho(jp))*w
                            wsum=wsum+w*ams(jp)/rho(jp)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
       avel(ip)%x=avel(ip)%x*0.3d0/wsum
       avel(ip)%y=avel(ip)%y*0.3d0/wsum
       avel(ip)%z=avel(ip)%z*0.3d0/wsum
    end do
    !$omp end do
  end subroutine xsph

  subroutine xsph2(np,hsml,ams,rho,avel,vel,pos,npair,pair_list,pair_w)
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:),rho(:)
    integer,intent(in) :: npair(:),pair_list(:,:)
    real(8),intent(in) :: pair_w(:,:)
    type(vec3d),intent(in) :: vel(:),pos(:)
    type(vec3d),intent(inout) :: avel(:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r
    real(8) :: dx(3),w,dvx(3),wsum
    real(8),parameter :: eps=0.3d0

    !$omp do private(ip,jp,m,dx,r,w,dvx,wsum)
    do ip=1,np
       w=cubic_spline_kernel_w(0.0d0,hsml)
       wsum=w*ams(ip)/rho(ip)
       do m=1,npair(ip)
          jp=pair_list(m,ip)
          avel(ip)%x=0.0d0
          avel(ip)%y=0.0d0
          avel(ip)%z=0.0d0
          r=0.0d0
          dx=0.0d0
          dx(1)=pos(ip)%x-pos(jp)%x
          dx(2)=pos(ip)%y-pos(jp)%y
          dx(3)=pos(ip)%z-pos(jp)%z
          dvx(1)=vel(ip)%x-vel(jp)%x
          dvx(2)=vel(ip)%y-vel(jp)%y
          dvx(3)=vel(ip)%z-vel(jp)%z
          w=pair_w(m,ip)
          avel(ip)%x = avel(ip)%x - ams(jp)*dvx(1)/(rho(jp))*w
          avel(ip)%y = avel(ip)%y - ams(jp)*dvx(2)/(rho(jp))*w
          avel(ip)%z = avel(ip)%z - ams(jp)*dvx(3)/(rho(jp))*w
          wsum=wsum+w*ams(jp)/rho(jp)
       end do
       avel(ip)%x=avel(ip)%x*0.3d0/wsum
       avel(ip)%y=avel(ip)%y*0.3d0/wsum
       avel(ip)%z=avel(ip)%z*0.3d0/wsum
    end do
    !$omp end do
  end subroutine xsph2


end module mod_xsph
