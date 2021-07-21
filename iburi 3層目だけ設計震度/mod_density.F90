module mod_density
  !$use omp_lib
  use mod_type,only : vec3d
  !use mod_box
  use mod_kernel
  implicit none
contains
    subroutine con_density2(np,ams,pos,vel,drho,npair,pair_list,pair_dwdx,rho)
    integer,intent(in) :: np
    real(8),intent(in) :: ams(:)
    type(vec3d),intent(in) :: pos(:),vel(:)
    integer,intent(in) :: npair(:),pair_list(:,:)
    real(8),intent(in) :: pair_dwdx(:,:,:),rho(:)
    real(8),intent(inout) :: drho(:)
    integer :: ip,jp,j
    real(8) :: amsj
    real(8) :: vcc
    type(vec3d) :: velip,veljp,posip,posjp

    !$omp do private(ip,velip,posip,j,jp,amsj,veljp,posjp,vcc)
    do ip=1,np
       velip=vel(ip)
       posip=pos(ip)
       drho(ip)=0.0d0
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          amsj=ams(jp)
          veljp=vel(jp)
          posjp=pos(jp)
          vcc = (velip%x-veljp%x)*pair_dwdx(1,j,ip) &
               & +(velip%y-veljp%y)*pair_dwdx(2,j,ip)+(velip%z-veljp%z)*pair_dwdx(3,j,ip)
          drho(ip) = drho(ip) + amsj*vcc/rho(jp)
       end do
       drho(ip) = drho(ip) * rho(ip)
    end do
    !$omp end do
  end subroutine con_density2

  subroutine con_density(np,hsml,xsp,ysp,zsp,ams,pos,vel,nbox,np_in_box,box,drho)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: ams(:)
    type(vec3d),intent(in) :: pos(:),vel(:)
    integer,intent(in) :: nbox(:),np_in_box(:),box(:,:)
    real(8),intent(inout) :: drho(:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r,amsj
    real(8) :: dx(3),w,dwdx(3)
    real(8) :: vcc
    type(vec3d) :: velip,veljp,posip,posjp
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2

    !$omp do
    do ip=1,np
       drho(ip)=0.0d0
    end do
    !$omp end do

    !$omp do private(ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,w,dwdx,amsj,vcc,velip,veljp,posip,posjp,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       velip=vel(ip)
       posip=pos(ip)
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
                      amsj=ams(jp)
                      veljp=vel(jp)
                      posjp=pos(jp)
                      dx(1)=posip%x-posjp%x
                      dx(2)=posip%y-posjp%y
                      dx(3)=posip%z-posjp%z
                      r=dsqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                      if ( r < (kh*hsml) ) then
!                         call cubic_spline_kernel(r,dx,hsml,w,dwdx)
                         dwdx = cubic_spline_kernel_dwdx(r,dx,hsml)
                         vcc = (velip%x-veljp%x)*dwdx(1)+(velip%y-veljp%y)*dwdx(2)+(velip%z-veljp%z)*dwdx(3)
                         drho(ip) = drho(ip) + amsj*vcc
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do
    !$omp end do
  end subroutine con_density

  subroutine sum_density(np,hsml,xsp,ysp,zsp,ams,pos,nbox,np_in_box,box,rho)
    integer,intent(in) :: np
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:)
    type(vec3d),intent(in) :: pos(:)
    integer,intent(in) :: nbox(:),np_in_box(:),box(:,:)
    real(8),intent(inout) :: rho(:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r,r2
    real(8) :: dx(3),w,dwdx(3),wsum(np)
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2
    !$omp do private(ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,w,dwdx,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       r=0.0d0
       w=cubic_spline_kernel_w(r,hsml)
       wsum(ip)=w*ams(ip)/rho(ip)
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
                            w=cubic_spline_kernel_w(r,hsml)
                            wsum(ip)=wsum(ip)+w*ams(jp)/rho(jp)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do
    !$omp end do

    !$omp do private(ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,w,dwdx,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       r=0.0d0
       dx=0.0d0
       w=cubic_spline_kernel_w(r,hsml)
       rho(ip)=w*ams(ip)
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
                            w = cubic_spline_kernel_w(r,hsml)
                            rho(ip)=rho(ip)+ams(jp)*w
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
       !rho(ip)=rho(ip)/wsum(ip)
    end do
    !$omp end do
  end subroutine sum_density

  subroutine sum_density2(np,hsml,ams,npair,pair_list,pair_w,wsum,rho)
    implicit none
    integer,intent(in) :: np,npair(:),pair_list(:,:)
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:),pair_w(:,:)
    real(8),intent(inout) :: rho(:)
    real(8),intent(out) :: wsum(:)
    integer :: ip,jp,j
    real(8) :: w
    !$omp do private(j,jp,w)
    do ip=1,np
       w=cubic_spline_kernel_w(0.0d0,hsml)
       wsum(ip)=w*ams(ip)/rho(ip)
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          wsum(ip)=wsum(ip)+pair_w(j,ip)*ams(jp)/rho(jp)
       end do
    end do
    !$omp end do

    !$omp do private(j,jp,w)
    do ip=1,np
       w=cubic_spline_kernel_w(0.0d0,hsml)
       rho(ip)=ams(ip)*w
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          rho(ip)=rho(ip)+ams(jp)*pair_w(j,ip)
       end do
    end do
    !$omp end do
  end subroutine sum_density2

  subroutine update_density(np,dt,rho,drho)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: dt
    real(8),intent(in) :: drho(:)
    real(8),intent(inout) :: rho(:)
    integer :: ip
    !$omp do
    do ip=1,np
       rho(ip)=rho(ip)+drho(ip)*dt
    end do
    !$omp end do
  end subroutine update_density

  subroutine particle_density(np,hsml,xsp,ysp,zsp,ams,pos,nbox,np_in_box,box,rho)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:)
    type(vec3d),intent(in) :: pos(:)
    integer,intent(in) :: nbox(:),np_in_box(:),box(:,:)
    real(8),intent(inout) :: rho(:)
    integer :: ip,jp
    integer :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    integer :: m
    real(8) :: r,r2
    real(8) :: dx(3),w,dwdx(3),wsum(np)
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2
    !$omp do private(ip,jp,ix,iy,iz,jx,jy,jz,kx,ky,kz,m,dx,r,w,dwdx,indx1)
    do ip=1,np
       ix=int((pos(ip)%x-xsp(1))/(hsml*kh))+2
       iy=int((pos(ip)%y-ysp(1))/(hsml*kh))+2
       iz=int((pos(ip)%z-zsp(1))/(hsml*kh))+2
       r=0.0d0
       w=cubic_spline_kernel_w(r,hsml)
       wsum(ip)=w
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
                            w=cubic_spline_kernel_w(r,hsml)
                            wsum(ip)=wsum(ip)+w
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
       rho(ip)=ams(ip)*wsum(ip)
    end do
    !$omp end do
  end subroutine particle_density

  subroutine particle_density2(np,hsml,ams,npair,pair_list,pair_w,rho)
    implicit none
    integer,intent(in) :: np,npair(:),pair_list(:,:)
    real(8),intent(in) :: hsml
    real(8),intent(in) :: ams(:),pair_w(:,:)
    real(8),intent(inout) :: rho(:)
    integer :: ip,jp,j
    real(8) :: wsum
    !$omp do private(j,jp,wsum)
    do ip=1,np
       wsum=cubic_spline_kernel_w(0.0d0,hsml)
       do j=1,npair(ip)
          jp=pair_list(j,ip)
          wsum=wsum+pair_w(j,ip)
       end do
       rho(ip)=wsum*ams(ip)
    end do
    !$omp end do
  end subroutine particle_density2
end module mod_density
