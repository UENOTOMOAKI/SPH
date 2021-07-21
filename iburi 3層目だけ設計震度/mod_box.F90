module mod_box
  use mod_vtk
  use mod_kernel
  implicit none
  integer,parameter :: max_np_in_box=100
contains

  subroutine init_box(xsp,ysp,zsp,hsml,nbox,np_in_box,box)
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2),hsml
    integer :: i,j,k
    integer,intent(out) :: nbox(:)
    integer,allocatable,intent(out) :: np_in_box(:),box(:,:)
    integer :: nbx,nby,nbz
    nbox(1)=int((xsp(2)-xsp(1))/(hsml*kh))+3
    nbox(2)=int((ysp(2)-ysp(1))/(hsml*kh))+3
    nbox(3)=int((zsp(2)-zsp(1))/(hsml*kh))+3
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2
    allocate(np_in_box(nbx*nby*nbz))
    allocate(box(max_np_in_box,nbx*nby*nbz))
  end subroutine init_box

  subroutine set_box(np,xsp,ysp,zsp,hsml,pos,nbox,np_in_box,box)
    type(vec3d),intent(in) :: pos(:)
    integer,intent(in) :: np
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),intent(in) :: hsml
    integer,intent(in) :: nbox(:)
    integer,intent(out) :: np_in_box(:),box(:,:)
    integer ip,ix,iy,iz,m
    real(8) :: hsml_inv
    integer :: nbx,nby,nbz,indx1
    nbx=nbox(1)!+2
    nby=nbox(2)!+2
    nbz=nbox(3)!+2
    !$omp do private(ix,iy,iz,indx1)
    do iz=1,nbz
       do iy=1,nby
          do ix=1,nbx
             indx1=nbx*nby*(iz-1)+nbx*(iy-1)+ix
             np_in_box(indx1)=0
          end do
       end do
    end do
    !$omp end do

    !$omp do private(ix,iy,iz,m,indx1)
    do iz=1,nbz
       do iy=1,nby
          do ix=1,nbx
             indx1=nbx*nby*(iz-1)+nbx*(iy-1)+ix
             do m=1,max_np_in_box
                box(m,indx1)=0
             end do
          end do
       end do
    end do
    !$omp end do

    !$omp single
    hsml_inv=1.0d0/(hsml*kh)
    do ip=1,np
       !if ( ip.eq.19897 ) write(20,*) pos(ip)%x,pos(ip)%y,pos(ip)%z
       ix=int((pos(ip)%x-xsp(1))*hsml_inv)+2
       iy=int((pos(ip)%y-ysp(1))*hsml_inv)+2
       iz=int((pos(ip)%z-zsp(1))*hsml_inv)+2
       indx1=nbx*nby*(iz-1)+nbx*(iy-1)+ix
       np_in_box(indx1)=np_in_box(indx1)+1
       if ( np_in_box(indx1) > max_np_in_box ) then
          write(0,*) 'Too many particle are pushed into a box.'
          write(0,*) 'Particle ID=',ip
          stop
       end if
       box(np_in_box(indx1),indx1)=ip
    end do
    !$omp end single
  end subroutine set_box

  subroutine vtk_out_box(hsml,xsp,ysp,zsp,nbox)
    implicit none
    real(8),intent(in) :: hsml
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    real(8),allocatable :: x(:),y(:),z(:)
    integer,intent(in) :: nbox(:)
    integer :: nx,ny,nz,np,i,j,k
    integer :: vunit
    nx=nbox(1)+1
    ny=nbox(2)+1
    nz=nbox(3)+1
    np=nx*ny*nz
    allocate(x(np),y(np),z(np))
    np=1
    do k=1,nz
       do j=1,ny
          do i=1,nx
             x(np)=xsp(1)+hsml*kh*(i-2)
             y(np)=ysp(1)+hsml*kh*(j-2)
             z(np)=zsp(1)+hsml*kh*(k-2)
             np=np+1
          end do
       end do
    end do
    np=np-1
    vunit=getUnit()
    call vtk_structured_header(vunit,'box.vtk',nx,ny,nz)
    call vtk_structured_point(vunit,np,x,y,z)
    close(vunit)
    deallocate(x,y,z)
  end subroutine vtk_out_box

end module mod_box
