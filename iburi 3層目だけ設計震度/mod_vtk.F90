module mod_vtk
  use mod_type,only : vec3d
  use mod_const,only : end_rec,output_format
  use mod_utils,only : getUnit,str  
contains
  subroutine write_header_vtk_ascii(vunit,np1,np2)
    implicit none
    integer,intent(in) :: np1,np2
    integer,intent(in) :: vunit
    write(vunit,'(a)') '<?xml version="1.0"?>'
    write(vunit,'(a)') '<VTKFile type="UnstructuredGrid">'
    write(vunit,'(a)') '<UnstructuredGrid>'
    write(vunit,'(a,i8,a,i8,a)') '<Piece NumberOfPoints="',np1,'" NumberOfCells="',np2,'">'
  end subroutine write_header_vtk_ascii

  subroutine write_footer_vtk_ascii(vunit)
    implicit none
    integer,intent(in) :: vunit
    write(vunit,'(a)') '</Piece>'
    write(vunit,'(a)') '</UnstructuredGrid>'
    write(vunit,'(a)') '</VTKFile>'
  end subroutine write_footer_vtk_ascii

  subroutine write_cells_vtk_ascii(vunit,np)
    implicit none
    integer,intent(in) :: np,vunit
    integer :: i
    write(vunit,'(a)') '<Cells>'
    write(vunit,'(a)') '<DataArray Name="connectivity" format="ascii" type="Int32">'
    do i=1,np
       write(vunit,'(i8)') i-1
    end do
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '<DataArray Name="offsets" format="ascii" type="Int32">'
    do i=1,np
       write(vunit,'(i8)') i
    end do
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '<DataArray Name="types" format="ascii" type="UInt8">'
    do i=1,np
       write(vunit,'(i8)') 1
    end do
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '</Cells>'
  end subroutine write_cells_vtk_ascii

  subroutine write_header_vtk_binary(vunit,np)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: np,vunit
    character(len=maxlen) :: s_buffer
    write(vunit) '<?xml version="1.0"?>'//end_rec
    write(vunit) '<VTKFile type="UnstructuredGrid">'//end_rec
    write(vunit) '<UnstructuredGrid>'//end_rec
    s_buffer='<Piece NumberOfPoints="'//trim(str(np))//'" NumberOfCells="'//trim(str(np))//'">'
  end subroutine write_header_vtk_binary

  subroutine output_model_vtk(np,id,pos,mat,ams)
    type(vec3d),intent(in) :: pos(:)
    integer,intent(in) :: id(:),mat(:)
    real(8),intent(in) :: ams(:)
    integer,parameter :: maxlen=100
    integer,intent(in) :: np
    integer :: i,vunit
    character(len=maxlen) :: s_buffer
    vunit=getUnit()
    select case(trim(output_format))
    case('ASCII')
       open(vunit,file='input.vtu',form='formatted')
       call write_header_vtk_ascii(vunit,np,np)
       call write_points_vtk_ascii(vunit,np,pos)
       call write_cells_vtk_ascii(vunit,np)
       write(vunit,'(a)') '<CellData>'
       write(vunit,'(a)') '<DataArray Name="id" NumberOfComponents="1" format="ascii" type="Int32">'
       do i=1,np
          write(vunit,'(i8)') id(i)
       end do
       write(vunit,'(a)') '</DataArray>'
       write(vunit,'(a)') '<DataArray Name="mat" NumberOfComponents="1" format="ascii" type="Int32">'
       do i=1,np
          write(vunit,'(i8)') mat(i)
       end do
       write(vunit,'(a)') '</DataArray>'
       write(vunit,'(a)') '<DataArray Name="ams" NumberOfComponents="1" format="ascii" type="Float32">'
       do i=1,np
          write(vunit,'(f12.5)') ams(i)
       end do
       write(vunit,'(a)') '</DataArray>'
       write(vunit,'(a)') '</CellData>'
       call write_footer_vtk_ascii(vunit)
       close(vunit)

    case('BINARY')
       open(vunit,file='input.vtu',form='unformatted')
       call write_header_vtk_binary(vunit,np)
       close(vunit)
    endselect
  end subroutine output_model_vtk

  subroutine output_vtk(nout,np,pos)
    type(vec3d),intent(in) :: pos(:)
    integer,parameter :: maxlen=100
    integer,intent(in) :: nout
    integer,intent(in) :: np
    character*8 step
    integer :: i
    character(len=maxlen) :: s_buffer
    integer :: vunit

    write(step,'(i8.8)') nout
    vunit=getUnit()

    select case(trim(output_format))
    case('ASCII')
       open(vunit,file='particles'//trim(step)//'.vtu',form='formatted')
       call write_header_vtk_ascii(vunit,np,np)
       call write_points_vtk_ascii(vunit,np,pos)
       call write_cells_vtk_ascii(vunit,np)
       write(vunit,'(a)') '<CellData>'
       write(vunit,'(a)') '<DataArray Name="id" NumberOfComponents="1" format="ascii" type="Int32">'
       do i=1,np
          write(vunit,'(i8)') i
       end do
       write(vunit,'(a)') '</DataArray>'
       write(vunit,'(a)') '</CellData>'
       call write_footer_vtk_ascii(vunit)
       close(vunit)

    case('BINARY')
       open(vunit,file='particles'//trim(step)//'.vtu',form='unformatted')
       call write_header_vtk_binary(vunit,np)
       close(vunit)
    endselect
  end subroutine output_vtk

  subroutine output_space_vtk(xsp,ysp,zsp)
    real(8),intent(in) :: xsp(2),ysp(2),zsp(2)
    integer :: i,vunit
    vunit=getUnit()
    open(vunit,file='space.vtu')
    call write_header_vtk_ascii(vunit,8,1)
    write(vunit,'(a)') '<Points>'
    write(vunit,'(a)') '<DataArray NumberOfComponents="3" format="ascii" type="Float32">'
    write(vunit,'(3f10.3)') xsp(1),ysp(1),zsp(1)
    write(vunit,'(3f10.3)') xsp(2),ysp(1),zsp(1)
    write(vunit,'(3f10.3)') xsp(2),ysp(2),zsp(1)
    write(vunit,'(3f10.3)') xsp(1),ysp(2),zsp(1)
    write(vunit,'(3f10.3)') xsp(1),ysp(1),zsp(2)
    write(vunit,'(3f10.3)') xsp(2),ysp(1),zsp(2)
    write(vunit,'(3f10.3)') xsp(2),ysp(2),zsp(2)
    write(vunit,'(3f10.3)') xsp(1),ysp(2),zsp(2)
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '</Points>'
    write(vunit,'(a)') '<Cells>'
    write(vunit,'(a)') '<DataArray Name="connectivity" format="ascii" type="Int32">'
    write(vunit,'(8i2)') (i,i=0,7)
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '<DataArray Name="offsets" format="ascii" type="Int32">'
    write(vunit,'(i1)') 8
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '<DataArray Name="types" format="ascii" type="UInt8">'
    write(vunit,'(i2)') 12
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '</Cells>'
    write(vunit,'(a)') '<CellData>'
    write(vunit,'(a)') '<DataArray Name="id" NumberOfComponents="1" format="ascii" type="Int32">'
    write(vunit,'(i1)') 1
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '</CellData>'
    call write_footer_vtk_ascii(vunit)
    close(vunit)
  end subroutine output_space_vtk



  subroutine vtk_var_tensor(vtkf,np,varname,var)
    implicit none
    integer,intent(in) :: vtkf,np
    character(*),intent(in) :: varname
    real(8),intent(in) :: var(:,:)
    integer :: ip
    write(vtkf) 'TENSORS '//trim(varname)//' double'//end_rec
    write(vtkf) (&
         var(1,ip),var(4,ip),var(6,ip),&
         var(4,ip),var(2,ip),var(5,ip),&
         var(6,ip),var(5,ip),var(6,ip),ip=1,np)
    write(vtkf) end_rec
  end subroutine vtk_var_tensor

  subroutine vtk_var_scalar_int(vtkf,np,varname,var)
    implicit none
    integer,intent(in) :: vtkf,np
    character(*),intent(in) :: varname
    integer,intent(in) :: var(:)
    integer :: ip
    write(vtkf) 'SCALARS '//trim(varname)//' int 1'//end_rec
    write(vtkf) 'LOOKUP_TABLE default'//end_rec
    write(vtkf) var !(var(ip),ip=1,np)
    write(vtkf) end_rec
  end subroutine vtk_var_scalar_int

  subroutine vtk_var_scalar(vtkf,np,varname,var)
    implicit none
    integer,intent(in) :: vtkf,np
    character(*),intent(in) :: varname
    real(8),intent(in) :: var(:)
    integer :: ip
    write(vtkf) 'SCALARS '//trim(varname)//' double 1'//end_rec
    write(vtkf) 'LOOKUP_TABLE default'//end_rec
    write(vtkf) (var(ip),ip=1,np)
    write(vtkf) end_rec
  end subroutine vtk_var_scalar

  subroutine vtk_data_point(vtkf,np)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: vtkf,np
    character(len=maxlen) :: s_buffer
    write(s_buffer,'(i8)') np
    write(vtkf) 'POINT_DATA '//trim(s_buffer)//end_rec
  end subroutine vtk_data_point

  subroutine vtk_unstructured_points(vtkf,np,pos)
    type(vec3d),intent(in) :: pos(:)
    integer,parameter :: maxlen=100
    integer,intent(in) :: vtkf,np
    integer :: ip
    character(len=maxlen) :: s_buffer
    write(s_buffer,'(i8)') np
    write(vtkf) 'POINTS '//trim(s_buffer)//' double'//end_rec
    write(vtkf) (pos(ip)%x,pos(ip)%y,pos(ip)%z,ip=1,np)
    write(vtkf) end_rec
  end subroutine vtk_unstructured_points

  subroutine vtk_unstructured_cells(vtkf,np)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: vtkf,np
    integer :: ip
    character(len=maxlen) :: s_buffer
    integer :: connect(2,np),ctype(np)
    write(vtkf) 'CELLS '
    write(s_buffer,'(i8)') np
    write(vtkf) trim(s_buffer)
    write(s_buffer,'(i8)') np*2
    write(vtkf) trim(s_buffer)//end_rec
    do ip=1,np
       connect(1,ip)=1
       connect(2,ip)=ip-1
    end do
    write(vtkf) connect
    write(vtkf) end_rec
    write(vtkf) 'CELL_TYPES '
    write(s_buffer,'(i8)') np
    write(vtkf) trim(s_buffer)//end_rec
    do ip=1,np
       ctype(ip)=1
    end do
    write(vtkf) ctype
    write(vtkf) end_rec
  end subroutine vtk_unstructured_cells

    subroutine vtk_var_vector(vtkf,np,varname,var)
    implicit none
    integer,intent(in) :: vtkf,np
    character(*),intent(in) :: varname
    type(vec3d),intent(in) :: var(:)
    integer :: ip
    write(vtkf) 'VECTORS '//trim(varname)//' double'//end_rec
    write(vtkf) (var(ip)%x,var(ip)%y,var(ip)%z,ip=1,np)
    write(vtkf) end_rec
  end subroutine vtk_var_vector

  subroutine vtk_unstructured_header(vunit,vtkf)
    implicit none
    integer,intent(in) :: vunit
    character(*),intent(in) :: vtkf
    write(*,*) vtkf
    open(vunit,file=trim(vtkf),form='UNFORMATTED',access='STREAM',action='WRITE',status='REPLACE')
    write(vunit) '# vtk DataFile Version 3.0'//end_rec
    write(vunit) trim(vtkf)//end_rec
    write(vunit) 'BINARY'//end_rec
    write(vunit) 'DATASET UNSTRUCTURED_GRID'//end_rec
  end subroutine vtk_unstructured_header

  subroutine vtk_structured_point(vunit,np,x,y,z)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: vunit,np
    real(8),intent(in) :: x(*),y(*),z(*)
    character(len=maxlen) :: s_buffer
    integer :: ip
    write(s_buffer,'(i8)') np
    write(vunit) 'POINTS '//trim(s_buffer)//' double'//end_rec
    write(vunit) (x(ip),y(ip),z(ip),ip=1,np)
  end subroutine vtk_structured_point

  subroutine write_points_vtk_ascii(vunit,np,pos)
    type(vec3d),intent(in) :: pos(:)
    integer,intent(in) :: np,vunit
    integer :: i
    write(vunit,'(a)') '<Points>'
    write(vunit,'(a)') '<DataArray NumberOfComponents="3" format="ascii" type="Float32">'
    do i=1,np
       write(vunit,'(3f10.3)') pos(i)%x,pos(i)%y,pos(i)%z
    end do
    write(vunit,'(a)') '</DataArray>'
    write(vunit,'(a)') '</Points>'
  end subroutine write_points_vtk_ascii
  

  subroutine vtk_structured_header(vunit,vtkf,nx,ny,nz)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: vunit,nx,ny,nz
    character(*),intent(in) :: vtkf
    character(len=maxlen) :: s_buffer
    write(*,*) vtkf
    open(vunit,file=trim(vtkf),form='UNFORMATTED',access='STREAM',action='WRITE',status='REPLACE')
    write(vunit) '# vtk DataFile Version 3.0'//end_rec
    write(vunit) trim(vtkf)//end_rec
    write(vunit) 'BINARY'//end_rec
    write(vunit) 'DATASET STRUCTURED_GRID'//end_rec
    write(s_buffer,'(3i8)') nx,ny,nz
    write(vunit) 'DIMENSIONS '//trim(s_buffer)//end_rec
  end subroutine vtk_structured_header


end module mod_vtk
