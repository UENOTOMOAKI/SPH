module mod_utils
  implicit none
  
  interface str
     module procedure str_i
  endinterface str
  
contains
    function getUnit() result(freeunit)
    implicit none
    integer :: freeunit
    integer :: n
    logical(4) :: used
    n=10
    do
       inquire(unit=n,opened=used)
       if (.NOT. used) then
          freeunit=n
          return
       end if
       n=n+1
    end do
  end function getUnit

  function str_i(n) result(c)
    implicit none
    integer,parameter :: maxlen=100
    integer,intent(in) :: n
    character(len=maxlen):: c
    write(c,'(i8)') n
  end function str_i
  
end module mod_utils
