module mod_const
  implicit none
  real(8),parameter :: grav = 9.80665d0
  character(1),parameter :: end_rec=char(10)
  character(6),parameter :: output_format = 'ASCII'  ! ASCII or BINARY(NA)
  real(8),parameter :: pi=3.14159265359d0
  integer,parameter :: max_pair = 200
end module mod_const
