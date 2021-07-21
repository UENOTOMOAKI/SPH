module mod_type
  implicit none
  type vec3d
     real(8) x,y,z
  end type vec3d

  type tensor
     real(8) xx,yy,zz,xy,yx,yz,zy,zx,xz
  end type tensor
end module mod_type
