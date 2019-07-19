module mod_dist

implicit none

! C Interface
interface

subroutine dist_getrefpara(energy0,mass0,a0,z0) bind(C, name="getrefpara")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE), intent(inout) :: energy0, mass0
  integer(kind=C_INT), intent(inout) :: a0, z0
end subroutine dist_getrefpara

subroutine dist_writefile(fileName, strlen) bind(C, name="writefile_f")
  use, intrinsic :: iso_c_binding
  character(kind=C_CHAR,len=1), intent(in) :: fileName
  integer(kind=C_INT), value, intent(in) :: strlen
end subroutine dist_writefile

subroutine dist_readfile(fileName, strlen) bind(C, name="readfile_f")
  use, intrinsic :: iso_c_binding
  character(kind=C_CHAR,len=1), intent(in) :: fileName
  integer(kind=C_INT), value, intent(in) :: strlen
end subroutine dist_readfile

subroutine dist_initializedistribution(numdist) bind(C, name="initializedistribution")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: numdist
end subroutine dist_initializedistribution

subroutine dist_setphysicalcut(variable, min, max) bind(C, name="setphysicalcut")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable
  real(kind=C_DOUBLE), value, intent(in) :: min, max
end subroutine dist_setphysicalcut

subroutine dist_setnormalizedcut(variable, min, max) bind(C, name="setnormalizedcut")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable
  real(kind=C_DOUBLE), value, intent(in) :: min, max
end subroutine dist_setnormalizedcut

subroutine dist_setdistribution(distcount) bind(C, name="setdistribution")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: distcount
end subroutine dist_setdistribution

subroutine dist_addclosedorbit(orbit) bind(C, name="addclosedorbit")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(6), intent(in) :: orbit
end subroutine dist_addclosedorbit

subroutine dist_get6trackcoord(x,xp,y,yp,sigma,deltap,npart) bind(C, name="get6trackcoord")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(*), intent(inout) :: x,xp,y,yp,sigma,deltap
  integer(kind=C_INT), intent(inout) :: npart
end subroutine dist_get6trackcoord

subroutine dist_sete0andmass0(energy0, mass0) bind(C, name="sete0andmass0")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),value, intent(in) :: energy0, mass0
end subroutine dist_sete0andmass0

subroutine dist_setemitt12(emitt1, emitt2) bind(C, name="setemitt12")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),value, intent(in) :: emitt1, emitt2
end subroutine dist_setemitt12

subroutine dist_setemitt3(emitt3) bind(C, name="setemitt3")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),value, intent(in) :: emitt3
end subroutine dist_setemitt3

subroutine dist_settasmatrix(tas) bind(C, name="settasmatrix")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(36), intent(in) :: tas
end subroutine dist_settasmatrix

subroutine dist_setcooords(xn,xnp,yn,ynp,zn,znp,npart, coordtype) bind(C, name="setcoords")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(*), intent(in) :: xn,xnp,yn,ynp,zn,znp
  integer(kind=C_INT),value,  intent(in) :: npart, coordtype
end subroutine dist_setcooords

subroutine dist_getarraylength(totlength) bind(C, name="getarraylength")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), intent(out) :: totlength
end subroutine dist_getarraylength

subroutine dist_settotalsteps(totlength) bind(C, name="settotalsteps")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: totlength
end subroutine dist_settotalsteps

subroutine dist_setscan_para_diagonal(variable, variable_type, space_type, start_value, stop_value) &
 bind(C, name="setscan_para_diagonal")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable, variable_type, space_type
  real(kind=C_DOUBLE), value, intent(in) :: start_value, stop_value
end subroutine dist_setscan_para_diagonal

subroutine dist_setscan_para_grid(variable, variable_type, space_type, start_value, stop_value, length) &
bind(C, name="setscan_para_diagonal")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable, variable_type, space_type, length
  real(kind=C_DOUBLE), value, intent(in) :: start_value, stop_value
end subroutine dist_setscan_para_grid


end interface 
end module