module st_params
  implicit none
  ! Common parameters for Suzuki–Trotter tests
  integer, parameter :: n = 200 !ngrid pt
  real(8), parameter :: beta  = 1.0d0
  real(8), parameter :: a     = -11.0d0
  real(8), parameter :: b     =  11.0d0
  integer, parameter :: nwpot = 1
end module st_params
