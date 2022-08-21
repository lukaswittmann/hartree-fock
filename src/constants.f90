module constants

use types, only: dp
implicit none

private
public pi, e_, i_, kB, R_, NA

real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_   = 2.7182818284590452353602874713527_dp
real(dp), parameter :: kB   = 1.380649E-23_dp
real(dp), parameter :: NA   = 6.02214076E-23_dp
real(dp), parameter :: R_   = 8.314472_dp
complex(dp), parameter :: i_ = (0, 1)

end module
