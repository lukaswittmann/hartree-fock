module basis

use types
use constants

! Define primitive gaussian for basis functions
type primitive_gaussian
    real(dp) :: alpha, coeff
    real(dp), dimension(3) :: coords
    contains
        procedure :: norm
end type primitive_gaussian

type orbital_1pg
    type(primitive_gaussian) :: a
end type orbital_1pg

type orbital_3pg
    type(primitive_gaussian) :: a, b, c
end type orbital_3pg

type molecule_structure
    type(orbital_3pg) :: H1_1s
    type(orbital_1pg) :: H1_2s
    type(orbital_3pg) :: H2_1s
    type(orbital_1pg) :: H2_2s
end type molecule_structure

contains

! Normalization of primitive_gaussian
real(dp) function norm(self)
    class(primitive_gaussian), intent(in) :: self
    norm = ((2.0 * self%alpha) / pi ** (0.75))
end function norm

subroutine molecule_6_21G_Basis()
 
   ! tbd

end subroutine molecule_6_21G_Basis

end module basis
