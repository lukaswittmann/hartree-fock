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

contains

! Normalization of primitive_gaussian
real(dp) function norm(self)
    class(primitive_gaussian), intent(in) :: self
    norm = ((2.0 * self%alpha) / pi) ** (3.0/4.0)
end function norm

! 6-31G Basis
function def_molecule(molecule_coords)
    real(dp), dimension(2,3), intent(in) :: molecule_coords
    
    type(primitive_gaussian) :: H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c ! Basis funcions H1, built from primitive gaussians
    type(primitive_gaussian) :: H2_1a, H2_1b, H2_1c, H2_2a, H2_2b, H2_2c ! Basis funcions H2 (same as H1)

    type(primitive_gaussian), allocatable :: H1(:,:),  H2(:,:), def_molecule2(:,:)
    type(primitive_gaussian), dimension(4,3) :: def_molecule

    integer :: i, j, k

    H1_1a%alpha = 0.1873113696E+02
    H1_1a%coeff =  0.3349460434E-01
    H1_1a%coords = molecule_coords(1,:)

    H1_1b%alpha = 0.2825394365E+01
    H1_1b%coeff =  0.2347269535E+00
    H1_1b%coords = molecule_coords(1,:)

    H1_1c%alpha = 0.6401216923E+00
    H1_1c%coeff =  0.8137573261E+00
    H1_1c%coords = molecule_coords(1,:)
            
    H1_2a%alpha = 0.1612777588E+00
    H1_2a%coeff =  1.0000000000E+00
    H1_2a%coords = molecule_coords(1,:)

    ! Coeff = 0 primitive gaussians, quick and dirty way of implementation, maybe someday
    ! I will explicitly define wrapper types for each sub-array type    
    H1_2b%alpha = 1.0
    H1_2b%coeff =  0.0
    H1_2b%coords = molecule_coords(1,:)

    H1_2c%alpha = 1.0
    H1_2c%coeff =  0.0
    H1_2c%coords = molecule_coords(1,:)

        
    H2_1a%alpha = 0.1873113696E+02
    H2_1a%coeff =  0.3349460434E-01
    H2_1a%coords = molecule_coords(2,:)

    H2_1b%alpha = 0.2825394365E+01
    H2_1b%coeff =  0.2347269535E+00
    H2_1b%coords = molecule_coords(2,:)

    H2_1c%alpha = 0.6401216923E+00
    H2_1c%coeff =  0.8137573261E+00
    H2_1c%coords = molecule_coords(2,:)
            
    H2_2a%alpha = 0.1612777588E+00
    H2_2a%coeff =  1.0000000000E+00
    H2_2a%coords = molecule_coords(2,:)

    ! Coeff = 0 primitive gaussians, quick and dirty way of implementation, maybe someday
    ! I will explicitly define wrapper types for each sub-array type
    H2_2b%alpha = 1.0
    H2_2b%coeff =  0.0
    H2_2b%coords = molecule_coords(2,:)

    H2_2c%alpha = 1.0
    H2_2c%coeff =  0.0
    H2_2c%coords = molecule_coords(2,:)
        
    H1 = reshape([H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c], (/3,2/)) ! H1_1s and H1_2s
    H2 = reshape([H2_1a, H2_1b, H2_1c, H2_2a, H2_2b, H2_2c], (/3,2/)) ! H2_1s and H2_2s

    def_molecule2 = reshape([H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c,H2_1a, H2_1b, H2_1c, H2_2a, H2_2b, H2_2c], (/3,4/))
    
    ! Resort array to make it more intuitive
    do j = 1,3
        do i = 1,4
            def_molecule(i,j) = def_molecule2(j,i)
        end do
    end do

end function def_molecule

end module basis
