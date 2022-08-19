program main

use types
use constants
use integrals
use basis


real(dp), dimension(2,3) :: molecule_coords
    
type(primitive_gaussian) :: H1_1a, H1_1b, H1_1c, H1_2a ! Basis funcions H1, built from primitive gaussians
type(primitive_gaussian) :: H2_1a, H2_1b, H2_1c, H2_2a ! Basis funcions H2 (same as H1)

type(orbital_3pg) :: H1_1s, H2_1s
type(orbital_1pg) :: H1_2s, H2_2s

type(molecule_structure) :: molecule

molecule_coords(1,:) = (/0.,0.,0./)  ! Position of H1
molecule_coords(2,:) = (/0.,0.,1.6/) ! Position of H2 in Bohr

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

H1_1s%a = H1_1a
H1_1s%b = H1_1b
H1_1s%c = H1_1c
H1_2s%a = H1_2a
    
molecule%H1_1s = H1_1s
molecule%H1_2s = H1_2s
    
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
    
H2_1s%a = H2_1a
H2_1s%b = H2_1b
H2_1s%c = H2_1c
H2_2s%a = H2_2a

molecule%H2_1s = H2_1s
molecule%H2_2s = H2_2s

call overlap(molecule)

end program main