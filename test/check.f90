program check

use types
use constants
use basis
implicit none

type(primitive_gaussian) :: H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c

type(primitive_gaussian), allocatable :: basisfunctions_H1(:,:),  basisfunctions_H2(:,:), molecule(:,:,:)

H1_1a%alpha = 0.1873113696E+02
H1_1a%coeff =  0.3349460434E-01

H1_1b%alpha = 0.2825394365E+01
H1_1b%coeff =  0.2347269535E+00

H1_1c%alpha = 0.6401216923E+00
H1_1c%coeff =  0.8137573261E+00
        
H1_2a%alpha = 0.1612777588E+00
H1_2a%coeff =  1.0000000000E+00

H1_2b%alpha = 0.0
H1_2b%coeff =  0.0

H1_2c%alpha = 0.0
H1_2c%coeff =  0.0

basisfunctions_H1 = reshape([H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c], (/2,3/))
basisfunctions_H2 = reshape([H1_1a, H1_1b, H1_1c, H1_2a, H1_2b, H1_2c], (/2,3/))

molecule = reshape([basisfunctions_H1, basisfunctions_H2], (/2,2,3/))

print *, molecule

print *, shape(molecule)

end program check
