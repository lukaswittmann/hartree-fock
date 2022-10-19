program main

use types
use constants
use integrals
use basis

real(dp), dimension(2,3) :: molecule_coords
real(dp), dimension(2) :: z
type(primitive_gaussian), allocatable :: molecule(:,:)

molecule_coords(1,:) = (/0.,0.,0./)  ! Position of H1
molecule_coords(2,:) = (/0.,0.,1.6/) ! Position of H2 in Bohr

molecule = def_molecule(molecule_coords)
z = (/1.0, 1.0/)

call overlap(molecule)
print *
call kinetic_energy(molecule)
print *
call en_interaction(molecule, z)

end program main