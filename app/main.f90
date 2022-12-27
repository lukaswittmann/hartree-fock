program main

use types
use constants
use integrals
use basis
use cycles

real(dp), dimension(2,3) :: molecule_coords
real(dp), dimension(2) :: z
type(primitive_gaussian), allocatable :: molecule(:,:)

integer :: i

! Def needed integrals
real(dp), dimension(4,4) :: S, T, V_ne
real(dp), dimension(4,4,4,4) :: V_ee
real(dp) :: V_nn

! SCF parameters
real(dp) :: conv_crit = 1.0e-6
integer :: max_iter = 100

molecule_coords(1,:) = (/0.,0.,0./)  ! Position of H1
molecule_coords(2,:) = (/0.,0.,1.6/) ! Position of H2 in Bohr

molecule = def_molecule(molecule_coords)
z = (/1.0, 1.0/)

call overlap(molecule, S)
call kinetic_energy(molecule, T)
call en_interaction(molecule, molecule_coords, z, V_ne)
call ee_interaction(molecule, V_ee)
call nn_interaction(molecule_coords, z, V_nn)

call scf_cycle(S, T, V_ne, V_ee, conv_crit, max_iter, molecule)



end program main