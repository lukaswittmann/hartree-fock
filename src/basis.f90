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

! for molecule_coordinate in molecule_coordinates:
!     '''6-31G Basis:'''
!     H1_pg1a = primitive_gaussian(0.1873113696E+02, 0.3349460434E-01, molecule_coordinate[0], 0, 0, 0)
!     H1_pg1b = primitive_gaussian(0.2825394365E+01, 0.2347269535E+00, molecule_coordinate[0], 0, 0, 0)
!     H1_pg1c = primitive_gaussian(0.6401216923E+00, 0.8137573261E+00, molecule_coordinate[0], 0, 0, 0)
!     H1_pg2a = primitive_gaussian(0.1612777588E+00, 1.0000000000E+00, molecule_coordinate[0], 0, 0, 0)

!     H2_pg1a = primitive_gaussian(0.1873113696E+02, 0.3349460434E-01, molecule_coordinate[1], 0, 0, 0)
!     H2_pg1b = primitive_gaussian(0.2825394365E+01, 0.2347269535E+00, molecule_coordinate[1], 0, 0, 0)
!     H2_pg1c = primitive_gaussian(0.6401216923E+00, 0.8137573261E+00, molecule_coordinate[1], 0, 0, 0)
!     H2_pg2a = primitive_gaussian(0.1612777588E+00, 1.0000000000E+00, molecule_coordinate[1], 0, 0, 0)

!     H1_1s = [H1_pg1a, H1_pg1b, H1_pg1c]
!     H1_2s = [H1_pg2a]
!     H2_1s = [H2_pg1a, H2_pg1b, H2_pg1c]
!     H2_2s = [H2_pg2a]

!     Z = [1.0, 1.0]
!     atom_coordinates = [np.array(molecule_coordinate[0]), np.array(molecule_coordinate[1])]
!     molecule = [H1_1s, H1_2s, H2_1s, H2_2s]

!     number_occupied_orbitals = 1

!     S = overlap(molecule)
!     T = kinetic(molecule)
!     V_ne= electron_nuclear_interaction(molecule, atom_coordinates, Z)
!     V_ee= electron_electron_interaction(molecule)
!     molecular_terms = [S, T, V_ne, V_ee]

!     V_nn= nuclear_nuclear_interaction(atom_coordinates, Z)

!     electronic_energy = scf_cycle(molecular_terms,  scf_parameters, molecule)

!     total_energy = electronic_energy + V_nn
!     total_energies.append(total_energy)

! end = time.time()
! print("6-31G done in ", round(end-start,1), "s")
! total_energies_631G = total_energies
! total_energies = []

subroutine molecule_6_21G_Basis()
 
   ! tbd

end subroutine molecule_6_21G_Basis

end module basis
