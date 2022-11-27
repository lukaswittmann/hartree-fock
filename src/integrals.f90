module integrals
use types
use constants
use basis
use functions

use stdlib_specialfunctions_gamma!, only: lig => lower_incomplete_gamma


contains


subroutine overlap(molecule)
    implicit none
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l
    real(dp) :: norm, p, q, coeff, Kab
    real(dp), dimension(3) :: Q_xyz
    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: S

    S = 0

    nbasis = size(molecule,1)

    ! Helgaker, Modern Electronic Structure Theory
    do i = 1, nbasis
        do j = 1, nbasis
                
        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:)) ! Number of primitives in i..
            do l = 1, size(molecule(j,:))

            ! SO A.9 with p exponent of new gaussian, norm and coeff needed due to multiple gaussians per bf
            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)
            S(i,j) = S(i,j) + (pi / p) ** (1.5) * Kab * norm * coeff

            end do
        end do   
        end do
    end do

    print *, S
    
end subroutine overlap


subroutine kinetic_energy(molecule)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l, m
    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: T

    real(dp) :: norm, p, coeff, Kab, S, ab, tmp, tmp2
    real(dp), dimension(3) :: Q_xyz, gP, Pp, PG

    S = 0
    T = 0

    nbasis = size(molecule,1)

    do i = 1, nbasis
        do j = 1, nbasis

        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:))
            do l = 1, size(molecule(j,:))

            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)
            S = (pi / p) ** (1.5) * Kab * norm * coeff

            ! OZ A.11
            ab = molecule(i, k)%alpha * molecule(j, l)%alpha
            T(i, j) = T(i, j) + (ab/p) * (3 + (2 * log(Kab))) * s

            end do
        end do   
        end do
    end do

    print *, T

end subroutine kinetic_energy


subroutine en_interaction(molecule, molecule_coords, z)
    implicit none
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    real(dp), dimension(2), intent(in) :: z
    real(dp), dimension(2,3), intent(in) :: molecule_coords
    integer :: atom, natoms, nbasis, i, j, k, l

    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: V_ne
    real(dp) :: boys, norm, p, coeff, Kab, S, x, n = 0.0
    real(dp), dimension(3) :: Rp, Rnp
    
    V_ne = 0

    natoms = size(molecule, 1) / 2
    nbasis = size(molecule, 1)

    ! Integral found in https://www.mathematica-journal.com/2014/12/08/evaluation-of-gaussian-molecular-integrals-4/
    do atom = 1, natoms
        do i = 1, nbasis
            do j = 1, nbasis

            do k = 1, size(molecule(i,:))
                do l = 1, size(molecule(j,:))

                call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)

                ! Center of gaussian product (als in gauss_product, will be shortened later)
                Rp = (molecule(i, k)%alpha * molecule(i, k)%coords + molecule(j, l)%alpha * molecule(j, l)%coords) / p
                ! Vector of R nuclear to R gaussian product
                Rnp = Rp - molecule_coords(atom,:) !! Rp - Rc

                ! OZ A.33 but with boys instead of erf
                V_ne(i,j) =  V_ne(i,j) - (norm * coeff) * (2.0 * pi * z(atom) * Kab / p  * calc_boys(p * (dot_product(Rnp,Rnp)), n))
                
                end do
            end do
            end do
        end do
    end do

    print *, V_ne

end subroutine en_interaction


subroutine ee_interaction(molecule)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l, pi, pj, pk, pl

    real(dp) :: norm, coeff, pij, pkl
    real(dp), dimension(3) :: gPij, gPkl

    nbasis = size(molecule, 1)

    ! Sum over Basisfunctions
    do i = 1, nbasis
        do j = 1, nbasis
            do k = 1, nbasis
                do l = 1, nbasis

                ! Sum over Primitives in each Basisfunction
                do pi = 1, size(molecule(i,:))
                    do pj = 1, size(molecule(j,:))
                        do pk = 1, size(molecule(k,:))
                            do pl = 1, size(molecule(l,:))

                            ! No intent for better workflow
                            norm = molecule(i, pi)%norm() * molecule(j, pj)%norm() * molecule(k, pk)%norm() * molecule(l, pl)%norm()
                            coeff = molecule(i, pi)%coeff * molecule(j, pj)%coeff * molecule(k, pk)%coeff * molecule(l, pl)%coeff

                            ! Use gaussian product theorem to make one center integral out of two center integral
                            pij = molecule(i, pi)%alpha + molecule(j, pj)%alpha
                            pkl = molecule(k, pk)%alpha + molecule(l, pl)%alpha

                            gPij = molecule(i, pi)%alpha * molecule(i, pi)%coords + molecule(j, pj)%alpha * molecule(j, pj)%coords
                            gPkl = molecule(k, pk)%alpha * molecule(k, pk)%coords + molecule(l, pl)%alpha * molecule(l, pl)%coords


                            end do
                        end do            
                    end do
                end do

                end do
            end do
        end do
    end do

end subroutine ee_interaction





end module integrals
