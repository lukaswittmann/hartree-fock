module integrals
use types
use constants
use basis
use functions

use stdlib_specialfunctions_gamma!, only: lig => lower_incomplete_gamma


contains


subroutine overlap(molecule, S)
    implicit none
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l
    real(dp) :: norm, p, coeff, Kab
    real(dp), dimension(3) :: Q_xyz
    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: S
    real(dp), dimension(3) :: Rp

    S = 0

    nbasis = size(molecule,1)

    ! Helgaker, Modern Electronic Structure Theory
    do i = 1, nbasis
        do j = 1, nbasis
                
        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:)) ! Number of primitives in i..
            do l = 1, size(molecule(j,:))

            ! SO A.9 with p exponent of new gaussian, norm and coeff needed due to multiple gaussians per bf
            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab, Rp)
            S(i,j) = S(i,j) + (pi / p) ** (1.5) * Kab * norm * coeff

            end do
        end do   
        end do
    end do

    !print *, S
    
end subroutine overlap


subroutine kinetic_energy(molecule, T)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l, m
    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: T

    real(dp) :: norm, p, coeff, Kab, S, ab
    real(dp), dimension(3) ::  Rp

    S = 0
    T = 0

    nbasis = size(molecule,1)

    do i = 1, nbasis
        do j = 1, nbasis

        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:))
            do l = 1, size(molecule(j,:))

            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab, Rp)
            S = (pi / p) ** (1.5) * Kab * norm * coeff

            ! OZ A.11
            ab = molecule(i, k)%alpha * molecule(j, l)%alpha
            T(i, j) = T(i, j) + (ab/p) * (3 + (2 * log(Kab))) * s

            end do
        end do   
        end do
    end do

    !print *, T

end subroutine kinetic_energy


subroutine en_interaction(molecule, molecule_coords, z, V_ne)
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

                call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab, Rp)

                ! Vector of R nuclear to R gaussian product
                Rnp = Rp - molecule_coords(atom,:) !! Rp - Rc

                ! OZ A.33 but with boys instead of erf
                V_ne(i,j) =  V_ne(i,j) - (norm * coeff) * (2.0 * pi * z(atom) * Kab / p  * calc_boys(p * (dot_product(Rnp,Rnp)), n))
                
                end do
            end do
            end do
        end do
    end do

    !print *, V_ne

end subroutine en_interaction


subroutine ee_interaction(molecule, V_ee)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, i, j, k, l, pi, pj, pk, pl

    real(dp) :: normGlob, coeffGlob ! global
    real(dp) :: normGp1, coeffGp1, pGp1, KabGp1, normGp2, coeffGp2, pGp2, KabGp2 ! gaussian product 1 & 2
    real(dp), dimension(3) :: RpGp1, RpGp2, Rpp

    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1)),INT(size(molecule,1)),INT(size(molecule,1))) :: V_ee

    real(dp) :: fact = 1.0, n = 0.

    integer :: x, y

    nbasis = size(molecule, 1)

    V_ee(:,:,:,:) = 0

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

                            !normGlob = molecule(i, pi)%norm() * molecule(j, pj)%norm() * molecule(k, pk)%norm() * molecule(l, pl)%norm()
                            !coeffGlob = molecule(i, pi)%coeff * molecule(j, pj)%coeff * molecule(k, pk)%coeff * molecule(l, pl)%coeff

                            ! ! Use gaussian product theorem to make one center integral out of two center integral
                            ! pij = molecule(i, pi)%alpha + molecule(j, pj)%alpha
                            ! pkl = molecule(k, pk)%alpha + molecule(l, pl)%alpha

                            ! gPij = molecule(i, pi)%alpha * molecule(i, pi)%coords + molecule(j, pj)%alpha * molecule(j, pj)%coords
                            ! gPkl = molecule(k, pk)%alpha * molecule(k, pk)%coords + molecule(l, pl)%alpha * molecule(l, pl)%coords


call gauss_product(molecule, i, pi, j, pj, normGp1, coeffGp1, pGp1, KabGp1, RpGp1)
call gauss_product(molecule, k, pk, l, pl, normGp2, coeffGp2, pGp2, KabGp2, RpGp2)

normGlob = normGp1 * normGp2
coeffGlob = coeffGp1 * coeffGp2

! Due to empty basisfunctions
if (coeffGlob /= 0.0) then


Rpp = RpGp1 - RpGp2 ! Rpp = Rij - Rkl

!!! for i = 2, 4 its wrong !!!
V_ee(i,j,k,l) = V_ee(i,j,k,l) + normGlob * coeffGlob &
                * (2. * pi ** (5./2.)) / (pGp1 * pGp2 * sqrt(pGp1 + pGp2)) & ! term1
                * KabGp1 * KabGp2 & !term 2
                * calc_boys(((pGp1 * pGp2 / (pGp1 + pGp2)) * (dot_product(Rpp,Rpp))), n)

end if


                            end do
                        end do            
                    end do
                end do

                end do
            end do
        end do
    end do

    !print *, V_ee

end subroutine ee_interaction


subroutine nn_interaction(molecule_coords, z, V_nn)
    real(dp), dimension(2,3), intent(in) :: molecule_coords
    real(dp), dimension(2), intent(in) :: z

    integer :: natoms, i, j
    real(dp) :: V_nn, Rpx, Rpy, Rpz, absRp
    real(dp), dimension(3) :: Rp

    V_nn = 0
    natoms = size(molecule_coords, 1)

    do i = 1, natoms
        do j = 1, natoms
            if (j > i) then

            ! Vector of R nuclear to R gaussian product
            Rpx = molecule_coords(i,1) - molecule_coords(j,1)
            Rpy = molecule_coords(i,2) - molecule_coords(j,2)
            Rpz = molecule_coords(i,3) - molecule_coords(j,3)
                
            absRp = sqrt(Rpx**2 + Rpy**2 + Rpz**2)
            V_nn = V_nn + z(i) * z(j) / absRp

            end if
        end do
    end do

    !print *, V_nn

end subroutine nn_interaction


end module integrals
