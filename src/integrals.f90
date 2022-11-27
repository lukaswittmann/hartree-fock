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

    ! Integral found in TYygve Helgaker, Modern Electronic Structure Theory, C H A P T E R 12,
    ! G A U S S I A N BASIS SETS A N D MOLECULAR I N T E G R A L S
    do i = 1, nbasis
        do j = 1, nbasis
                
        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:)) ! Number of primitives in i..
            do l = 1, size(molecule(j,:))

            ! norm = molecule(i, k)%norm() * molecule(j, l)%norm()
            ! ! Eq. 63
            ! Q_xyz = (molecule(i, k)%coords - molecule(j, l)%coords)
            ! ! Eq. 64, 65
            ! p = (molecule(i, k)%alpha + molecule(j, l)%alpha)
            ! q = (molecule(i, k)%alpha * molecule(j, l)%alpha) / p 
            ! ! Eq. 66
            ! Kab = exp(-q * dot_product(Q_xyz,Q_xyz))
                
            ! coeff = molecule(i, k)%coeff * molecule(j, l)%coeff

            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)

            S(i,j) = S(i,j) + norm * coeff * Kab * (pi / p) ** (1.5)

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

    real(dp) :: norm, p, q, coeff, Kab, S
    real(dp), dimension(3) :: Q_xyz, gP, Pp, PG

    S = 0
    T = 0

    nbasis = size(molecule,1)

    ! Integral found in TYygve Helgaker, Modern Electronic Structure Theory, C H A P T E R 12,
    ! G A U S S I A N BASIS SETS A N D MOLECULAR I N T E G R A L S
    do i = 1, nbasis
        do j = 1, nbasis

        ! Iterate over l and m primitives in basis
        do k = 1, size(molecule(i,:))
            do l = 1, size(molecule(j,:))

            ! coeff = molecule(i, k)%coeff * molecule(j, l)%coeff
            ! norm = molecule(i, k)%norm() * molecule(j, l)%norm()

            ! ! Eq. 63
            ! Q_xyz = (molecule(i, k)%coords - molecule(j, l)%coords)
            ! ! Eq. 64, 65
            ! p = (molecule(i, k)%alpha + molecule(j, l)%alpha)
            ! q = (molecule(i, k)%alpha * molecule(j, l)%alpha) / p 
            ! ! Eq. 66
            ! Kab = exp(-q * dot_product(Q_xyz,Q_xyz))

            call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)

            S = norm * coeff * Kab * (pi / p) ** (1.5)

            gP = (molecule(i, k)%alpha * molecule(i, k)%coords + molecule(j, l)%alpha * molecule(j, l)%coords)
            Pp = gP / p
            PG = Pp - molecule(j, l)%coords

            T(i, j) = T(i, j) + 3 * molecule(j, l)%alpha * S 
            do m = 1, 3
                T(i, j) = T(i, j) - 2 * molecule(j, l)%alpha ** 2 * S * ((PG(m) ** 2) + 0.5 / p)
            end do
                
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
    real(dp) :: boys, norm, p, q, coeff, Kab, S, x, n = 0.0
    real(dp), dimension(3) :: Q_xyz, gP, Pp, PG
    
    V_ne = 0

    natoms = size(molecule, 1) / 2
    nbasis = size(molecule, 1)

    ! Integral found in https://www.mathematica-journal.com/2014/12/08/evaluation-of-gaussian-molecular-integrals-4/
    do atom = 1, natoms
        do i = 1, nbasis
            do j = 1, nbasis

            do k = 1, size(molecule(i,:))
                do l = 1, size(molecule(j,:))

                ! norm = molecule(i, k)%norm() * molecule(j, l)%norm()

                ! ! Eq. 63
                ! Q_xyz = (molecule(i, k)%coords - molecule(j, l)%coords)
                ! ! Eq. 64, 65
                ! p = (molecule(i, k)%alpha + molecule(j, l)%alpha)
                ! q = (molecule(i, k)%alpha * molecule(j, l)%alpha) / p 
                ! ! Eq. 66
                ! Kab = exp(-q * dot_product(Q_xyz,Q_xyz))
                    
                ! coeff = molecule(i, k)%coeff * molecule(j, l)%coeff

                call gauss_product(molecule, i, k, j, l, norm, coeff, p, Kab)

                gP = (molecule(i, k)%alpha * molecule(i, k)%coords + molecule(j, l)%alpha * molecule(j, l)%coords)
                Pp = gP / p
                PG = Pp - molecule_coords(atom,:)
                
                ! ! Calculate Boys Function
                ! x = (p * (dot_product(PG,PG)))
                ! if (x == 0.) then
                !     boys = 1. / (2. * n + 1)
                ! else 
                !     boys = regularized_gamma_p(n + 0.5, x) * gamma(n + 0.5) * (1. / (2.  * x ** (n + 0.5)))
                ! end if

                V_ne(i,j) =  V_ne(i,j) - z(atom) * norm * coeff * Kab * (2.0 * pi / p) * calc_boys((p * (dot_product(PG,PG))), n)
                
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
