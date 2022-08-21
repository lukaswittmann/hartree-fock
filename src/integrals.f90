module integrals

use types
use constants
use basis

contains

subroutine overlap(molecule)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer :: nbasis, natoms, nprimitives_i, nprimitives_j, i, j, k, l
    real(dp) :: norm, p, q, coeff, Kab
    real(dp), dimension(3) :: Q_xyz
    real(dp), dimension(INT(size(molecule,1)),INT(size(molecule,1))) :: S

    nbasis = size(molecule,1)

    ! Integral found in TYygve Helgaker, Modern Electronic Structure Theory, C H A P T E R 12,
    ! G A U S S I A N BASIS SETS A N D MOLECULAR I N T E G R A L S
    do i = 1, nbasis
        do j = 1, nbasis
                
        nprimitives_i = size(molecule(i,:))
        nprimitives_j = size(molecule(j,:))

        ! Iterate over l and m primitives in basis
        do k = 1, nprimitives_i
            do l = 1, nprimitives_j

            norm = molecule(i, k)%norm() * molecule(j, l)%norm()
            ! Eq. 63
            Q_xyz = (molecule(i, k)%coords - molecule(j, l)%coords)
            ! Eq. 64, 65
            p = (molecule(i, k)%alpha + molecule(j, l)%alpha)
            q = (molecule(i, k)%alpha * molecule(j, l)%alpha) / p 
            ! Eq. 66
            Kab = exp(-q * dot_product(Q_xyz,Q_xyz))
                
            coeff = molecule(i, k)%coeff * molecule(j, l)%coeff

            S(i,j) = S(i,j) + norm * coeff * Kab * (pi / p) ** (1.5)

            end do
        end do   
        end do
    end do

    print *, S
    
end subroutine overlap

end module integrals
