module functions
use types
use constants
use basis    


    
contains
    
subroutine gauss_product(molecule, a, b, c, d, norm, coeff, p, Kab) !result(output)
    type(primitive_gaussian), intent(in) :: molecule(:,:)
    integer, intent(in) :: a, b, c, d
    real(dp) :: diff
    real(dp) :: p, norm, Kab, coeff
    ! real(dp), dimension(3) :: Rp

    ! Szabo Ostlund p. 311
    coeff = molecule(a, b)%coeff * molecule(c, d)%coeff
    ! Product exponent ! ! Eq. 64, 65
    p = molecule(a, b)%alpha + molecule(c, d)%alpha
    ! Normalization
    norm = ((4 * molecule(a, b)%alpha * molecule(c, d)%alpha) / (pi ** 2)) ** (3.0 / 4.0)
    ! Product prefactor ! ! Eq. 63, 66
    Kab = exp(- molecule(a, b)%alpha * molecule(c, d)%alpha/ p &
    * dot_product(molecule(a, b)%coords - molecule(c, d)%coords,molecule(a, b)%coords - molecule(c, d)%coords))
    ! Product center
    ! Rp = (molecule(a, b)%alpha * molecule(a, b)%coords + molecule(c, d)%alpha * molecule(c, d)%coords) / p

end subroutine gauss_product


function calc_boys(x, n) result(res)
    use stdlib_specialfunctions_gamma!, only: lig => lower_incomplete_gamma

    ! SZ p. 412, electrostatics of spherical Gaussian charge distributions are determined by the Boys function
    real(dp), intent(in) :: x, n
    real(dp) :: res

    if (x == 0.) then
        res = 1. / (2. * n + 1)
    else 
        res = regularized_gamma_p(n + 0.5, x) * gamma(n + 0.5) * (1. / (2.  * x ** (n + 0.5)))
    end if

end function calc_boys

end module functions