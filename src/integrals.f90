module integrals

use types
use constants
use basis

contains

subroutine overlap(molecule)
    type(molecule_structure), intent(in) :: molecule

    print *, molecule

end subroutine overlap

end module integrals
