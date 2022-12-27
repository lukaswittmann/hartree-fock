module cycles
    use types
    use constants
    use basis
    use functions
    
    use stdlib_specialfunctions_gamma!, only: lig => lower_incomplete_gamma
    
    
    contains
    
    
    
    subroutine compute_G(density_matrix, V_ee, G)
        real(dp), dimension(:,:), intent(in) :: density_matrix
        real(dp), dimension(:,:,:,:), intent(in) :: V_ee
    
        real(dp), dimension(size(density_matrix, 1),size(density_matrix, 1)) :: G
        integer :: nbasis, i, j, k, l
    
        nbasis = size(density_matrix, 1)
    
        G = 0
    
        do i = 1, nbasis
            do j = 1, nbasis
                do k = 1, nbasis
                    do l = 1, nbasis
                        
                    ! V_ee(i,j,k,l) is Hartree Term, V_ee(i,l,k,j) is Fock-Slater Term
                    G(i,j) = G(i,j) + density_matrix(k,l) * (V_ee(i,j,k,l) - 0.5_dp * V_ee(i,l,k,j))
    
                    end do
                end do
            end do
        end do
    
        !print *, G
    
    end subroutine compute_G
    
    
    subroutine compute_density_matrix(molecular_orbitals, density_matrix)
        real(dp), dimension(:,:), intent(in) :: molecular_orbitals
        real(dp), dimension(size(molecular_orbitals, 1),size(molecular_orbitals, 1)), intent(out) :: density_matrix
    
        integer :: nbasis, i, j, k, l, occ = 2 ! num of occupied orbitals
        
        nbasis = size(molecular_orbitals, 1)
    
        do i = 1, nbasis
            do j = 1, nbasis
                do k = 1, occ
                    density_matrix(i,j) = density_matrix(i,j) + molecular_orbitals(i,k) * molecular_orbitals(j,k)
                end do
            end do
        end do
    
        !print *, density_matrix
    
    end subroutine compute_density_matrix
    
    
    subroutine compute_electronic_energy(density_matrix, T, V_ne, G, electronic_energy)
        real(dp), dimension(:,:), intent(in) :: density_matrix
        real(dp), dimension(:,:), intent(in) :: T
        real(dp), dimension(:,:), intent(in) :: V_ne
        real(dp), dimension(:,:), intent(in) :: G
        real(dp) :: electronic_energy
        real(dp), dimension(size(density_matrix, 1),size(density_matrix, 1)) :: H_core
    
        integer :: nbasis, i, j
    
        H_core = V_ne + T
        nbasis = size(density_matrix, 1)
        electronic_energy = 0
    
        do i = 1, nbasis
            do j = 1, nbasis
                !H_core = H_core + T(i,j) + V_ne(i,j)
    
                electronic_energy = electronic_energy + 0.5_dp * density_matrix(i,j) * (H_core(i,j) + G(i,j))
            end do
        end do
    
        !print *, electronic_energy
    
    end subroutine compute_electronic_energy
    
    
    subroutine scf_cycle(S, T, V_ne, V_ee, conv_crit, max_iter, molecule)
        real(dp), dimension(:,:), intent(in) :: S
        real(dp), dimension(:,:), intent(in) :: T
        real(dp), dimension(:,:), intent(in) :: V_ne
        real(dp), dimension(:,:,:,:), intent(in) :: V_ee
        real(dp), intent(in) :: conv_crit
        integer, intent(in) :: max_iter
        type(primitive_gaussian), intent(in) :: molecule(:,:)
        
        real(dp), dimension(size(S, 1),size(S, 1)) :: S_inv
        real(dp), dimension(size(S, 1),size(S, 1)) :: F
        real(dp), dimension(size(S, 1),size(S, 1)) :: F_unit_S
        real(dp), dimension(size(S, 1),size(S, 1)) :: G
        real(dp), dimension(size(S, 1),size(S, 1)) :: density_matrix
        real(dp), dimension(size(S, 1),size(S, 1)) :: molecular_orbitals
        real(dp), dimension(size(S, 1),size(S, 1)) :: molecular_orbitals_old
        real(dp), dimension(size(S, 1),size(S, 1)) :: C
        real(dp), dimension(size(S, 1),size(S, 1)) :: C_old
        real(dp), dimension(size(S, 1),size(S, 1)) :: F_prime
        real(dp), dimension(size(S, 1),size(S, 1)) :: F_prime_old
    
        real(dp) :: electronic_energy, electronic_energy_old
        integer :: ipiv, info
        
        nbasis = size(S, 1)
    
        do scf_inter = 1, max_iter
    
            electronic_energy_old = electronic_energy
    
            call compute_G(density_matrix, V_ee, G)
            F = T + V_ne + G
               
    
    
    
            
    
        end do
    
    end subroutine scf_cycle
    
    
    end module cycles
    