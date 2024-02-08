!> This is your module to write your very own SCF program.
module scf_main
   !> Include standard Fortran environment for IO
   use iso_fortran_env, only : output_unit, error_unit

   !> interface to LAPACK's double precision symmetric eigenvalue solver (dspev)
   use linear_algebra, only : solve_spev

   !> expansion of slater-functions into contracted gaussians,
   use slater, only : expand_slater

   !> calculates one-electron integrals and two-electron integrals over
   use integrals, only : oneint, twoint

   !> prints a matrix quantity to screen
   use print_matrix, only : write_vector, write_matrix

   !> other tools that may help you jump ahead with I/O-heavy tasks
   use io_tools, only : read_line

   !> Personal tools
   use tools

   !> Always declare everything explicitly
   implicit none

   !> All subroutines within this module are not exported, except for scf_prog
   private
   public :: scf_prog

contains

   !> This is the entry point to your program, do not modify the dummy arguments
   subroutine scf_prog(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
                      ehf, MAX_SCF, TOL_SCF, do_mp2, print_level)
                      
      implicit none

      !> System specific data
      !> Number of atoms
      integer, intent(in) :: nat
      !> Number of electrons
      integer, intent(in) :: nel
      !> Number of basis functions
      integer, intent(in) :: nbf
      !> number of primitive Gaussians in slater expansion (STO-nG)
      integer, intent(in) :: ng

      !> Atom coordinates of the system, all distances in bohr
      real(wp), intent(in), allocatable :: xyz(:, :)

      !> Nuclear charges
      real(wp), intent(in), allocatable :: chrg(:)

      !> Slater expnts of basis functions
      real(wp), intent(in), allocatable :: zeta(:)

      !> maps basis function to atom
      !> e.g.: (1, 1, 2) -> 1st bf from A1, 2nd bf from A1, 3rd bf from A3
      integer, intent(in), allocatable :: bf_atom_map(:)

      !> convergence criteria
      integer :: iter
      real(wp), intent(in) :: TOL_SCF
      integer, intent(in) :: MAX_SCF

      !> amount of info printed during program run (0, 1, 2)
      integer, intent(in), optional :: print_level

      !> do MP2 energy calculation at end
      integer, intent(in), optional :: do_mp2

      !> HF energy
      real(wp), intent(out) :: ehf
      real(wp) :: emp2


      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> energies
      !> nuclear repulsion energy
      real(wp) :: enn
      real(wp) :: escf

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> expnts of primitive Gaussian functions
      real(wp), dimension(:, :), allocatable :: expnts

      !> coeffs of primitive Gaussian functions
      real(wp), dimension(:, :), allocatable :: coeffs

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(wp), dimension(:, :), allocatable :: S, T, V
      real(wp), dimension(:), allocatable :: S_packed, T_packed, V_packed

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(wp), dimension(:, :), allocatable :: X        ! orthonormalizer
      real(wp), dimension(:, :), allocatable :: P        ! density matrix
      real(wp), dimension(:, :), allocatable :: C        ! coeffs matrix
      real(wp), dimension(:), allocatable :: C_packed    ! coeffs matrix
      real(wp), dimension(:, :), allocatable :: C_prime        ! coeffs matrix
      real(wp), dimension(:), allocatable :: C_prime_packed    ! coeffs matrix
      real(wp), allocatable, dimension(:) :: eps         ! eigenvalues
      integer, dimension(:, :), allocatable :: n_occ     ! occupation matrix
      real(wp), dimension(:, :), allocatable :: F        ! Fock matrix
      real(wp), dimension(:, :), allocatable :: F_prime  ! Fock matrix
      real(wp), dimension(:), allocatable :: F_packed    ! Fock matrix
      real(wp), dimension(:), allocatable :: F_prime_packed  ! Fock matrix

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: alloc_stat
      real(wp), dimension(:, :, :, :), allocatable :: two_ints ! 4D array of tei

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: i, j, k, l

      101 format(//, "---------------------------------------------------------------")
      102 format("---------------------------------------------------------------")

      !> Print xyz
      if (present(print_level)) then
         call write_matrix(xyz, name='xyz')
      end if

      !*********************************************************
      !******************* BASIS SET SETUP  ********************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Basis set setup"
      write (*, 102)

      !> Basis coefficients
      allocate (coeffs(ng, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of coeffs failed."

      !> Basis exponents
      allocate (expnts(ng, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of expnts failed."

      !> Expand slater functions into primitive gaussians
      do i = 1, nbf
         call expand_slater(ng, zeta(i), expnts(:, i), coeffs(:, i))
      end do

      !> Print the expanded slater functions
      if (present(print_level)) then
         call write_matrix(coeffs, name='coeffs')
         call write_matrix(expnts, name='expnts')
      end if

      !*********************************************************
      !******************* NUCLEAR REPULSION *******************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Nuclear repulsion"
      write (*, 102)

      !> Nuclear repulsion energy
      enn = 0.0_wp
      do i = 1, nat
         do j = i+1, nat
            enn = enn + chrg(i) * chrg(j) / sqrt(sum((xyz(:, i) - xyz(:, j))**2))
         end do
      end do

      !> Print the nuclear repulsion energy
      if (present(print_level)) then
         write(*, '(A, F20.14)') 'enn:', enn
      end if

      !*********************************************************
      !******************* ONE-ELECTRON INTEGRALS **************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "One-electron integrals"
      write (*, 102)

      !> Overlap matrix
      allocate(S(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of S failed."
      allocate(S_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of S_packed failed."

      !> Kinetic energy matrix
      allocate(T(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of T failed."
      allocate(T_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of T_packed failed."

      !> Nuclear attraction matrix
      allocate(V(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of V failed."
      allocate(V_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of V_packed failed."

      !> Iterate over all basis functions
      do i = 1, nbf
         do j = 1, nbf
            call oneint( & 
               xyz, chrg, &
               xyz(:,bf_atom_map(i)), xyz(:,bf_atom_map(j)), &
               expnts(:, i), expnts(:, j), &
               coeffs(:, i), coeffs(:, j), &
               S(i, j), T(i, j), V(i, j))
         end do
      end do

      !> Print the one-electron integrals
      if (present(print_level)) then
         call write_matrix(S, name='S')
         call write_matrix(T, name='T')
         call write_matrix(V, name='V')
      end if

      !> Pack the matrices
      call pack_matrix(S, S_packed)
      call pack_matrix(T, T_packed)
      call pack_matrix(V, V_packed)

      !*********************************************************
      !******************* SYMMETRIC ORTHOGONALIZER ************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Symmetric orthonormalizer"
      write (*, 102)

      !> Symmetric orthonormalizer
      allocate(X(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of X failed."

      !> Calculate the symmetric orthonormalizer
      call calc_symmetric_orthonormalizer(S_packed, nbf, X)

      !> Print the symmetric orthonormalizer
      if (present(print_level)) then
         call write_matrix(X, name='X')
      end if

      !*********************************************************
      !******************* INITIAL GUESS ***********************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Initial guess"
      write (*, 102)
      
      !> Set the occupation matrix, for RHF this is the identity matrix
      call set_n_occ(nbf, nel, n_occ)

      !> Print the occupation matrix manually
      write(*, '(A)') '', 'matrix : Initial n_occ:', ''
      if (present(print_level)) then
         do i = 1, nbf
               write(output_unit, '(I6, I7)') i, n_occ(i, i)
         end do
      end if

      !> Fock matrix
      allocate(F(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F failed."
      allocate(F_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_packed failed."
      allocate(F_prime(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime failed."
      allocate(F_prime_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime_packed failed."

      !>  Use Hcore as initial guess for Hamiltonian: F = H_0 = T + V
      F = T + V
      call pack_matrix(F, F_packed)

      !> Transform F to orthonormal basis
      F_prime = matmul(matmul(transpose(X), F), X)
      call pack_matrix(F_prime, F_prime_packed)

      !> Coeffs matrix
      allocate(C(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C failed."
      allocate(C_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_packed failed."
      allocate(C_prime(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_prime failed."

      !> Eigenvalues
      allocate(eps(nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of eps failed."

      !> Diagonalize initial Fock matrix to obtain initial orbital coefficients
      call solve_spev(F_prime_packed, eps, C_prime, alloc_stat)
      if (alloc_stat /= 0) error stop "Diagonalization of F failed."
     
      !> Print the initial guesses
      if (present(print_level)) then
         call write_matrix(F, name='Initial guess F')
         call write_matrix(F_prime, name='Initial guess F_prime')
         call write_matrix(C, name='Initial guess C')
         call write_matrix(C_prime, name='Initial guess C_prime')
         call write_vector(eps, name='Initial guess eps')
      end if

      !> Transform C to original basis
      C = matmul(X, C_prime)

      !> Density matrix
      allocate(P(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of P failed."

      !> Form density matrix
      P = matmul(matmul(C, n_occ), transpose(C))

      !> Print the density matrix
      if (present(print_level)) then
         call write_matrix(P, name='Initial P')
      end if

      !> Calculate initial HF energy
      call calc_hf_energy(T+V, F, P, nbf, escf)

      !> Print the initial HF energy
      if (present(print_level)) then
         write(*, *) ''
         write(*, '(A, F20.14)') 'E_HF_init :', escf
      end if

      !*********************************************************
      !******************* TWO-ELECTRON INTEGRALS **************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Two-electron integrals"
      write (*, 102)

      !> Two-electron integrals
      allocate(two_ints(nbf, nbf, nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of two_ints failed."

      !> Calculate two-electron integrals
      do i = 1, nbf
         do j = 1, i
         do k = 1, i
            do l = 1, merge(j, k, i == k)
                  call twoint( &
                     xyz(:, bf_atom_map(i)), xyz(:, bf_atom_map(i)), &
                     xyz(:, bf_atom_map(k)), xyz(:, bf_atom_map(k)), &
                     expnts(:, i), expnts(:, j), expnts(:, k), expnts(:, l), &
                     coeffs(:, i), coeffs(:, j), coeffs(:, k), coeffs(:, l), &
                     two_ints(i, j, k, l))

                  !> Symmetry
                  two_ints(j, i, k, l) = two_ints(i, j, k, l)
                  two_ints(i, j, l, k) = two_ints(i, j, k, l)
                  two_ints(j, i, l, k) = two_ints(i, j, k, l)
                  two_ints(k, l, i, j) = two_ints(i, j, k, l)
                  two_ints(l, k, i, j) = two_ints(i, j, k, l)
                  two_ints(k, l, j, i) = two_ints(i, j, k, l)
                  two_ints(l, k, j, i) = two_ints(i, j, k, l)
               end do
            end do
         end do
      end do

      !> Print the two-electron integrals
      if (present(print_level)) then
         write(*, '(A)') '', 'We dont print two electron integrals, they are too thicc'
      end if

      !*********************************************************
      !******************* SCF ITERATIONS **********************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "SCF iterations"
      write (*, 102)
      write(*, '(A)') '', ''
      
      !> Header
      if (present(print_level)) then
         write(*, '(A6, A20,A20, A20)') 'Iter', 'E_scf', 'E_tot', 'Delta'
      end if

      !> SCF iterations
      do iter = 1, MAX_SCF

         !> Construct the Fock matrix from the density matrix
         F = T + V
         do i = 1, nbf
            do j = 1, nbf
               do k = 1, nbf
                  do l = 1, nbf
                     ! F(i, j) = F(i, j) + P(k, l) * (two_ints(i, j, k, l) - 0.5_wp *  two_ints(i, l, k, j))
                     F(i, j) = F(i, j) + P(k, l)*( &
                        two_ints(i, j, l, k) - 0.5_wp*two_ints(i, l, j, k))
                  end do
               end do
            end do
         end do
         call pack_matrix(F, F_packed)

         !> Transform F to orthonormal basis
         F_prime = matmul(matmul(transpose(X), F), X)
         call pack_matrix(F_prime, F_prime_packed)

         !> Diagonalize Fock matrix to obtain new orbital coefficients
         call solve_spev(F_prime_packed, eps, C_prime, alloc_stat)
         if (alloc_stat /= 0) error stop "Diagonalization of F failed."

         !> Transform C to original basis
         C = matmul(X, C_prime)

         !> Form density matrix
         P = matmul(matmul(C, n_occ), transpose(C))

         !> Calculate the HF energy
         call calc_hf_energy(T+V , F, P, nbf, escf)

         !> Print the HF energy
         write(*, '(I6, F20.14, F20.14, F20.14)') iter, escf, escf + enn, escf - ehf

         !> Check for convergence
         if (abs(escf - ehf) < TOL_SCF) then
            ehf = escf
            write (*, 101)
            write (*, "(A, I3, A, F18.14)") "SCF converged in ", iter, " iterations with energy ", ehf
            write (*, 102)
            exit
         end if

         ! !> Debug print
         ! if (present(print_level)) then
         !    write (*, 101)
         !    write (*, *) "Scf iteration ", iter
         !    call write_matrix(F, name='F')
         !    call write_matrix(F_prime, name='F_prime')
         !    call write_vector(eps, name='eps')
         !    call write_matrix(P, name='P')
         !    write (*, 102)
         ! end if

         !> Update the HF energy
         ehf = escf

      end do

      !*********************************************************
      !******************* PROPERTY CALCULATIONS ***************
      !*********************************************************



   end subroutine scf_prog

   !> Calculate the symmetric orthonormalizer
   subroutine calc_symmetric_orthonormalizer(S_packed, nbf, X)
      
      implicit none

      !> Packed overlap matrix
      real(wp), intent(inout) :: S_packed(:)

      !> Number of basis functions
      integer, intent(in) :: nbf

      !> Symmetric orthonormalizer
      real(wp), dimension(:,:), allocatable, intent(out) :: X

      !> Eigenvalues and eigenvectors
      real(wp), dimension(:), allocatable :: S_diag_packed
      real(wp), dimension(:,:), allocatable :: S
      real(wp), dimension(:,:), allocatable :: s_to_the_minus_one_halfs
      real(wp), dimension(:,:), allocatable :: eigenvec
      real(wp), dimension(:,:), allocatable :: unitary
      real(wp), dimension(:), allocatable :: eigenval

      !> Stat
      integer :: alloc_stat, lapack_stat

      !> Loop variables
      integer :: i, j, k

      !> Allocate memory
      allocate(S_diag_packed(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of S_diag_packed failed."
      allocate(eigenval(nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of eigenval failed."
      allocate(eigenvec(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of eigenvec failed."

      !> Diagnalize S
      S_diag_packed = S_packed
      call solve_spev(S_diag_packed, eigenval, eigenvec, lapack_stat)
      if (lapack_stat /= 0) error stop "Diagonalization of S failed."
      
      !> Calculate s**−1/2
      allocate(s_to_the_minus_one_halfs(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of s_to_the_minus_one_halfs failed."

      s_to_the_minus_one_halfs = 0.0_wp
      do i = 1, nbf
         s_to_the_minus_one_halfs(i,i) = 1.0_wp / sqrt(eigenval(i))
      end do

      !> Calculate X = U s**−1/2 U**⊤
      allocate(X(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of X failed."

      X = matmul(matmul(eigenvec, s_to_the_minus_one_halfs), transpose(eigenvec))

      !> Check if X**T S X = 1
      allocate(S(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of S failed."
      allocate(unitary(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of X failed."

      !> Unpack S manually
      do i = 1, nbf
         do j = 1, i
            S(i, j) = S_packed((i-1)*i/2+j)
            S(j, i) = S(i, j)
         end do
      end do

      !> Calculate X**T S X
      unitary = matmul(matmul(transpose(X), S), X)

      if (is_mat_identity(unitary) .eqv. .false.) then
         call write_matrix(unitary, name='X^T S X')
         error stop "X**T S X is not the identity matrix!"
      end if

   end subroutine calc_symmetric_orthonormalizer

   !> Set the occupation matrix
   subroutine set_n_occ(nbf, nel, n_occ)
      implicit none
      integer, intent(in) :: nbf
      integer, intent(in) :: nel
      integer, intent(out), allocatable, dimension(:, :) :: n_occ

      ! variables for routine
      integer :: alloc_stat, i, temp

      allocate (n_occ(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 1

      n_occ = 0.0_wp
      temp = nel
      do i = 1, nbf
         n_occ(i, i) = 2
         temp = temp - 2
         if (temp == 0) then
         exit
         end if
      end do
   
   end subroutine set_n_occ

   !> Calculate HF energy
   subroutine calc_hf_energy(H_0, F, P, nbf, ehf)
      implicit none
      real(wp), intent(in) :: H_0(:,:)
      real(wp), intent(in) :: F(:,:)
      real(wp), intent(in) :: P(:,:)
      integer, intent(in) :: nbf
      real(wp), intent(out) :: ehf

      real(wp), allocatable :: temp(:,:)

      integer :: i, j

      !> Allocate memory
      allocate(temp(nbf, nbf))

      !> Calculate F = H_0 + F
      temp = matmul(H_0 + F, P)

      !> Calculate HF energy
      ehf = 0.0_wp
      do i = 1, nbf
         ehf = ehf + temp(i, i) * 0.5_wp
      end do

   end subroutine calc_hf_energy

end module scf_main