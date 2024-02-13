!> This is your module to write your very own SCF program.
module uhf_main
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
   public :: uhf_prog

contains

   !> This is the entry point to your program, do not modify the dummy arguments
   subroutine uhf_prog(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
                      ehf, MAX_SCF, TOL_SCF, print_level)
                      
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

      !> HF energy
      real(wp), intent(out) :: ehf
      real(wp) :: ehf_a, ehf_b

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
      real(wp) :: damp

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> alpha & beta
      real(wp), allocatable, dimension(:,:) :: S, T, V
      real(wp), allocatable, dimension(:) :: S_packed, T_packed, V_packed
      real(wp), allocatable, dimension(:,:) :: X

      !> alpha
      real(wp), allocatable, dimension(:) :: eps_a
      integer, dimension(:, :), allocatable :: n_occ_a
      real(wp), dimension(:, :), allocatable :: F_a
      real(wp), dimension(:, :), allocatable :: F_prime_a
      real(wp), dimension(:), allocatable :: F_packed_a
      real(wp), dimension(:), allocatable :: F_prime_packed_a
      real(wp), dimension(:, :), allocatable :: C_a
      real(wp), dimension(:), allocatable :: C_packed_a
      real(wp), dimension(:, :), allocatable :: C_prime_a
      real(wp), dimension(:), allocatable :: C_prime_packed_a
      real(wp), dimension(:, :), allocatable :: P_a

      !> beta
      real(wp), allocatable, dimension(:) :: eps_b
      integer, dimension(:, :), allocatable :: n_occ_b
      real(wp), dimension(:, :), allocatable :: F_b
      real(wp), dimension(:, :), allocatable :: F_prime_b
      real(wp), dimension(:), allocatable :: F_packed_b
      real(wp), dimension(:), allocatable :: F_prime_packed_b
      real(wp), dimension(:, :), allocatable :: C_b
      real(wp), dimension(:), allocatable :: C_packed_b
      real(wp), dimension(:, :), allocatable :: C_prime_b
      real(wp), dimension(:), allocatable :: C_prime_packed_b
      real(wp), dimension(:, :), allocatable :: P_b

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: alloc_stat
      real(wp), dimension(:, :, :, :), allocatable :: two_ints

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: i, j, k, l
      integer :: spin

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

      !> Nuclear repulsion energy enn = ZA * ZB / | R1 - R2 |
      enn = 0.0_wp
      do i = 1, nat - 1
         do j = i + 1, nat
         enn = enn + chrg(i)*chrg(j)/norm2(xyz(:, i) - xyz(:, j))
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
                     xyz(:, bf_atom_map(i)), xyz(:, bf_atom_map(j)), &
                     xyz(:, bf_atom_map(k)), xyz(:, bf_atom_map(l)), &
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
      !******************* INITIAL GUESS ***********************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "Initial guess"
      write (*, 102)
      
      !> Set the occupation matrix, for RHF this is the identity matrix
      call set_n_occ(nbf, nel, n_occ_a, n_occ_b)

      !> Print the occupation matrix manually
      write(*, '(A)') '', 'matrix : Initial n_occ:', ''
      if (present(print_level)) then
         do i = 1, nbf
               write(output_unit, '(I6, A3, I7)') i, "a", n_occ_a(i, i)
               write(output_unit, '(I6, A3, I7)') i, "b", n_occ_b(i, i)
         end do
      end if

      !> Density matrix
      allocate(P_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of P failed."
      allocate(P_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of P failed."

      !> Form density matrix
      P_a = matmul(matmul(X, n_occ_a), transpose(X))
      !> Break the symmetry
      P_b = 0.0_wp

      !> Fock matrix
      allocate(F_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F failed."
      allocate(F_prime_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime failed."
      allocate(F_prime_packed_a(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime_packed failed."
      allocate(F_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F failed."
      allocate(F_prime_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime failed."
      allocate(F_prime_packed_b(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of F_prime_packed failed."
      
      !>  Use Hcore as initial guess for Hamiltonian: F = H_0 = T + V
      F_a = T + V
      F_b = T + V
      do i = 1, nbf
         do j = 1, nbf
            do k = 1, nbf
               do l = 1, nbf
                  F_a(i, j) = F_a(i, j) + &
                              (P_a(k, l) + P_b(k, l))*two_ints(i, j, l, k) - &
                              P_a(k, l)*two_ints(i, l, j, k)
                  F_b(i, j) = F_b(i, j) + &
                              (P_a(k, l) + P_b(k, l))*two_ints(i, j, l, k) - &
                              P_b(k, l)*two_ints(i, l, j, k)
               end do
            end do
         end do
      end do

      !> Transform F to orthonormal basis
      F_prime_a = matmul(matmul(transpose(X), F_a), X)
      F_prime_b = matmul(matmul(transpose(X), F_b), X)
      call pack_matrix(F_prime_a, F_prime_packed_a)
      call pack_matrix(F_prime_b, F_prime_packed_b)

      !> Coeffs matrix
      allocate(C_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C failed."
      allocate(C_packed_a(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_packed failed."
      allocate(C_prime_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_prime failed."
      allocate(C_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C failed."
      allocate(C_packed_b(nbf*(nbf+1)/2), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_packed failed."
      allocate(C_prime_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of C_prime failed."

      !> Eigenvalues
      allocate(eps_a(nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of eps failed."
      allocate(eps_b(nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop "Allocation of eps failed."

      !> Diagonalize initial Fock matrix to obtain initial orbital coefficients
      call solve_spev(F_prime_packed_a, eps_a, C_prime_a, alloc_stat)
      if (alloc_stat /= 0) error stop "Diagonalization of F_a failed."
      call solve_spev(F_prime_packed_b, eps_b, C_prime_b, alloc_stat)
      if (alloc_stat /= 0) error stop "Diagonalization of F_b failed."
     
      !> Print the initial guesses
      if (present(print_level)) then
         call write_matrix(F_a, name='Initial guess F_a')
         call write_matrix(F_b, name='Initial guess F_b')
         call write_matrix(F_prime_a, name='Initial guess F_prime_a')
         call write_matrix(F_prime_b, name='Initial guess F_prime_b')
         call write_matrix(C_a, name='Initial guess C_a')
         call write_matrix(C_b, name='Initial guess C_b')
         call write_matrix(C_prime_a, name='Initial guess C_prime_a')
         call write_matrix(C_prime_b, name='Initial guess C_prime_b')
         call write_vector(eps_a, name='Initial guess eps_a')
         call write_vector(eps_b, name='Initial guess eps_b')
      end if

      !> Transform C to original basis
      C_a = matmul(X, C_prime_a)
      C_b = matmul(X, C_prime_b)

      !> Form density matrix
      P_a = matmul(matmul(C_a, n_occ_a), transpose(C_a))
      P_b = matmul(matmul(C_b, n_occ_b), transpose(C_b))

      !> Print the density matrix
      if (present(print_level)) then
         call write_matrix(P_a, name='Initial P_a')
         call write_matrix(P_b, name='Initial P_b')
      end if

      !> Calculate initial HF energy
      call calc_hf_energy(T+V, F_a, F_b, P_a, P_b, nbf, ehf)

      !> Print the initial HF energy
      if (present(print_level)) then
         write(*, *) ''
         write(*, '(A, 1F20.14)') 'E_HF_init :', ehf
      end if

   
      !*********************************************************
      !******************* SCF ITERATIONS **********************
      !*********************************************************

      write (*, 101)
      write (*, "(A)") "SCF iterations"
      write (*, 102)
      write(*, '(A)') '', ''
      
      !> Damping factor
      damp = 0.1_wp

      !> Header
      write(*, '(A6, 2A20, A8)') 'Iter', 'E_scf', 'Delta', 'Damp'
      write(*, '(A6, 2A20, A8)') '----', '------------------', '------------------', '------'
      write(*, '(I6, 2F20.14, F8.4)') 0, ehf, 0.0_wp, damp

      !> SCF iterations
      do iter = 1, MAX_SCF

         !> Construct the Fock matrices from the density matrices
         F_a = T + V
         F_b = T + V
         do i = 1, nbf
            do j = 1, nbf
               do k = 1, nbf
                  do l = 1, nbf
                     F_a(i, j) = F_a(i, j) + &
                                 (P_a(k, l) + P_b(k, l))*two_ints(i, j, l, k) - &
                                 P_a(k, l)*two_ints(i, l, j, k)
                     F_b(i, j) = F_b(i, j) + &
                                 (P_a(k, l) + P_b(k, l))*two_ints(i, j, l, k) - &
                                 P_b(k, l)*two_ints(i, l, j, k)
                  end do
               end do
            end do
         end do

         !> Transform F to orthonormal basis with damping
         F_prime_a = damp * F_prime_a + (1 - damp) * matmul(matmul(transpose(X), F_a), X)
         F_prime_b = damp * F_prime_b + (1 - damp) * matmul(matmul(transpose(X), F_b), X)
         call pack_matrix(F_prime_a, F_prime_packed_a)
         call pack_matrix(F_prime_b, F_prime_packed_b)

         !> Diagonalize Fock matrix to obtain new orbital coefficients
         call solve_spev(F_prime_packed_a, eps_a, C_prime_a, alloc_stat)
         if (alloc_stat /= 0) error stop "Diagonalization of F failed."
         call solve_spev(F_prime_packed_b, eps_b, C_prime_b, alloc_stat)
         if (alloc_stat /= 0) error stop "Diagonalization of F failed."

         !> Transform C to original basis
         C_a = matmul(X, C_prime_a)
         C_b = matmul(X, C_prime_b)

         !> Form density matrix
         P_a = matmul(matmul(C_a, n_occ_a), transpose(C_a))
         P_b = matmul(matmul(C_b, n_occ_b), transpose(C_b))

         !> Calculate the HF energy
         call calc_hf_energy(T+V , F_a, F_b, P_a, P_b, nbf, escf)

         !> Print the HF energy
         write(*, '(I6, 2F20.14, F8.4)') iter, escf, escf - ehf, damp

         !> Check for convergence
         if (abs(escf - ehf) < TOL_SCF) then
            ehf = escf
            write (*, 101)
            write (*, "(A, I3, A)") "SCF converged in ", iter, " iterations"
            write (*, 102)
            exit
         end if

         !> Check if we have reached the maximum number of iterations
         if (iter == MAX_SCF) then
            write (*, 101)
            write (*, "(A)") "SCF did not converge! Maximum number of iterations reached."
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

      !> Print the final energies
      write (*, 101)
      write (*, "(A)") "Final energies:"
      write(*, '(A, F20.14)') '  E_HF    :', ehf
      write(*, '(A, F20.14)') '  V_nn    :', enn
      write(*, '(A, F20.14)') '  E_total :', ehf + enn
      write (*, 102)

      !*********************************************************
      !******************* SPIN CONTAMINATION ******************
      !*********************************************************

   end subroutine uhf_prog

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
   subroutine set_n_occ(nbf, nel, n_occ_a, n_occ_b)
      implicit none
      integer, intent(in) :: nbf
      integer, intent(in) :: nel
      integer, intent(out), allocatable, dimension(:, :) :: n_occ_a, n_occ_b

      ! variables for routine
      integer :: alloc_stat, i, temp

      allocate (n_occ_a(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 1
      allocate (n_occ_b(nbf, nbf), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 1

      n_occ_a = 0.0_wp
      n_occ_b = 0.0_wp

      temp = nel

      do i = 1, nbf

         n_occ_a(i, i) = 1
         temp = temp - 1
         
         if (temp == 0) then
         exit
         end if
      
         n_occ_b(i, i) = 1
         temp = temp - 1
         
         if (temp == 0) then
         exit
         end if

      end do
   
   end subroutine set_n_occ

   !> Calculate HF energy
   subroutine calc_hf_energy(H_0, F_a, F_b, P_a, P_b, nbf, ehf)
      implicit none
      real(wp), intent(in) :: H_0(:,:)
      real(wp), intent(in) :: F_a(:,:)
      real(wp), intent(in) :: F_b(:,:)
      real(wp), intent(in) :: P_a(:,:)
      real(wp), intent(in) :: P_b(:,:)
      integer, intent(in) :: nbf
      real(wp), intent(out) :: ehf

      real(wp), allocatable :: temp_a(:,:)
      real(wp), allocatable :: temp_b(:,:)

      integer :: i, j

      !> Allocate memory
      allocate(temp_a(nbf, nbf))
      allocate(temp_b(nbf, nbf))

      !> Calculate F = H_0 + F
      temp_a = matmul(H_0 + F_a, P_a)
      temp_b = matmul(H_0 + F_b, P_b)

      !> Calculate HF energy
      ehf = 0.0_wp
      do i = 1, nbf
         ehf = ehf + 0.5_wp * ( temp_a(i, i) + temp_b(i, i) )
      end do

   end subroutine calc_hf_energy

end module uhf_main