!> This is the actual main program, we provide a wrapper around the code you
!  are writing here, so you can skip parts of the necessary IO.
program main_prog
  !> Include standard Fortran environment for IO
  use iso_fortran_env, only: input_unit, error_unit
  !> We use our own read_argument function to obtain a command line arguments
  use io_tools, only: read_argument

  !> Import input file reader
  use input_reader, only: read_file

  !> prints a matrix quantity to screen
  use print_matrix

  !> Import rhf program
  use scf_main, only: scf_prog
!   use scf_main, only: opt_coords
!   use scf_main, only: opt_expnts

!   !> Import uhf program
!   use uhf, only: uhf_prog
  
!   !> convenient wrapper for root mean square
!   use array_funcs, only: rms, rms_vec

  !> Always declare everything explicitly
  implicit none

  !> amount of info printed during program run (0, 1, 2, 3)
  integer, parameter :: print_level = 3

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  !> Name of the input file
  character(len=:), allocatable :: input_file

  !> Does the provided input file exist
  logical :: exist

  !> Unit for reading the input from
  integer :: io_unit

  !> Error status from opening file for IO
  integer :: error

  !> Error message from opening file for IO
  character(len=256) :: message

  !> Number of atoms
  integer :: nat
  !> Number of electrons
  integer :: nel
  !> Number of basis functions
  integer :: nbf
  !> boolean for optimization (1 -> yes, 0 -> no opt)
  integer :: do_opt
  !> boolean for mp2 (1 -> yes, 0 -> no opt)
  integer :: do_mp2
  !> boolean for UHF (1 -> yes, 0 -> no opt)
  integer :: do_uhf

  !> Atom coordinates of the system, all distances in bohr
  real(wp), allocatable, dimension(:, :) :: xyz

  !> Nuclear charges
  real(wp), allocatable :: chrg(:)

  !> Slater expnts of basis functions
  real(wp), allocatable, dimension(:) :: zeta

  !> maps basis function to atom
  !> e.g.: (1, 1, 2) -> 1st bf from A1, 2nd bf from A1, 3rd bf from A3
  integer, allocatable :: bf_atom_map(:)

  !> number of primitive Gaussians in slater expansion (STO-nG)
  integer, parameter :: ng = 6

  !> HF energy
  real(wp) :: ehf

  !> convergence criteria
  real(wp), parameter :: TOL_SCF = 1.0e-8_wp
  real(wp), parameter :: TOL_OPT = 1.0e-4_wp
  integer, parameter :: MAX_SCF = 100
  integer, parameter :: MAX_OPT_GEOM = 1000
  integer, parameter :: MAX_OPT_EXP = 1000
  real(wp), parameter :: COORD_STEP_SIZE = 0.1_wp
  real(wp), parameter :: EXPNTS_STEP_SIZE = 0.01_wp
  real(wp), parameter :: ETA = 0.5_wp

  !> This code snippet optionally opens a file or allows reading from STDIN
  !> In case of no command line arguments, you read from the STDIN
  if (command_argument_count() == 0) then
    io_unit = input_unit
  else
    !> In case there are command line arguments we obtain the first one
    call read_argument(1, input_file)
    !> Check if the argument corresponds to an existing file
    inquire (file=input_file, exist=exist)
    !> The file does not exist, we will return with a meaningful error message
    if (.not. exist) then
      write (error_unit, '("ERROR:", 1x, a)') &
          & "The input file '"//input_file//"' does not exist"
      error stop 1
    end if
    !> If the file exist, we open it to a new unit an pass it to the scf
    open (file=input_file, newunit=io_unit, status='old', &
        & iomsg=message, iostat=error)
    if (error /= 0) then
      write (error_unit, '("ERROR:", 1x, a)') trim(message)
      error stop 1
    end if
  end if

  write (*, "(A24,A,A2)") "Running calculation on '", input_file, "'."

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

101 format(//, "---------------------------------------------------------------")
102 format("---------------------------------------------------------------")

  !*********************************************************
  !********************* READING INPUT *********************
  !*********************************************************

  !> init all variables from input file
  write (*, 101)
  write (*, "(A)") "Reading input file"
  call read_file(io_unit, xyz, chrg, zeta, bf_atom_map, nat, nel, nbf, &
                 do_opt, do_mp2, do_uhf)
  write (*, 102)

  write (*, 101)
  write (*, "(A)") "Settings"
  write (*, 102)

  write (*, "(A)") "Molecule"
  write (*, "(A20,15X,I1)") " Number of atoms ...", nat
  write (*, "(A24,11X,I1)") " Number of electrons ...", nel

  write (*, "(A)") "Basis set"
  write (*, "(A31,4X,I1)") " Number of Slater functions ...", nbf
  write (*, "(A18,17X,A4,I1,A1)") " Gaussian Type ...", "STO-", ng, "G"

  write (*, "(A)") "SCF"
  write (*, "(A22,13X)", advance="no") " Wavefunction type ..."
  if (do_uhf == 1 .or. modulo(nel, 2) /= 0) then
    write (*, "(A)") "UHF"
  else
    write (*, "(A)") "RHF"
  end if
  write (*, "(A26,9X,A)") " Convergence criterion ...", "Energy"
  write (*, "(A26,9X,ES7.1)") " Convergence threshold ...", TOL_SCF
  write (*, "(A33,2X,I3)") " Maximal number of iterations ...", MAX_SCF

  if (do_opt == 1) then
    write (*, "(A)") "GEOMETRY OPTIMIZATION"
    write (*, "(A32,3X,ES7.1)") " Convergence threshold (RMS) ...", TOL_OPT
    write (*, "(A33,2X,I4)") " Maximal number of iterations ...", MAX_OPT_GEOM
    write (*, "(A26,9X,F4.2)") " Step size coordinates ...", COORD_STEP_SIZE
    write (*, "(A18,17X,F4.2)") " Learprint_levelning rate ...", eta
    write (*, "(A18,17X,F4.2)") " Learprint_levelning rate ...", eta
  end if
  if (do_opt == 2) then
    write (*, "(A)") "EXPONENT OPTIMIZATION"
    write (*, "(A32,3X,ES7.1)") " Convergence threshold (RMS) ...", TOL_OPT
    write (*, "(A33,2X,I4)") " Maximal number of iterations ...", MAX_OPT_EXP
    write (*, "(A25,10X,F4.2)") " Step size expontents ...", EXPNTS_STEP_SIZE
    write (*, "(A18,17X,F4.2)") " Learning rate ...", eta
  end if

  if (print_level >= 3) then
    write (*, 101)
    write (*, "(A)") "Settings: extended output"
    write (*, 102)
    call write_matrix(xyz, name='Coordinates')
    call write_vector(chrg, name="Charges")
    call write_vector(real(bf_atom_map, wp), name="Mapping: basis function -> atom")
  end if

  !*********************************************************
  !********************* MAIN PROGRAM **********************
  !*********************************************************

  !> UHF if specified or number of electrons odd
  if (do_uhf == 1 .or. modulo(nel, 2) /= 0) then
    ! call uhf_prog(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
                !   ehf, MAX_SCF, TOL_SCF, print_level)
    write (*, 101)
    write (*, "(A)") "UHF not implemented yet"
    write (*, 102)
    error stop 1

    !> RHF
  else
    ! !> geometry optimization
    ! if (do_opt == 1 .or. do_opt == 3) then
    !   write (*, 101)
    !   write (*, "(A)") "Geometry optimization"
    !   write (*, 102)

    !   !> xyz is modified (inout!)
    !   call opt_coords(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
    !                   MAX_SCF, TOL_SCF, MAX_OPT_GEOM, TOL_OPT, ETA, &
    !                   COORD_STEP_SIZE, print_level)
    ! end if

!     !> Slater exponent optimization
!     if (do_opt == 2 .or. do_opt == 3) then
!       write (*, 101)
!       write (*, "(A)") "Exponent optimization"
!       write (*, 102)

!       !> zeta is modified (inout!)
!       call opt_expnts(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
!                       MAX_SCF, TOL_SCF, MAX_OPT_EXP, TOL_OPT, ETA, &
!                       EXPNTS_STEP_SIZE, print_level)
!     end if

    !> final energy calculation (uses optimized values)
    write (*, 101)
    write (*, "(A)") "RHF calculation"
    write (*, 102)
    call scf_prog(nat, nel, nbf, ng, xyz, chrg, zeta, bf_atom_map, &
                  ehf, MAX_SCF, TOL_SCF, do_mp2, print_level)

  end if

  !*********************************************************
  !************************ CLEAN UP ***********************
  !*********************************************************

  deallocate (zeta)
  deallocate (bf_atom_map)
  deallocate (chrg)
  deallocate (xyz)

  !> If we have opened a file, we have to cleanup now, so we close it again
  if (io_unit /= input_unit) then
    close (io_unit)
  end if

end program main_prog