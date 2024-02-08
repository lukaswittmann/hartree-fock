module input_reader
  implicit none

  !> Selecting double precision real number
  integer, parameter :: wp = selected_real_kind(15)

contains
  subroutine read_file(input, xyz, chrg, zeta, bf_atom_map, nat, nel, nbf,&
     do_opt, do_mp2, do_uhf)
    implicit none

    ! in and out variables
    integer, intent(in) :: input
    real(wp), allocatable, intent(out) :: xyz(:, :)
    real(wp), allocatable, intent(out) :: chrg(:)
    real(wp), allocatable, intent(out) :: zeta(:)
    integer, allocatable, intent(out) :: bf_atom_map(:)
    integer, intent(out) :: nat
    integer, intent(out) :: nel
    integer, intent(out) :: nbf
    integer, intent(out) :: do_opt, do_mp2, do_uhf

    ! variables for calculation
    integer :: io_stat, alloc_stat, i, j
    integer :: nbf_of_atom
    integer :: counter
    counter = 1


    ! ! read first line containing info for tasks
    ! read (input, *, iostat=io_stat) do_opt, do_mp2, do_uhf
    ! if (io_stat /= 0) then
    !   write(*,*) 'Error reading first line of input file' // &
    !     ' (do_opt, do_mp2, do_uhf)'
    !   error stop 1
    ! end if

    ! read second line containing info for allocation
    read (input, *, iostat=io_stat) nat, nel, nbf
    if (io_stat /= 0) then
      write(*,*) 'Error reading second line of input file' // &
      ' (nat, nel, nbf)'
      error stop 1
    end if

    ! memory allocation of arrays
    allocate (zeta(nbf), stat=alloc_stat)
    if (alloc_stat /= 0) then
      write(*,*) 'Error allocating memory for zeta'
      error stop 1
    end if
    allocate (chrg(nat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      write(*,*) 'Error allocating memory for chrg'
      error stop 1
    end if
    allocate (xyz(3, nat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      write(*,*) 'Error allocating memory for xyz'
      error stop 1
    end if
    allocate (bf_atom_map(nbf), stat=alloc_stat)
    if (alloc_stat /= 0) then
      write(*,*) 'Error allocating memory for bf_atom_map'
      error stop 1
    end if

    ! read rest of file
    ! iterate over number of atoms
    do i = 1, nat
      read (input, *, iostat=io_stat) xyz(:, i), chrg(i), nbf_of_atom
      if (io_stat /= 0) exit

      ! iterate over number of basis functions of this atom
      do j = 1, nbf_of_atom
        read (input, *, iostat=io_stat) zeta(counter)
        if (io_stat /= 0) exit
        bf_atom_map(counter) = i
        counter = counter + 1
      end do
    end do
  end subroutine

end module input_reader