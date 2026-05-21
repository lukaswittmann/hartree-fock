module tools

   use print_matrix

   implicit none

   integer, parameter :: wp = selected_real_kind(15)

   contains

   subroutine pack_matrix(mat, mat_packed)
      implicit none
      real(wp), intent(in), dimension(:, :), allocatable :: mat
      real(wp), intent(out), dimension(:), allocatable :: mat_packed

      ! variables for calculation
      integer :: i, j, counter, alloc_stat
      integer :: dim_packed, dim
      integer, dimension(2) :: dims

      !> get dim of matrix
      dims = shape(mat)
      dim = dims(1)

      if (is_mat_symmetric(mat) .eqv. .false.) then
      write (*, *) "failed"
      write (*, "(A,A)") "", "-> Could not pack matrix. Not symmetric!"
      call write_matrix(mat, name="mat")
      error stop 1
      end if

      !> formula for dimension of packed matrix
      dim_packed = dim*(1 + dim)/2

      !> allocate with calculated dimension
      allocate (mat_packed(dim_packed), stat=alloc_stat)
      if (alloc_stat /= 0) error stop 1

      !> packing
      counter = 1
      do i = 1, dim
      do j = 1, i
         mat_packed(counter) = mat(i, j)
         counter = counter + 1
      end do
      end do

   end subroutine pack_matrix

   logical function is_mat_symmetric(mat) result(is_symmetric)
      implicit none
      intrinsic :: abs

      real(wp), intent(in), dimension(:, :):: mat
      integer :: num_cols, num_rows
      integer, dimension(2) :: shape_mat
      integer :: i, j

      !> thresholds for two reals being equal
      real(wp), parameter :: TOL_EQ = 1.0e-8_wp

      !> get shape of matrix
      shape_mat = shape(mat)
      num_rows = shape_mat(1)
      num_cols = shape_mat(2)

      if (num_cols /= num_rows) then
      is_symmetric = .false.
      return
      end if

      is_symmetric = .true.
      outer: do i = 1, num_rows
      inner: do j = 1, i
         if (i /= j) then
            if (abs(mat(i, j) - mat(j, i)) > TOL_EQ) then
            is_symmetric = .false.
            exit outer
            end if
         end if
      end do inner
      end do outer
   end function is_mat_symmetric

   logical function is_mat_identity(mat) result(is_identity)
      implicit none
      intrinsic :: abs

      real(wp), intent(in), dimension(:, :), allocatable :: mat
      integer :: num_cols, num_rows
      integer, dimension(2) :: shape_mat
      integer :: i, j

      !> thresholds for two reals being equal
      real(wp), parameter :: TOL_EQ = 1.0e-8_wp

      !> get shape of matrix
      shape_mat = shape(mat)
      num_rows = shape_mat(1)
      num_cols = shape_mat(2)

      if (num_cols /= num_rows) then
      is_identity = .false.
      else
      is_identity = .true.

      outer: do i = 1, num_rows
         inner: do j = 1, num_rows
            if (i == j) then
            if (abs(mat(i, i) - 1.0_wp) > TOL_EQ) then
               is_identity = .false.
               exit outer
            end if
            else
            if (abs(mat(i, j)) > TOL_EQ) then
               is_identity = .false.
               exit outer
            end if
            end if
         end do inner
      end do outer

      end if
   end function is_mat_identity

   integer function get_dim_mat(mat) result(dim)
      implicit none
      intrinsic :: abs

      real(wp), intent(in), dimension(:, :), allocatable :: mat
      integer :: num_cols, num_rows
      integer, dimension(2) :: shape_mat

      !> get shape of matrix
      shape_mat = shape(mat)
      num_rows = shape_mat(1)
      num_cols = shape_mat(2)

      !> only for symmetric matrices
      if (num_cols /= num_rows) then
      write (*, *) "Matrix not symmetric. Aborting..."
      error stop 1
      end if

      dim = num_cols
   end function get_dim_mat

   real(wp) function rms(mat) result(val)
      implicit none
      intrinsic :: sqrt, sum

      real(wp), intent(in), dimension(:, :), allocatable :: mat

      val = sqrt(sum(mat**2)/size(mat))
   end function rms

   real(wp) function rms_vec(vec) result(val)
      implicit none
      intrinsic :: sqrt, sum

      real(wp), intent(in), dimension(:), allocatable :: vec

      val = sqrt(sum(vec**2)/size(vec))
   end function rms_vec

end module tools
