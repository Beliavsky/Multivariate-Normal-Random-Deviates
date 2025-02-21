module linear_algebra_mod
use kind_mod, only: dp
use constants_mod, only: pi
implicit none
private
public :: cholesky, eigenvalues_jacobi
real(kind=dp), parameter :: TOL = 1.0e-12_dp  ! Tolerance for checks
contains
subroutine cholesky(A, L, info)
!> Compute the Cholesky decomposition of a symmetric positive definite matrix.
!! Returns a lower triangular matrix L such that A = L * L^T.
real(kind=dp), intent(in) :: A(:,:)
real(kind=dp), intent(out) :: L(:,:)
integer, intent(out) :: info
integer :: n, i, j, k
real(kind=dp) :: sum_val

n = size(A, 1)
L = 0.0_dp
info = 0

do j = 1, n
   sum_val = A(j,j)
   do k = 1, j - 1
      sum_val = sum_val - L(j,k)**2
   end do
   if (sum_val <= 0.0_dp) then
      info = 1
      return
   end if
   L(j,j) = sqrt(sum_val)
   do i = j + 1, n
      sum_val = A(i,j)
      do k = 1, j - 1
         sum_val = sum_val - L(i,k) * L(j,k)
      end do
      L(i,j) = sum_val / L(j,j)
   end do
end do
end subroutine cholesky

subroutine eigenvalues_jacobi(matrix, eigenvalues, max_iter, info)
!> Compute eigenvalues of a symmetric matrix using the Jacobi method.
!! Iteratively applies Givens rotations to diagonalize the matrix, with
!! eigenvalues extracted from the diagonal upon convergence.
real(kind=dp), intent(inout) :: matrix(:,:)
real(kind=dp), intent(out) :: eigenvalues(:)
integer, intent(in) :: max_iter
integer, intent(out) :: info
integer :: n, i, j, p, q, iter
real(kind=dp) :: theta, t, c, s, tau, off_diagonal_sum, prev_sum
real(kind=dp), allocatable :: temp_mat(:,:)

n = size(matrix, 1)
if (size(matrix, 2) /= n .or. size(eigenvalues) /= n) then
   info = -1
   return
end if
temp_mat = matrix
info = 0
iter = 0
do while (iter < max_iter)
   iter = iter + 1
   off_diagonal_sum = 0.0_dp
   do i = 1, n
      do j = i+1, n
         off_diagonal_sum = off_diagonal_sum + temp_mat(i,j)**2
      end do
   end do
   off_diagonal_sum = sqrt(2.0_dp * off_diagonal_sum)
   if (iter > 1 .and. abs(off_diagonal_sum - prev_sum) < TOL) exit
   prev_sum = off_diagonal_sum
   p = 1
   q = 2
   do i = 1, n
      do j = i+1, n
         if (abs(temp_mat(i,j)) > abs(temp_mat(p,q))) then
            p = i
            q = j
         end if
      end do
   end do
   if (off_diagonal_sum < TOL) exit
   if (abs(temp_mat(p,p) - temp_mat(q,q)) < TOL) then
      theta = pi / 4.0_dp
   else
      tau = (temp_mat(q,q) - temp_mat(p,p)) / (2.0_dp * temp_mat(p,q))
      t = sign(1.0_dp, tau) / (abs(tau) + sqrt(1.0_dp + tau**2))
      theta = atan(t)
   end if
   c = cos(theta)
   s = sin(theta)
   do i = 1, n
      if (i /= p .and. i /= q) then
         temp_mat(i,p) = c * temp_mat(i,p) + s * temp_mat(i,q)
         temp_mat(p,i) = temp_mat(i,p)
         temp_mat(i,q) = -s * temp_mat(i,p) + c * temp_mat(i,q)
         temp_mat(q,i) = temp_mat(i,q)
      end if
   end do
   temp_mat(p,p) = c**2 * temp_mat(p,p) + 2.0_dp * c * s * temp_mat(p,q) + s**2 * temp_mat(q,q)
   temp_mat(q,q) = s**2 * temp_mat(p,p) - 2.0_dp * c * s * temp_mat(p,q) + c**2 * temp_mat(q,q)
   temp_mat(p,q) = 0.0_dp
   temp_mat(q,p) = 0.0_dp
end do
if (iter >= max_iter) info = 1
do i = 1, n
   eigenvalues(i) = temp_mat(i,i)
end do
matrix = temp_mat
end subroutine eigenvalues_jacobi
end module linear_algebra_mod
