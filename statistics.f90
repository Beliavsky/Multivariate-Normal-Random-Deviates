module statistics_mod
use kind_mod, only: dp
implicit none
private
public :: corr_mat, mean, sd, print_means_and_sds
real(kind=dp), parameter :: TOL = 1.0e-12_dp  ! Tolerance for checks

contains

function mean(x) result(m)
! Compute mean of a vector
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: m
m = sum(x) / size(x)
end function mean

function sd(x) result(s)
! Compute standard deviation of a vector
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: s, m
integer :: n
n = size(x)
if (n < 2) then
   s = 0.0_dp
   return
end if
m = mean(x)
s = sqrt(sum((x - m)**2) / (n - 1))
end function sd

function corr_mat(x) result(corr)
! Compute correlation matrix from data (cols = variables, rows = samples)
real(kind=dp), intent(in) :: x(:,:)
real(kind=dp), allocatable :: corr(:,:)
integer :: n_vars, n_samples, i, j
real(kind=dp) :: m_i, m_j, s_i, s_j, cov
n_vars = size(x, 2)     ! Variables are now columns
n_samples = size(x, 1)  ! Samples are now rows
allocate(corr(n_vars, n_vars))
do i = 1, n_vars
   m_i = mean(x(:,i))  ! Mean over column i
   s_i = sd(x(:,i))    ! Std dev over column i
   if (abs(s_i) < TOL) s_i = 1.0_dp  ! Avoid division by zero
   do j = 1, n_vars
      if (i == j) then
         corr(i,j) = 1.0_dp
      else
         m_j = mean(x(:,j))  ! Mean over column j
         s_j = sd(x(:,j))    ! Std dev over column j
         if (abs(s_j) < TOL) s_j = 1.0_dp
         cov = sum((x(:,i) - m_i) * (x(:,j) - m_j)) / (n_samples - 1)
         corr(i,j) = cov / (s_i * s_j)
      end if
   end do
end do
end function corr_mat

subroutine print_means_and_sds(x)
! print means and standard deviations of columns of x(:,:)
real(kind=dp), intent(in) :: x(:,:)
integer :: i
print "(/,a)", "empirical means and standard deviations"
do i = 1, size(x, 2)
   print "(A,I2,A,F12.8,A,F12.8)", "Variable ", i, &
   ": Mean = ", mean(x(:,i)), ", StdDev = ", sd(x(:,i))
end do
end subroutine print_means_and_sds

end module statistics_mod
