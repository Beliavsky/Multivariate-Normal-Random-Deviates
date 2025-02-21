program test_mv_normal
!> Test program to generate and verify multivariate normal deviates.
!! Generates samples from a multivariate normal distribution
!! with zero mean, unit standard deviation, and a specified correlation matrix.
!! Verifies that empirical means, standard deviations, and correlations match
!! theoretical values
use       kind_mod, only: dp
use  mv_normal_mod, only: generate_mv_normal, cov_to_ltrans, check_cov_info
use     random_mod, only: random_normal
use statistics_mod, only: mean, sd, corr_mat, print_means_and_sds
implicit none
real(kind=dp), allocatable :: xcorr(:,:), xcov(:,:), deviates(:,:), &
emp_corr(:,:), means(:), sds(:), ltrans(:,:)
integer :: n_vars, n_samples, info, i, j

! Test parameters
n_vars = 3
n_samples = 10**6
print *, "#obs:", n_samples
allocate(xcorr(n_vars, n_vars), xcov(n_vars, n_vars), deviates(n_samples, n_vars), &
ltrans(n_vars, n_vars))
means = real([-5, 0, 5], kind=dp)
sds = real([10, 20, 30], kind=dp)
! Sample correlation matrix (true correlation)
xcorr = reshape([1.0_dp, 0.5_dp, 0.3_dp,  &
                 0.5_dp, 1.0_dp, 0.4_dp,  &
                 0.3_dp, 0.4_dp, 1.0_dp], &
                 [n_vars, n_vars], order=[2,1])
do i=1,n_vars
   do j=1,n_vars
      xcov(i,j) = xcorr(i,j) * sds(i) * sds(j)
   end do
end do
! Generate multivariate normal deviates
call generate_mv_normal(xcov, deviates, info, means=means)
call check_cov_info(info) ! Check generation status

print "(/,a)", "true means and standard deviations"
do i = 1, n_vars
   print "(A,I2,A,F12.8,A,F12.8)", "Variable ", i, &
         ": Mean = ", means(i), ", StdDev = ", sds(i)
end do

call print_means_and_sds(deviates) ! Verify means and standard deviations

! Compare true and empirical correlations
print "(/,a)", "true correlation matrix:"
do i = 1, n_vars
   print "(*(F12.6))", xcorr(i,:)
end do

emp_corr = corr_mat(deviates)     ! Compute empirical correlation matrix
print "(/,a)", "empirical correlation matrix:"
do i = 1, n_vars
   print "(*(F12.6))", emp_corr(i,:)
end do

! Compute maximum absolute difference
print "(/, A,F12.6)", "maxval(abs(xcorr - emp_corr)): ", maxval(abs(xcorr - emp_corr))

! Generate multivariate normal deviates by multiplying uncorrelated deviates by the Cholesky factor
call cov_to_ltrans(xcov, ltrans, info)
do i=1,n_samples
   deviates(i, :) = matmul(random_normal(n_vars), ltrans)
end do

call print_means_and_sds(deviates) ! Verify means and standard deviations

emp_corr = corr_mat(deviates)     ! Compute empirical correlation matrix
print "(/,a)", "empirical correlation matrix:"
do i = 1, n_vars
   print "(*(F12.6))", emp_corr(i,:)
end do
! Compute maximum absolute difference
print "(/, A,F12.6)", "maxval(abs(xcorr - emp_corr)): ", maxval(abs(xcorr - emp_corr))

print "(/,a)", "(4) finished xmv_normal_test.f90"
end program test_mv_normal
