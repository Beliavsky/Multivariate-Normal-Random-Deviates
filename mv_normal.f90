module mv_normal_mod
    use      kind_mod, only: dp
    use    random_mod, only: random_normal
    use constants_mod, only: pi
    use linear_algebra_mod, only: cholesky, eigenvalues_jacobi
    implicit none
    private
    public :: generate_mv_normal_corr, generate_mv_normal, cov_to_ltrans, &
              check_cov_info
    real(kind=dp), parameter :: TOL = 1.0e-12_dp  ! Tolerance for checks

contains

    subroutine generate_mv_normal_corr(corr_mat, deviates, info, means, sds)
        !> Generate multivariate normal deviates with specified correlation.
        !! Validates the correlation matrix, computes its Cholesky decomposition, and
        !! transforms standard normal deviates to match the correlation structure.
        !! Output deviates have shape (n_samples, n_vars), with variables in columns.
        !! Optional means and sds set column means and standard deviations if provided.
        real(kind=dp), intent(in) :: corr_mat(:,:)
        real(kind=dp), intent(out) :: deviates(:,:)  ! (n_samples, n_vars)
        integer, intent(out) :: info
        real(kind=dp), intent(in), optional :: means(:)  ! Column means
        real(kind=dp), intent(in), optional :: sds(:)    ! Column standard deviations
        integer :: n_samples, n_vars, i, j
        real(kind=dp), allocatable :: L(:,:), z(:,:), temp_corr(:,:), eigenvalues(:), LLT(:,:)
        real(kind=dp), allocatable :: sd_factor(:)

        n_samples = size(deviates, 1)
        n_vars = size(corr_mat, 1)
        if (size(corr_mat, 2) /= n_vars .or. size(deviates, 2) /= n_vars) then
            info = 1
            return
        end if
        if (present(means) .and. size(means) /= n_vars) then
            info = 9
            return
        end if
        if (present(sds)) then
            if (size(sds) /= n_vars) then
                info = 10
                return
            end if
            if (any(sds <= 0.0_dp)) then
                info = 11
                return
            end if
        end if

        info = 0
        do i = 1, n_vars
            do j = 1, i - 1
                if (abs(corr_mat(i,j) - corr_mat(j,i)) > TOL) then
                    info = 2
                    return
                end if
            end do
        end do
        do i = 1, n_vars
            if (abs(corr_mat(i,i) - 1.0_dp) > TOL) then
                info = 3
                return
            end if
        end do
        if (any(abs(corr_mat) > 1.0_dp)) then
            info = 4
            return
        end if

        allocate(temp_corr(n_vars, n_vars), eigenvalues(n_vars))
        temp_corr = corr_mat
        call eigenvalues_jacobi(temp_corr, eigenvalues, 1000, info)
        if (info /= 0) then
            info = 5
            deallocate(temp_corr, eigenvalues)
            return
        end if
        do i = 1, n_vars
            if (eigenvalues(i) < -TOL) then
                info = 6
                deallocate(temp_corr, eigenvalues)
                return
            end if
        end do
        deallocate(eigenvalues, temp_corr)

        allocate(L(n_vars, n_vars))
        call cholesky(corr_mat, L, info)
        if (info /= 0) then
            info = 7
            deallocate(L)
            return
        end if

        allocate(LLT(n_vars, n_vars))
        LLT = matmul(L, transpose(L))
        do i = 1, n_vars
            do j = 1, n_vars
                if (abs(LLT(i,j) - corr_mat(i,j)) > 1.0e-6_dp) then
                    print *, "Cholesky error at (", i, ",", j, "): L*L^T =", LLT(i,j), &
                             " vs corr_mat =", corr_mat(i,j)
                    info = 8
                    deallocate(L, LLT)
                    return
                end if
            end do
        end do
        deallocate(LLT)

        allocate(z(n_samples, n_vars))
        z = random_normal(n_samples, n_vars)
        deviates = matmul(z, transpose(L))

        allocate(sd_factor(n_vars))
        sd_factor = 1.0_dp
        if (present(sds)) sd_factor = sds
        do i = 1, n_vars
            deviates(:,i) = deviates(:,i) * sd_factor(i)
        end do
        deallocate(sd_factor)

        if (present(means)) then
            do i = 1, n_vars
                deviates(:,i) = deviates(:,i) + means(i)
            end do
        end if

        deallocate(L, z)
    end subroutine generate_mv_normal_corr

    subroutine generate_mv_normal(cov_mat, deviates, info, means)
        !> Generate multivariate normal deviates with zero mean and specified covariance.
        !! Uses cov_to_ltrans to get transpose(L), then transforms standard normal deviates.
        !! Output deviates have shape (n_samples, n_vars), with variables in columns.
        !! Optional means sets column means if provided.
        real(kind=dp), intent(in) :: cov_mat(:,:)
        real(kind=dp), intent(out) :: deviates(:,:)  ! (n_samples, n_vars)
        integer, intent(out) :: info
        real(kind=dp), intent(in), optional :: means(:)  ! Column means
        integer :: n_samples, n_vars, i
        real(kind=dp), allocatable :: ltrans(:,:), z(:,:)

        n_samples = size(deviates, 1)
        n_vars = size(cov_mat, 1)
        if (size(cov_mat, 2) /= n_vars .or. size(deviates, 2) /= n_vars) then
            info = 1  ! Invalid dimensions
            return
        end if
        if (present(means) .and. size(means) /= n_vars) then
            info = 9  ! Invalid means size
            return
        end if

        ! Get transpose(L) from cov_to_ltrans
        allocate(ltrans(n_vars, n_vars))
        call cov_to_ltrans(cov_mat, ltrans, info)
        if (info /= 0) then
            deallocate(ltrans)
            return  ! Error already set by cov_to_ltrans
        end if

        ! Generate and transform deviates
        allocate(z(n_samples, n_vars))
        z = random_normal(n_samples, n_vars)
        deviates = matmul(z, ltrans)

        ! Apply means if provided
        if (present(means)) then
            do i = 1, n_vars
                deviates(:,i) = deviates(:,i) + means(i)
            end do
        end if

        deallocate(ltrans, z)
    end subroutine generate_mv_normal

    subroutine cov_to_ltrans(cov_mat, ltrans, info)
        !> Compute the transpose of the Cholesky factor L for a covariance matrix.
        !! Validates the covariance matrix (symmetric, positive semi-definite) and
        !! returns ltrans = transpose(L) where cov_mat = L * L^T. Sets info to signal errors.
        real(kind=dp), intent(in) :: cov_mat(:,:)
        real(kind=dp), intent(out) :: ltrans(:,:)  ! transpose(L), (n_vars, n_vars)
        integer, intent(out) :: info
        integer :: n_vars, i, j
        real(kind=dp), allocatable :: L(:,:), temp_cov(:,:), eigenvalues(:), LLT(:,:)

        n_vars = size(cov_mat, 1)
        if (size(cov_mat, 2) /= n_vars .or. size(ltrans, 1) /= n_vars .or. &
            size(ltrans, 2) /= n_vars) then
            info = 1  ! Invalid dimensions
            return
        end if

        info = 0
        do i = 1, n_vars
            do j = 1, i - 1
                if (abs(cov_mat(i,j) - cov_mat(j,i)) > TOL) then
                    info = 2  ! Not symmetric
                    return
                end if
            end do
        end do
        do i = 1, n_vars
            if (cov_mat(i,i) <= 0.0_dp) then
                info = 3  ! Non-positive diagonal
                return
            end if
        end do

        allocate(temp_cov(n_vars, n_vars), eigenvalues(n_vars))
        temp_cov = cov_mat
        call eigenvalues_jacobi(temp_cov, eigenvalues, 1000, info)
        if (info /= 0) then
            info = 5  ! Eigenvalue computation failed
            deallocate(temp_cov, eigenvalues)
            return
        end if
        do i = 1, n_vars
            if (eigenvalues(i) < -TOL) then
                info = 6  ! Not positive semi-definite
                deallocate(temp_cov, eigenvalues)
                return
            end if
        end do
        deallocate(eigenvalues, temp_cov)

        allocate(L(n_vars, n_vars))
        call cholesky(cov_mat, L, info)
        if (info /= 0) then
            info = 7  ! Cholesky decomposition failed
            deallocate(L)
            return
        end if

        allocate(LLT(n_vars, n_vars))
        LLT = matmul(L, transpose(L))
        do i = 1, n_vars
            do j = 1, n_vars
                if (abs(LLT(i,j) - cov_mat(i,j)) > 1.0e-6_dp) then
                    print *, "Cholesky error at (", i, ",", j, "): L*L^T =", LLT(i,j), &
                             " vs cov_mat =", cov_mat(i,j)
                    info = 8  ! Cholesky verification failed
                    deallocate(L, LLT)
                    return
                end if
            end do
        end do
        deallocate(LLT)

        ltrans = transpose(L)
        deallocate(L)
    end subroutine cov_to_ltrans

subroutine check_cov_info(info)
integer, intent(in) :: info
    if (info /= 0) then
        print *, "Error: INFO =", info
        select case (info)
            case (1) ; print *, "Invalid dimensions"
            case (2) ; print *, "Matrix not symmetric"
            case (3) ; print *, "Diagonal not 1"
            case (4) ; print *, "Off-diagonal out of [-1, 1]"
            case (5) ; print *, "Eigenvalue computation failed"
            case (6) ; print *, "Not positive semi-definite"
            case (7) ; print *, "Cholesky decomposition failed"
        end select
        error stop
    end if
end subroutine check_cov_info

end module mv_normal_mod
