! pdasymv: row 1-norm of |alpha|*|A_sym|*|X|, plus |beta*Y|.
! A is symmetric (only the UPLO triangle is referenced); the
! reference must reflect the symmetric entries when summing.
program test_pdasymv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix, gen_distrib_vector, &
                             gather_vector
    use target_pblas,  only: target_name, target_eps, target_pdasymv
    implicit none

    integer, parameter :: ns(*) = [40, 80]
    character(len=1), parameter :: uplos(*) = ['U', 'L']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, iu, n, ii, jj
    integer :: desca(9), descx(9), descy(9)
    character(len=1) :: uplo
    real(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    real(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol, acc, aij
    character(len=48) :: label

    call grid_init()
    call report_init('pdasymv', target_name, my_rank)

    alpha = 0.6_ep; beta = 0.3_ep
    do iu = 1, size(uplos)
        uplo = uplos(iu)
        do i = 1, size(ns)
            n = ns(i)

            call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, &
                                    seed = 6201 + 11 * i + 113 * iu)
            call gen_distrib_vector(n, nb, x_loc, x_glob, &
                                    seed = 6221 + 11 * i + 113 * iu)
            call gen_distrib_vector(n, mb, y_loc, y_glob, &
                                    seed = 6241 + 11 * i + 113 * iu)

            call local_desc(desca, n, n, mb, nb)
            call local_desc(descx, n, 1, nb, 1)
            call local_desc(descy, n, 1, mb, 1)

            call target_pdasymv(uplo, n, alpha, A_loc, 1, 1, desca, &
                                x_loc, 1, 1, descx, 1, beta, &
                                y_loc, 1, 1, descy, 1)
            call gather_vector(n, mb, y_loc, y_got)

            if (my_rank == 0) then
                allocate(y_ref(n))
                ! Reference: A is symmetric with only the UPLO triangle
                ! referenced; mirror via aij = A(max(i,j), min(i,j)) for L,
                ! aij = A(min(i,j), max(i,j)) for U.
                do ii = 1, n
                    acc = 0.0_ep
                    do jj = 1, n
                        if (uplo == 'U') then
                            if (ii <= jj) then
                                aij = A_glob(ii, jj)
                            else
                                aij = A_glob(jj, ii)
                            end if
                        else
                            if (ii >= jj) then
                                aij = A_glob(ii, jj)
                            else
                                aij = A_glob(jj, ii)
                            end if
                        end if
                        acc = acc + abs(aij) * abs(x_glob(jj))
                    end do
                    y_ref(ii) = abs(alpha) * acc + abs(beta * y_glob(ii))
                end do
                err = max_rel_err_vec(y_got, y_ref)
                tol = 32.0_ep * 2.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,i0)') 'uplo=', uplo, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(y_ref, y_got)
            end if
            deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pdasymv
