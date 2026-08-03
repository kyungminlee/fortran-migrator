program test_pdsyr2k
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dsyr2k
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdsyr2k
    implicit none

    integer, parameter :: ns(*) = [32, 80, 128]
    integer, parameter :: ks(*) = [24, 48, 100]
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T']
    integer, parameter :: mb = 8
    integer :: i, ic, n, k
    integer :: ar, ac
    integer :: desca(9), descb(9), descc(9)
    character(len=1) :: uplo, trans
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), &
                             C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdsyr2k', target_name, my_rank)

    alpha = 0.5_ep; beta = 0.25_ep
    do ic = 1, size(uplos)
        uplo = uplos(ic); trans = transes(ic)
        do i = 1, size(ns)
            n = ns(i); k = ks(i)
            if (trans == 'N') then
                ar = n; ac = k
            else
                ar = k; ac = n
            end if

            call gen_distrib_matrix(ar, ac, mb, mb, A_loc, A_glob, &
                                    seed = 14001 + 41 * i + 211 * ic)
            call gen_distrib_matrix(ar, ac, mb, mb, B_loc, B_glob, &
                                    seed = 14011 + 41 * i + 211 * ic)
            call gen_distrib_matrix(n, n, mb, mb, C_loc, C0, &
                                    seed = 14021 + 41 * i + 211 * ic)

            call local_desc(desca, ar, ac, mb, mb)
            call local_desc(descb, ar, ac, mb, mb)
            call local_desc(descc, n, n, mb, mb)

            call target_pdsyr2k(uplo, trans, n, k, alpha, &
                                A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                                beta, C_loc, 1, 1, descc)
            call gather_matrix(n, n, mb, mb, C_loc, C_got)

            if (my_rank == 0) then
                allocate(C_ref(n, n))
                C_ref = C0
                call dsyr2k(uplo, trans, n, k, alpha, A_glob, ar, B_glob, ar, &
                            beta, C_ref, n)
                err = max_rel_err_mat(C_got, C_ref)
                tol = 32.0_ep * 2.0_ep * real(k, ep) * target_eps
                write(label, '(a,a,a,a,a,i0,a,i0)') &
                    'uplo=', uplo, ',trans=', trans, ',n=', n, ',k=', k
                call report_case(trim(label), err, tol)
                deallocate(C_ref, C_got)
            end if
            deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsyr2k
