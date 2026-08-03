program test_pzgemv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zgemv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzgemv
    implicit none

    integer, parameter :: ms(*) = [24, 64, 100]
    integer, parameter :: ns(*) = [32, 48, 120]
    character(len=1), parameter :: transes(*) = ['N', 'T', 'C']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, it, m, n
    integer :: lenx, leny
    integer :: desca(9), descx(9), descy(9)
    character(len=1) :: trans
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), y_got(:), y_ref(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzgemv', target_name, my_rank)

    alpha = cmplx(0.6_ep, 0.2_ep, ep); beta = cmplx(0.3_ep, -0.1_ep, ep)
    do it = 1, size(transes)
        trans = transes(it)
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (trans == 'N') then
                lenx = n; leny = m
            else
                lenx = m; leny = n
            end if

            call gen_distrib_matrix_z(m, n, mb, nb, A_loc, A_glob, &
                                      seed = 4101 + 23 * i + 113 * it)
            call gen_distrib_vector_z(lenx, nb, x_loc, x_glob, &
                                      seed = 4111 + 23 * i + 113 * it)
            call gen_distrib_vector_z(leny, mb, y_loc, y_glob, &
                                      seed = 4121 + 23 * i + 113 * it)

            call local_desc(desca, m, n, mb, nb)
            call local_desc(descx, lenx, 1, nb, 1)
            call local_desc(descy, leny, 1, mb, 1)

            call target_pzgemv(trans, m, n, alpha, A_loc, 1, 1, desca, &
                               x_loc, 1, 1, descx, 1, beta, &
                               y_loc, 1, 1, descy, 1)
            call gather_vector_z(leny, mb, y_loc, y_got)

            if (my_rank == 0) then
                allocate(y_ref(leny))
                y_ref = y_glob
                call zgemv(trans, m, n, alpha, A_glob, m, x_glob, 1, beta, y_ref, 1)
                err = max_rel_err_vec_z(y_got, y_ref)
                tol = 32.0_ep * 8.0_ep * real(max(m, n), ep) * target_eps
                write(label, '(a,a,a,i0,a,i0)') 'trans=', trans, ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(y_ref, y_got)
            end if
            deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgemv
