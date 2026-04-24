! Per-target PBLAS wrapper for the kind16 (qpblas) build.
!
! kind16 uses REAL(KIND=16) directly, so the wrappers are passthroughs
! to the migrated q-prefix (real) and x-prefix (complex) PBLAS
! routines. ep == KIND(1.0_ep) == 16 already, so no array conversion
! is needed — test code works in the same type the target library
! computes in.
!
! Migrated PBLAS entry points follow the scalapack-style prefix rule
! (single-character substitution at the classifier's type slot):
!   pdgemm → pqgemm, pdgemv → pqgemv, pddot → pqdot, ...
!   pzgemm → pxgemm, pzaxpy → pxaxpy, pzdotc → pxdotc, pzherk → pxherk
!
! Each C entry point exposes a Fortran-compatible symbol with a
! trailing underscore (pqgemm_), which gfortran resolves from
! `call pqgemm(...)` automatically. Descriptor arguments (DESCA, ...)
! are 9-element INTEGER arrays and pass through unchanged.

module target_pblas
    use prec_kinds, only: ep
    implicit none
    private

    public :: target_name, target_eps
    ! Level 1 — real
    public :: target_pddot, target_pdnrm2, target_pdasum
    public :: target_pdscal, target_pdaxpy, target_pdcopy
    ! Level 1 — complex
    public :: target_pzdotc, target_pzaxpy
    ! Level 2 — real
    public :: target_pdgemv, target_pdger, target_pdsymv, target_pdtrsv
    ! Level 2 — complex
    public :: target_pzgemv
    ! Level 3 — real
    public :: target_pdgemm, target_pdsymm, target_pdsyrk
    public :: target_pdtrmm, target_pdtrsm
    ! Level 3 — complex
    public :: target_pzgemm, target_pzherk

    character(len=*), parameter :: target_name = 'kind16'
    real(ep),         parameter :: target_eps  = epsilon(1.0_ep)

    ! Explicit interfaces to the migrated q/x-prefix PBLAS routines.
    ! Arrays are assumed-size (*) so the wrapper doesn't care about
    ! local leading dimension — that's encoded in the descriptor.
    interface
        ! ── Level 1 — real (Q) ───────────────────────────────────────
        subroutine pqdot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)    :: descx(9), descy(9)
            real(ep), intent(in)    :: x(*), y(*)
            real(ep), intent(out)   :: dot
        end subroutine
        subroutine pqnrm2(n, norm2, x, ix, jx, descx, incx)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx
            integer,  intent(in)    :: descx(9)
            real(ep), intent(in)    :: x(*)
            real(ep), intent(out)   :: norm2
        end subroutine
        subroutine pqasum(n, asum, x, ix, jx, descx, incx)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx
            integer,  intent(in)    :: descx(9)
            real(ep), intent(in)    :: x(*)
            real(ep), intent(out)   :: asum
        end subroutine
        subroutine pqscal(n, alpha, x, ix, jx, descx, incx)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx
            integer,  intent(in)    :: descx(9)
            real(ep), intent(in)    :: alpha
            real(ep), intent(inout) :: x(*)
        end subroutine
        subroutine pqaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)    :: descx(9), descy(9)
            real(ep), intent(in)    :: alpha, x(*)
            real(ep), intent(inout) :: y(*)
        end subroutine
        subroutine pqcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: ep
            integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)    :: descx(9), descy(9)
            real(ep), intent(in)    :: x(*)
            real(ep), intent(out)   :: y(*)
        end subroutine

        ! ── Level 1 — complex (X) ────────────────────────────────────
        subroutine pxdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: ep
            integer,     intent(in)  :: n, ix, jx, incx, iy, jy, incy
            integer,     intent(in)  :: descx(9), descy(9)
            complex(ep), intent(in)  :: x(*), y(*)
            complex(ep), intent(out) :: dot
        end subroutine
        subroutine pxaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: ep
            integer,     intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,     intent(in)    :: descx(9), descy(9)
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine

        ! ── Level 2 — real (Q) ───────────────────────────────────────
        subroutine pqgemv(trans, m, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in)    :: desca(9), descx(9), descy(9)
            real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine pqger(m, n, alpha, x, ix, jx, descx, incx, &
                         y, iy, jy, descy, incy, A, ia, ja, desca)
            import :: ep
            integer,  intent(in)    :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
            integer,  intent(in)    :: desca(9), descx(9), descy(9)
            real(ep), intent(in)    :: alpha, x(*), y(*)
            real(ep), intent(inout) :: A(*)
        end subroutine
        subroutine pqsymv(uplo, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in)    :: desca(9), descx(9), descy(9)
            real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine pqtrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                          x, ix, jx, descx, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, ia, ja, ix, jx, incx
            integer,   intent(in)    :: desca(9), descx(9)
            real(ep),  intent(in)    :: A(*)
            real(ep),  intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex (X) ────────────────────────────────────
        subroutine pxgemv(trans, m, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,     intent(in)    :: desca(9), descx(9), descy(9)
            complex(ep), intent(in)    :: alpha, beta, A(*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine

        ! ── Level 3 — real (Q) ───────────────────────────────────────
        subroutine pqgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: ep
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,   intent(in)    :: desca(9), descb(9), descc(9)
            real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
            real(ep),  intent(inout) :: C(*)
        end subroutine
        subroutine pqsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: ep
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, ia, ja, ib, jb, ic, jc
            integer,   intent(in)    :: desca(9), descb(9), descc(9)
            real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
            real(ep),  intent(inout) :: C(*)
        end subroutine
        subroutine pqsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                          beta, C, ic, jc, descc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, ia, ja, ic, jc
            integer,   intent(in)    :: desca(9), descc(9)
            real(ep),  intent(in)    :: alpha, beta, A(*)
            real(ep),  intent(inout) :: C(*)
        end subroutine
        subroutine pqtrmm(side, uplo, trans, diag, m, n, alpha, &
                          A, ia, ja, desca, B, ib, jb, descb)
            import :: ep
            character, intent(in)    :: side, uplo, trans, diag
            integer,   intent(in)    :: m, n, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            real(ep),  intent(in)    :: alpha, A(*)
            real(ep),  intent(inout) :: B(*)
        end subroutine
        subroutine pqtrsm(side, uplo, trans, diag, m, n, alpha, &
                          A, ia, ja, desca, B, ib, jb, descb)
            import :: ep
            character, intent(in)    :: side, uplo, trans, diag
            integer,   intent(in)    :: m, n, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            real(ep),  intent(in)    :: alpha, A(*)
            real(ep),  intent(inout) :: B(*)
        end subroutine

        ! ── Level 3 — complex (X) ────────────────────────────────────
        subroutine pxgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: ep
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,     intent(in)    :: desca(9), descb(9), descc(9)
            complex(ep), intent(in)    :: alpha, beta, A(*), B(*)
            complex(ep), intent(inout) :: C(*)
        end subroutine
        subroutine pxherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                          beta, C, ic, jc, descc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, ia, ja, ic, jc
            integer,     intent(in)    :: desca(9), descc(9)
            real(ep),    intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(*)
            complex(ep), intent(inout) :: C(*)
        end subroutine
    end interface

contains

    ! ── Level 1 — real ───────────────────────────────────────────────
    subroutine target_pddot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*), y(*)
        real(ep), intent(out) :: dot
        call pqdot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
    end subroutine

    subroutine target_pdnrm2(n, norm2, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: norm2
        call pqnrm2(n, norm2, x, ix, jx, descx, incx)
    end subroutine

    subroutine target_pdasum(n, asum, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: asum
        call pqasum(n, asum, x, ix, jx, descx, incx)
    end subroutine

    subroutine target_pdscal(n, alpha, x, ix, jx, descx, incx)
        integer,  intent(in)    :: n, ix, jx, incx
        integer,  intent(in)    :: descx(9)
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        call pqscal(n, alpha, x, ix, jx, descx, incx)
    end subroutine

    subroutine target_pdaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)    :: descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        call pqaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
    end subroutine

    subroutine target_pdcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        call pqcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    subroutine target_pzdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)  :: descx(9), descy(9)
        complex(ep), intent(in)  :: x(*), y(*)
        complex(ep), intent(out) :: dot
        call pxdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
    end subroutine

    subroutine target_pzaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)    :: descx(9), descy(9)
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        call pxaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
    end subroutine

    ! ── Level 2 — real ───────────────────────────────────────────────
    subroutine target_pdgemv(trans, m, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,   intent(in)    :: desca(9), descx(9), descy(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
        real(ep),  intent(inout) :: y(*)
        call pqgemv(trans, m, n, alpha, A, ia, ja, desca, &
                    x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
    end subroutine

    subroutine target_pdger(m, n, alpha, x, ix, jx, descx, incx, &
                            y, iy, jy, descy, incy, A, ia, ja, desca)
        integer,  intent(in)    :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
        integer,  intent(in)    :: desca(9), descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(*)
        call pqger(m, n, alpha, x, ix, jx, descx, incx, &
                   y, iy, jy, descy, incy, A, ia, ja, desca)
    end subroutine

    subroutine target_pdsymv(uplo, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,   intent(in)    :: desca(9), descx(9), descy(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
        real(ep),  intent(inout) :: y(*)
        call pqsymv(uplo, n, alpha, A, ia, ja, desca, &
                    x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
    end subroutine

    subroutine target_pdtrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                             x, ix, jx, descx, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx
        integer,   intent(in)    :: desca(9), descx(9)
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: x(*)
        call pqtrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                    x, ix, jx, descx, incx)
    end subroutine

    ! ── Level 2 — complex ────────────────────────────────────────────
    subroutine target_pzgemv(trans, m, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,     intent(in)    :: desca(9), descx(9), descy(9)
        complex(ep), intent(in)    :: alpha, beta, A(*), x(*)
        complex(ep), intent(inout) :: y(*)
        call pxgemv(trans, m, n, alpha, A, ia, ja, desca, &
                    x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
    end subroutine

    ! ── Level 3 — real ───────────────────────────────────────────────
    subroutine target_pdgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
        integer,   intent(in)    :: desca(9), descb(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
        real(ep),  intent(inout) :: C(*)
        call pqgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                    B, ib, jb, descb, beta, C, ic, jc, descc)
    end subroutine

    subroutine target_pdsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, ia, ja, ib, jb, ic, jc
        integer,   intent(in)    :: desca(9), descb(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
        real(ep),  intent(inout) :: C(*)
        call pqsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                    B, ib, jb, descb, beta, C, ic, jc, descc)
    end subroutine

    subroutine target_pdsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                             beta, C, ic, jc, descc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, ia, ja, ic, jc
        integer,   intent(in)    :: desca(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*)
        real(ep),  intent(inout) :: C(*)
        call pqsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                    beta, C, ic, jc, descc)
    end subroutine

    subroutine target_pdtrmm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        call pqtrmm(side, uplo, trans, diag, m, n, alpha, &
                    A, ia, ja, desca, B, ib, jb, descb)
    end subroutine

    subroutine target_pdtrsm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        call pqtrsm(side, uplo, trans, diag, m, n, alpha, &
                    A, ia, ja, desca, B, ib, jb, descb)
    end subroutine

    ! ── Level 3 — complex ────────────────────────────────────────────
    subroutine target_pzgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
        integer,     intent(in)    :: desca(9), descb(9), descc(9)
        complex(ep), intent(in)    :: alpha, beta, A(*), B(*)
        complex(ep), intent(inout) :: C(*)
        call pxgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                    B, ib, jb, descb, beta, C, ic, jc, descc)
    end subroutine

    subroutine target_pzherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                             beta, C, ic, jc, descc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, ia, ja, ic, jc
        integer,     intent(in)    :: desca(9), descc(9)
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(*)
        complex(ep), intent(inout) :: C(*)
        call pxherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                    beta, C, ic, jc, descc)
    end subroutine

end module target_pblas
