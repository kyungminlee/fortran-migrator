! Per-target wrapper for the kind16 (qblas) build.
!
! kind16 uses REAL(KIND=16) directly, so wrappers are passthroughs to
! the migrated q-prefix (real) and x-prefix (complex) routines. The
! wrapper interface lets test programs call target_<routine> without
! knowing the target's prefix or precision.

module target_blas
    use prec_kinds, only: ep
    implicit none
    private

    public :: target_name, target_eps
    ! Level 1 — real
    public :: target_dasum, target_daxpy, target_dcopy, target_ddot
    public :: target_dnrm2, target_drot, target_drotg, target_drotm
    public :: target_drotmg, target_dscal, target_dswap, target_idamax
    ! Level 1 — complex
    public :: target_zaxpy, target_zdotc, target_zdotu, target_zscal
    public :: target_dzasum
    ! Level 2 — real
    public :: target_dgemv, target_dgbmv, target_dger, target_dsymv
    public :: target_dspmv, target_dsbmv, target_dtrmv, target_dtbmv
    public :: target_dtpmv, target_dtrsv
    ! Level 2 — complex
    public :: target_zgemv, target_zhemv, target_zgerc
    ! Level 3 — real
    public :: target_dgemm, target_dsymm, target_dsyrk, target_dsyr2k
    public :: target_dtrmm, target_dtrsm
    ! Level 3 — complex
    public :: target_zgemm, target_zhemm, target_zherk, target_ztrsm

    character(len=*), parameter :: target_name = 'kind16'
    real(ep),         parameter :: target_eps  = epsilon(1.0_ep)

    ! Explicit interfaces to the migrated q-prefix (real) and x-prefix
    ! (complex) routines. These map onto the symbols in libqblas.a.
    interface
        ! ── Level 1 — real ───────────────────────────────────────────
        function qdot(n, x, incx, y, incy) result(r)
            import :: ep
            integer,  intent(in) :: n, incx, incy
            real(ep), intent(in) :: x(*), y(*)
            real(ep) :: r
        end function

        function qasum(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function

        function qnrm2(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function

        function iqamax(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            integer :: r
        end function

        subroutine qaxpy(n, alpha, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(in)    :: alpha, x(*)
            real(ep), intent(inout) :: y(*)
        end subroutine

        subroutine qcopy(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)  :: n, incx, incy
            real(ep), intent(in)  :: x(*)
            real(ep), intent(out) :: y(*)
        end subroutine

        subroutine qscal(n, alpha, x, incx)
            import :: ep
            integer,  intent(in)    :: n, incx
            real(ep), intent(in)    :: alpha
            real(ep), intent(inout) :: x(*)
        end subroutine

        subroutine qswap(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
        end subroutine

        subroutine qrot(n, x, incx, y, incy, c, s)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: c, s
        end subroutine

        subroutine qrotg(a, b, c, s)
            import :: ep
            real(ep), intent(inout) :: a, b
            real(ep), intent(out)   :: c, s
        end subroutine

        subroutine qrotm(n, x, incx, y, incy, param)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: param(5)
        end subroutine

        subroutine qrotmg(d1, d2, x1, y1, param)
            import :: ep
            real(ep), intent(inout) :: d1, d2, x1
            real(ep), intent(in)    :: y1
            real(ep), intent(out)   :: param(5)
        end subroutine

        ! ── Level 1 — complex ────────────────────────────────────────
        function xdotc(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function

        function xdotu(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function

        function qxasum(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function

        subroutine xaxpy(n, alpha, x, incx, y, incy)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine

        subroutine xscal(n, alpha, x, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — real ───────────────────────────────────────────
        subroutine qgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine

        subroutine qgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine

        subroutine qger(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,  intent(in)    :: m, n, incx, incy, lda
            real(ep), intent(in)    :: alpha, x(*), y(*)
            real(ep), intent(inout) :: A(lda,*)
        end subroutine

        subroutine qsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine

        subroutine qspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, incy
            real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine

        subroutine qsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, k, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine

        subroutine qtrmv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine

        subroutine qtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, k, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine

        subroutine qtpmv(uplo, trans, diag, n, ap, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, incx
            real(ep),  intent(in)    :: ap(*)
            real(ep),  intent(inout) :: x(*)
        end subroutine

        subroutine qtrsv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex ────────────────────────────────────────
        subroutine xgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine

        subroutine xhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine

        subroutine xgerc(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,     intent(in)    :: m, n, incx, incy, lda
            complex(ep), intent(in)    :: alpha, x(*), y(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine

        ! ── Level 3 — real ───────────────────────────────────────────
        subroutine qgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine

        subroutine qsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine

        subroutine qsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine

        subroutine qsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine

        subroutine qtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine

        subroutine qtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine

        ! ── Level 3 — complex ────────────────────────────────────────
        subroutine xgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine

        subroutine xhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: side, uplo
            integer,     intent(in)    :: m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine

        subroutine xherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldc
            real(ep),    intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine

        subroutine xtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character,   intent(in)    :: side, uplo, transa, diag
            integer,     intent(in)    :: m, n, lda, ldb
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
        end subroutine
    end interface

contains

    ! ── Level 1 — real ───────────────────────────────────────────────
    function target_ddot(n, x, incx, y, incy) result(r)
        integer,  intent(in) :: n, incx, incy
        real(ep), intent(in) :: x(*), y(*)
        real(ep) :: r
        r = qdot(n, x, incx, y, incy)
    end function

    function target_dasum(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        r = qasum(n, x, incx)
    end function

    function target_dnrm2(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        r = qnrm2(n, x, incx)
    end function

    function target_idamax(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        integer :: r
        r = iqamax(n, x, incx)
    end function

    subroutine target_daxpy(n, alpha, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        call qaxpy(n, alpha, x, incx, y, incy)
    end subroutine

    subroutine target_dcopy(n, x, incx, y, incy)
        integer,  intent(in)  :: n, incx, incy
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        call qcopy(n, x, incx, y, incy)
    end subroutine

    subroutine target_dscal(n, alpha, x, incx)
        integer,  intent(in)    :: n, incx
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        call qscal(n, alpha, x, incx)
    end subroutine

    subroutine target_dswap(n, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        call qswap(n, x, incx, y, incy)
    end subroutine

    subroutine target_drot(n, x, incx, y, incy, c, s)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: c, s
        call qrot(n, x, incx, y, incy, c, s)
    end subroutine

    subroutine target_drotg(a, b, c, s)
        real(ep), intent(inout) :: a, b
        real(ep), intent(out)   :: c, s
        call qrotg(a, b, c, s)
    end subroutine

    subroutine target_drotm(n, x, incx, y, incy, param)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: param(5)
        call qrotm(n, x, incx, y, incy, param)
    end subroutine

    subroutine target_drotmg(d1, d2, x1, y1, param)
        real(ep), intent(inout) :: d1, d2, x1
        real(ep), intent(in)    :: y1
        real(ep), intent(out)   :: param(5)
        call qrotmg(d1, d2, x1, y1, param)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    function target_zdotc(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = xdotc(n, x, incx, y, incy)
    end function

    function target_zdotu(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = xdotu(n, x, incx, y, incy)
    end function

    function target_dzasum(n, x, incx) result(r)
        integer,     intent(in) :: n, incx
        complex(ep), intent(in) :: x(*)
        real(ep) :: r
        r = qxasum(n, x, incx)
    end function

    subroutine target_zaxpy(n, alpha, x, incx, y, incy)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        call xaxpy(n, alpha, x, incx, y, incy)
    end subroutine

    subroutine target_zscal(n, alpha, x, incx)
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: x(*)
        call xscal(n, alpha, x, incx)
    end subroutine

    ! ── Level 2 — real ───────────────────────────────────────────────
    subroutine target_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call qgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_dgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call qgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_dger(m, n, alpha, x, incx, y, incy, A, lda)
        integer,  intent(in)    :: m, n, incx, incy, lda
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(lda,*)
        call qger(m, n, alpha, x, incx, y, incy, A, lda)
    end subroutine

    subroutine target_dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call qsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, incx, incy
        real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
        real(ep),  intent(inout) :: y(*)
        call qspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    end subroutine

    subroutine target_dsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, k, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call qsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_dtrmv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call qtrmv(uplo, trans, diag, n, A, lda, x, incx)
    end subroutine

    subroutine target_dtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, k, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call qtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
    end subroutine

    subroutine target_dtpmv(uplo, trans, diag, n, ap, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, incx
        real(ep),  intent(in)    :: ap(*)
        real(ep),  intent(inout) :: x(*)
        call qtpmv(uplo, trans, diag, n, ap, x, incx)
    end subroutine

    subroutine target_dtrsv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call qtrsv(uplo, trans, diag, n, A, lda, x, incx)
    end subroutine

    ! ── Level 2 — complex ────────────────────────────────────────────
    subroutine target_zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call xgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call xhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine target_zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        integer,     intent(in)    :: m, n, incx, incy, lda
        complex(ep), intent(in)    :: alpha, x(*), y(*)
        complex(ep), intent(inout) :: A(lda,*)
        call xgerc(m, n, alpha, x, incx, y, incy, A, lda)
    end subroutine

    ! ── Level 3 — real ───────────────────────────────────────────────
    subroutine target_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call qgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine

    subroutine target_dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call qsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine

    subroutine target_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call qsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
    end subroutine

    subroutine target_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call qsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine

    subroutine target_dtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        call qtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine

    subroutine target_dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        call qtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine

    ! ── Level 3 — complex ────────────────────────────────────────────
    subroutine target_zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call xgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine

    subroutine target_zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: side, uplo
        integer,     intent(in)    :: m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call xhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine

    subroutine target_zherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, lda, ldc
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call xherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
    end subroutine

    subroutine target_ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character,   intent(in)    :: side, uplo, transa, diag
        integer,     intent(in)    :: m, n, lda, ldb
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        call xtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine

end module target_blas
