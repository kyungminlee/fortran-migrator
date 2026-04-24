! Per-target wrapper for the kind10 (eblas) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers cast quad
! inputs down to REAL(KIND=10) for the call into the migrated e/y
! prefix routines, then cast results back up to quad. The cast-down
! is the inherent precision floor of the target — comparison against
! the quad reference will then measure how closely the kind10 result
! reproduces the quad value (~17 decimal digits agreement expected).
!
! Migrated routine prefixes for kind10:
!   D → E   (egemm, edot, eaxpy, …)
!   Z → Y   (ygemm, ydotc, yaxpy, …)
!   I → I   (ieamax, iyamax)
!   DZ → EY (eyasum)

module target_blas
    use prec_kinds, only: ep
    implicit none
    private

    public :: target_name, target_eps
    public :: target_dasum, target_daxpy, target_dcopy, target_ddot
    public :: target_dnrm2, target_drot, target_drotg, target_drotm
    public :: target_drotmg, target_dscal, target_dswap, target_idamax
    public :: target_zaxpy, target_zdotc, target_zdotu, target_zscal
    public :: target_dzasum
    public :: target_dgemv, target_dgbmv, target_dger, target_dsymv
    public :: target_dspmv, target_dsbmv, target_dtrmv, target_dtbmv
    public :: target_dtpmv, target_dtrsv
    public :: target_zgemv, target_zhemv, target_zgerc
    public :: target_dgemm, target_dsymm, target_dsyrk, target_dsyr2k
    public :: target_dtrmm, target_dtrsm
    public :: target_zgemm, target_zhemm, target_zherk, target_ztrsm

    integer,          parameter :: tk = 10  ! target kind
    character(len=*), parameter :: target_name = 'kind10'
    real(ep),         parameter :: target_eps  = real(epsilon(1.0_10), ep)

    interface
        ! ── e-prefix (D family) ──────────────────────────────────────
        function edot(n, x, incx, y, incy) result(r)
            integer,  intent(in)  :: n, incx, incy
            real(10), intent(in)  :: x(*), y(*)
            real(10) :: r
        end function
        function easum(n, x, incx) result(r)
            integer,  intent(in) :: n, incx
            real(10), intent(in) :: x(*)
            real(10) :: r
        end function
        function enrm2(n, x, incx) result(r)
            integer,  intent(in) :: n, incx
            real(10), intent(in) :: x(*)
            real(10) :: r
        end function
        function ieamax(n, x, incx) result(r)
            integer,  intent(in) :: n, incx
            real(10), intent(in) :: x(*)
            integer :: r
        end function
        subroutine eaxpy(n, alpha, x, incx, y, incy)
            integer,  intent(in)    :: n, incx, incy
            real(10), intent(in)    :: alpha, x(*)
            real(10), intent(inout) :: y(*)
        end subroutine
        subroutine ecopy(n, x, incx, y, incy)
            integer,  intent(in)  :: n, incx, incy
            real(10), intent(in)  :: x(*)
            real(10), intent(out) :: y(*)
        end subroutine
        subroutine escal(n, alpha, x, incx)
            integer,  intent(in)    :: n, incx
            real(10), intent(in)    :: alpha
            real(10), intent(inout) :: x(*)
        end subroutine
        subroutine eswap(n, x, incx, y, incy)
            integer,  intent(in)    :: n, incx, incy
            real(10), intent(inout) :: x(*), y(*)
        end subroutine
        subroutine erot(n, x, incx, y, incy, c, s)
            integer,  intent(in)    :: n, incx, incy
            real(10), intent(inout) :: x(*), y(*)
            real(10), intent(in)    :: c, s
        end subroutine
        subroutine erotg(a, b, c, s)
            real(10), intent(inout) :: a, b
            real(10), intent(out)   :: c, s
        end subroutine
        subroutine erotm(n, x, incx, y, incy, param)
            integer,  intent(in)    :: n, incx, incy
            real(10), intent(inout) :: x(*), y(*)
            real(10), intent(in)    :: param(5)
        end subroutine
        subroutine erotmg(d1, d2, x1, y1, param)
            real(10), intent(inout) :: d1, d2, x1
            real(10), intent(in)    :: y1
            real(10), intent(out)   :: param(5)
        end subroutine

        ! ── y-prefix (Z family) ──────────────────────────────────────
        function ydotc(n, x, incx, y, incy) result(r)
            integer,     intent(in) :: n, incx, incy
            complex(10), intent(in) :: x(*), y(*)
            complex(10) :: r
        end function
        function ydotu(n, x, incx, y, incy) result(r)
            integer,     intent(in) :: n, incx, incy
            complex(10), intent(in) :: x(*), y(*)
            complex(10) :: r
        end function
        function eyasum(n, x, incx) result(r)
            integer,     intent(in) :: n, incx
            complex(10), intent(in) :: x(*)
            real(10) :: r
        end function
        subroutine yaxpy(n, alpha, x, incx, y, incy)
            integer,     intent(in)    :: n, incx, incy
            complex(10), intent(in)    :: alpha, x(*)
            complex(10), intent(inout) :: y(*)
        end subroutine
        subroutine yscal(n, alpha, x, incx)
            integer,     intent(in)    :: n, incx
            complex(10), intent(in)    :: alpha
            complex(10), intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — real (E) ───────────────────────────────────────
        subroutine egemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, lda, incx, incy
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine egbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine eger(m, n, alpha, x, incx, y, incy, A, lda)
            integer,  intent(in)    :: m, n, incx, incy, lda
            real(10), intent(in)    :: alpha, x(*), y(*)
            real(10), intent(inout) :: A(lda,*)
        end subroutine
        subroutine esymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, incx, incy
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine espmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, incy
            real(10),  intent(in)    :: alpha, beta, ap(*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine esbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, k, lda, incx, incy
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine etrmv(uplo, trans, diag, n, A, lda, x, incx)
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: x(*)
        end subroutine
        subroutine etbmv(uplo, trans, diag, n, k, A, lda, x, incx)
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, k, lda, incx
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: x(*)
        end subroutine
        subroutine etpmv(uplo, trans, diag, n, ap, x, incx)
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, incx
            real(10),  intent(in)    :: ap(*)
            real(10),  intent(inout) :: x(*)
        end subroutine
        subroutine etrsv(uplo, trans, diag, n, A, lda, x, incx)
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex (Y) ────────────────────────────────────
        subroutine ygemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, lda, incx, incy
            complex(10), intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(lda,*), x(*)
            complex(10), intent(inout) :: y(*)
        end subroutine
        subroutine yhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx, incy
            complex(10), intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(lda,*), x(*)
            complex(10), intent(inout) :: y(*)
        end subroutine
        subroutine ygerc(m, n, alpha, x, incx, y, incy, A, lda)
            integer,     intent(in)    :: m, n, incx, incy, lda
            complex(10), intent(in)    :: alpha, x(*), y(*)
            complex(10), intent(inout) :: A(lda,*)
        end subroutine

        ! ── Level 3 — real (E) ───────────────────────────────────────
        subroutine egemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, lda, ldb, ldc
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), B(ldb,*)
            real(10),  intent(inout) :: C(ldc,*)
        end subroutine
        subroutine esymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, lda, ldb, ldc
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), B(ldb,*)
            real(10),  intent(inout) :: C(ldc,*)
        end subroutine
        subroutine esyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldc
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: C(ldc,*)
        end subroutine
        subroutine esyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldb, ldc
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(in)    :: A(lda,*), B(ldb,*)
            real(10),  intent(inout) :: C(ldc,*)
        end subroutine
        subroutine etrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(10),  intent(in)    :: alpha
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: B(ldb,*)
        end subroutine
        subroutine etrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(10),  intent(in)    :: alpha
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: B(ldb,*)
        end subroutine

        ! ── Level 3 — complex (Y) ────────────────────────────────────
        subroutine ygemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, lda, ldb, ldc
            complex(10), intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(lda,*), B(ldb,*)
            complex(10), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine yhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            character,   intent(in)    :: side, uplo
            integer,     intent(in)    :: m, n, lda, ldb, ldc
            complex(10), intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(lda,*), B(ldb,*)
            complex(10), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine yherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldc
            real(10),    intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(lda,*)
            complex(10), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine ytrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            character,   intent(in)    :: side, uplo, transa, diag
            integer,     intent(in)    :: m, n, lda, ldb
            complex(10), intent(in)    :: alpha
            complex(10), intent(in)    :: A(lda,*)
            complex(10), intent(inout) :: B(ldb,*)
        end subroutine
    end interface

contains

    ! ── helpers ──────────────────────────────────────────────────────
    pure function span_real(arr, n, inc) result(m)
        real(ep), intent(in) :: arr(*)
        integer,  intent(in) :: n, inc
        integer :: m
        m = (n - 1) * abs(inc) + 1
        ! `arr` is just to enable an assumed-size pattern in callers.
        if (.false.) m = m + int(arr(1) * 0.0_ep)
    end function

    ! ── Level 1 — real ───────────────────────────────────────────────
    function target_ddot(n, x, incx, y, incy) result(r)
        integer,  intent(in) :: n, incx, incy
        real(ep), intent(in) :: x(*), y(*)
        real(ep) :: r
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        r = real(edot(n, xt, incx, yt, incy), ep)
    end function

    function target_dasum(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        real(tk), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        r = real(easum(n, xt, incx), ep)
    end function

    function target_dnrm2(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        real(tk), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        r = real(enrm2(n, xt, incx), ep)
    end function

    function target_idamax(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        integer :: r
        integer :: nx
        real(tk), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        r = ieamax(n, xt, incx)
    end function

    subroutine target_daxpy(n, alpha, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        real(tk) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        alpha_t = real(alpha, tk)
        call eaxpy(n, alpha_t, xt, incx, yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dcopy(n, x, incx, y, incy)
        integer,  intent(in)  :: n, incx, incy
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        call ecopy(n, xt, incx, yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dscal(n, alpha, x, incx)
        integer,  intent(in)    :: n, incx
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        integer :: nx
        real(tk), allocatable :: xt(:)
        real(tk) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        alpha_t = real(alpha, tk)
        call escal(n, alpha_t, xt, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    subroutine target_dswap(n, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        call eswap(n, xt, incx, yt, incy)
        x(1:nx) = real(xt, ep)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_drot(n, x, incx, y, incy, c, s)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: c, s
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        real(tk) :: c_t, s_t
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        c_t = real(c, tk); s_t = real(s, tk)
        call erot(n, xt, incx, yt, incy, c_t, s_t)
        x(1:nx) = real(xt, ep)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_drotg(a, b, c, s)
        real(ep), intent(inout) :: a, b
        real(ep), intent(out)   :: c, s
        real(tk) :: a_t, b_t, c_t, s_t
        a_t = real(a, tk); b_t = real(b, tk)
        call erotg(a_t, b_t, c_t, s_t)
        a = real(a_t, ep); b = real(b_t, ep)
        c = real(c_t, ep); s = real(s_t, ep)
    end subroutine

    subroutine target_drotm(n, x, incx, y, incy, param)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: param(5)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        real(tk) :: param_t(5)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        param_t = real(param, tk)
        call erotm(n, xt, incx, yt, incy, param_t)
        x(1:nx) = real(xt, ep)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_drotmg(d1, d2, x1, y1, param)
        real(ep), intent(inout) :: d1, d2, x1
        real(ep), intent(in)    :: y1
        real(ep), intent(out)   :: param(5)
        real(tk) :: d1_t, d2_t, x1_t, y1_t, param_t(5)
        d1_t = real(d1, tk); d2_t = real(d2, tk)
        x1_t = real(x1, tk); y1_t = real(y1, tk)
        call erotmg(d1_t, d2_t, x1_t, y1_t, param_t)
        d1 = real(d1_t, ep); d2 = real(d2_t, ep); x1 = real(x1_t, ep)
        param = real(param_t, ep)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    function target_zdotc(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        integer :: nx, ny
        complex(tk), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = cmplx(x(1:nx), kind=tk)
        yt = cmplx(y(1:ny), kind=tk)
        r = cmplx(ydotc(n, xt, incx, yt, incy), kind=ep)
    end function

    function target_zdotu(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        integer :: nx, ny
        complex(tk), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = cmplx(x(1:nx), kind=tk)
        yt = cmplx(y(1:ny), kind=tk)
        r = cmplx(ydotu(n, xt, incx, yt, incy), kind=ep)
    end function

    function target_dzasum(n, x, incx) result(r)
        integer,     intent(in) :: n, incx
        complex(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        complex(tk), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = cmplx(x(1:nx), kind=tk)
        r = real(eyasum(n, xt, incx), ep)
    end function

    subroutine target_zaxpy(n, alpha, x, incx, y, incy)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        complex(tk), allocatable :: xt(:), yt(:)
        complex(tk) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = cmplx(x(1:nx), kind=tk)
        yt = cmplx(y(1:ny), kind=tk)
        alpha_t = cmplx(alpha, kind=tk)
        call yaxpy(n, alpha_t, xt, incx, yt, incy)
        y(1:ny) = cmplx(yt, kind=ep)
    end subroutine

    subroutine target_zscal(n, alpha, x, incx)
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: x(*)
        integer :: nx
        complex(tk), allocatable :: xt(:)
        complex(tk) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = cmplx(x(1:nx), kind=tk)
        alpha_t = cmplx(alpha, kind=tk)
        call yscal(n, alpha_t, xt, incx)
        x(1:nx) = cmplx(xt, kind=ep)
    end subroutine

    ! ── Level 2 — real ───────────────────────────────────────────────
    subroutine target_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny, ar, ac
        real(tk), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
            ac = n
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
            ac = n  ! for 'T'/'C', A is m x n still in storage
        end if
        ar = lda
        allocate(At(ar, ac), xt(nx), yt(ny))
        At(1:ar, 1:ac) = real(A(1:ar, 1:ac), tk)
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        call egemv(trans, m, n, real(alpha, tk), At, lda, xt, incx, &
                   real(beta, tk), yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
        end if
        allocate(At(lda, n), xt(nx), yt(ny))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        call egbmv(trans, m, n, kl, ku, real(alpha, tk), At, lda, &
                   xt, incx, real(beta, tk), yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dger(m, n, alpha, x, incx, y, incy, A, lda)
        integer,  intent(in)    :: m, n, incx, incy, lda
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(lda,*)
        integer :: nx, ny
        real(tk), allocatable :: At(:,:), xt(:), yt(:)
        nx = (m - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk)
        yt = real(y(1:ny), tk)
        call eger(m, n, real(alpha, tk), xt, incx, yt, incy, At, lda)
        A(1:lda, 1:n) = real(At, ep)
    end subroutine

    subroutine target_dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call esymv(uplo, n, real(alpha, tk), At, lda, xt, incx, &
                   real(beta, tk), yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, incx, incy
        real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: aps, nx, ny
        real(tk), allocatable :: apt(:), xt(:), yt(:)
        aps = n * (n + 1) / 2
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(apt(aps), xt(nx), yt(ny))
        apt = real(ap(1:aps), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call espmv(uplo, n, real(alpha, tk), apt, xt, incx, &
                   real(beta, tk), yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, k, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call esbmv(uplo, n, k, real(alpha, tk), At, lda, xt, incx, &
                   real(beta, tk), yt, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_dtrmv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        real(tk), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk)
        call etrmv(uplo, trans, diag, n, At, lda, xt, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    subroutine target_dtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, k, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        real(tk), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk)
        call etbmv(uplo, trans, diag, n, k, At, lda, xt, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    subroutine target_dtpmv(uplo, trans, diag, n, ap, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, incx
        real(ep),  intent(in)    :: ap(*)
        real(ep),  intent(inout) :: x(*)
        integer :: aps, nx
        real(tk), allocatable :: apt(:), xt(:)
        aps = n * (n + 1) / 2
        nx = (n - 1) * abs(incx) + 1
        allocate(apt(aps), xt(nx))
        apt = real(ap(1:aps), tk)
        xt = real(x(1:nx), tk)
        call etpmv(uplo, trans, diag, n, apt, xt, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    subroutine target_dtrsv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        real(tk), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = real(A(1:lda, 1:n), tk)
        xt = real(x(1:nx), tk)
        call etrsv(uplo, trans, diag, n, At, lda, xt, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    ! ── Level 2 — complex ────────────────────────────────────────────
    subroutine target_zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        complex(tk), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
        end if
        allocate(At(lda, n), xt(nx), yt(ny))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call ygemv(trans, m, n, cmplx(alpha, kind=tk), At, lda, &
                   xt, incx, cmplx(beta, kind=tk), yt, incy)
        y(1:ny) = cmplx(yt, kind=ep)
    end subroutine

    subroutine target_zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        complex(tk), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call yhemv(uplo, n, cmplx(alpha, kind=tk), At, lda, xt, incx, &
                   cmplx(beta, kind=tk), yt, incy)
        y(1:ny) = cmplx(yt, kind=ep)
    end subroutine

    subroutine target_zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        integer,     intent(in)    :: m, n, incx, incy, lda
        complex(ep), intent(in)    :: alpha, x(*), y(*)
        complex(ep), intent(inout) :: A(lda,*)
        integer :: nx, ny
        complex(tk), allocatable :: At(:,:), xt(:), yt(:)
        nx = (m - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call ygerc(m, n, cmplx(alpha, kind=tk), xt, incx, yt, incy, At, lda)
        A(1:lda, 1:n) = cmplx(At, kind=ep)
    end subroutine

    ! ── Level 3 — real ───────────────────────────────────────────────
    subroutine target_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an, bn
        real(tk), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (transa == 'N' .or. transa == 'n') then; an = k; else; an = m; end if
        if (transb == 'N' .or. transb == 'n') then; bn = n; else; bn = k; end if
        allocate(At(lda, an), Bt(ldb, bn), Ct(ldc, n))
        At = real(A(1:lda, 1:an), tk)
        Bt = real(B(1:ldb, 1:bn), tk)
        Ct = real(C(1:ldc, 1:n), tk)
        call egemm(transa, transb, m, n, k, real(alpha, tk), At, lda, &
                   Bt, ldb, real(beta, tk), Ct, ldc)
        C(1:ldc, 1:n) = real(Ct, ep)
    end subroutine

    subroutine target_dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        real(tk), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n), Ct(ldc, n))
        At = real(A(1:lda, 1:an), tk)
        Bt = real(B(1:ldb, 1:n), tk)
        Ct = real(C(1:ldc, 1:n), tk)
        call esymm(side, uplo, m, n, real(alpha, tk), At, lda, &
                   Bt, ldb, real(beta, tk), Ct, ldc)
        C(1:ldc, 1:n) = real(Ct, ep)
    end subroutine

    subroutine target_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        real(tk), allocatable :: At(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Ct(ldc, n))
        At = real(A(1:lda, 1:an), tk)
        Ct = real(C(1:ldc, 1:n), tk)
        call esyrk(uplo, trans, n, k, real(alpha, tk), At, lda, &
                   real(beta, tk), Ct, ldc)
        C(1:ldc, 1:n) = real(Ct, ep)
    end subroutine

    subroutine target_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        real(tk), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, an), Ct(ldc, n))
        At = real(A(1:lda, 1:an), tk)
        Bt = real(B(1:ldb, 1:an), tk)
        Ct = real(C(1:ldc, 1:n), tk)
        call esyr2k(uplo, trans, n, k, real(alpha, tk), At, lda, &
                    Bt, ldb, real(beta, tk), Ct, ldc)
        C(1:ldc, 1:n) = real(Ct, ep)
    end subroutine

    subroutine target_dtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer :: an
        real(tk), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = real(A(1:lda, 1:an), tk)
        Bt = real(B(1:ldb, 1:n), tk)
        call etrmm(side, uplo, transa, diag, m, n, real(alpha, tk), &
                   At, lda, Bt, ldb)
        B(1:ldb, 1:n) = real(Bt, ep)
    end subroutine

    subroutine target_dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer :: an
        real(tk), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = real(A(1:lda, 1:an), tk)
        Bt = real(B(1:ldb, 1:n), tk)
        call etrsm(side, uplo, transa, diag, m, n, real(alpha, tk), &
                   At, lda, Bt, ldb)
        B(1:ldb, 1:n) = real(Bt, ep)
    end subroutine

    ! ── Level 3 — complex ────────────────────────────────────────────
    subroutine target_zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an, bn
        complex(tk), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (transa == 'N' .or. transa == 'n') then; an = k; else; an = m; end if
        if (transb == 'N' .or. transb == 'n') then; bn = n; else; bn = k; end if
        allocate(At(lda, an), Bt(ldb, bn), Ct(ldc, n))
        At = cmplx(A(1:lda, 1:an), kind=tk)
        Bt = cmplx(B(1:ldb, 1:bn), kind=tk)
        Ct = cmplx(C(1:ldc, 1:n), kind=tk)
        call ygemm(transa, transb, m, n, k, cmplx(alpha, kind=tk), At, lda, &
                   Bt, ldb, cmplx(beta, kind=tk), Ct, ldc)
        C(1:ldc, 1:n) = cmplx(Ct, kind=ep)
    end subroutine

    subroutine target_zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: side, uplo
        integer,     intent(in)    :: m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an
        complex(tk), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n), Ct(ldc, n))
        At = cmplx(A(1:lda, 1:an), kind=tk)
        Bt = cmplx(B(1:ldb, 1:n), kind=tk)
        Ct = cmplx(C(1:ldc, 1:n), kind=tk)
        call yhemm(side, uplo, m, n, cmplx(alpha, kind=tk), At, lda, &
                   Bt, ldb, cmplx(beta, kind=tk), Ct, ldc)
        C(1:ldc, 1:n) = cmplx(Ct, kind=ep)
    end subroutine

    subroutine target_zherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, lda, ldc
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an
        complex(tk), allocatable :: At(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Ct(ldc, n))
        At = cmplx(A(1:lda, 1:an), kind=tk)
        Ct = cmplx(C(1:ldc, 1:n), kind=tk)
        call yherk(uplo, trans, n, k, real(alpha, tk), At, lda, &
                   real(beta, tk), Ct, ldc)
        C(1:ldc, 1:n) = cmplx(Ct, kind=ep)
    end subroutine

    subroutine target_ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character,   intent(in)    :: side, uplo, transa, diag
        integer,     intent(in)    :: m, n, lda, ldb
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer :: an
        complex(tk), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = cmplx(A(1:lda, 1:an), kind=tk)
        Bt = cmplx(B(1:ldb, 1:n), kind=tk)
        call ytrsm(side, uplo, transa, diag, m, n, cmplx(alpha, kind=tk), &
                   At, lda, Bt, ldb)
        B(1:ldb, 1:n) = cmplx(Bt, kind=ep)
    end subroutine

end module target_blas
