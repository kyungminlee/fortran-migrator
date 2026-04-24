! Per-target wrapper for the multifloats (ddblas) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers split each
! quad value into the high (double approx) + low (double remainder)
! representation that TYPE(float64x2) expects, call the migrated
! dd*/zz* routine, then recombine the two limbs back into quad for
! the comparison. This preserves multifloats' full ~32-decimal-digit
! precision through the wrapper boundary.
!
! Migrated routine prefixes for multifloats:
!   D → DD   (ddgemm, dddot, ddaxpy, …)
!   Z → ZZ   (zzgemm, zzdotc, zzaxpy, …)
!   I → I    (iddamax, izzamax)
!   DZ → DDZZ (ddzzasum)

module target_blas
    use prec_kinds, only: ep, dp
    use multifloats, only: float64x2, complex64x2
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

    character(len=*), parameter :: target_name = 'multifloats'
    ! Double-double effective epsilon ≈ 2^-104 ≈ 4.93e-32.
    real(ep),         parameter :: target_eps  = real(2.0_dp**(-104), ep)

    interface
        ! ── dd-prefix (D family) ─────────────────────────────────────
        function dddot(n, x, incx, y, incy) result(r)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(in) :: x(*), y(*)
            type(float64x2) :: r
        end function
        function ddasum(n, x, incx) result(r)
            import :: float64x2
            integer, intent(in) :: n, incx
            type(float64x2), intent(in) :: x(*)
            type(float64x2) :: r
        end function
        function ddnrm2(n, x, incx) result(r)
            import :: float64x2
            integer, intent(in) :: n, incx
            type(float64x2), intent(in) :: x(*)
            type(float64x2) :: r
        end function
        function iddamax(n, x, incx) result(r)
            import :: float64x2
            integer, intent(in) :: n, incx
            type(float64x2), intent(in) :: x(*)
            integer :: r
        end function
        subroutine ddaxpy(n, alpha, x, incx, y, incy)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(in)    :: alpha, x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddcopy(n, x, incx, y, incy)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(in)  :: x(*)
            type(float64x2), intent(out) :: y(*)
        end subroutine
        subroutine ddscal(n, alpha, x, incx)
            import :: float64x2
            integer, intent(in) :: n, incx
            type(float64x2), intent(in)    :: alpha
            type(float64x2), intent(inout) :: x(*)
        end subroutine
        subroutine ddswap(n, x, incx, y, incy)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(inout) :: x(*), y(*)
        end subroutine
        subroutine ddrot(n, x, incx, y, incy, c, s)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(inout) :: x(*), y(*)
            type(float64x2), intent(in)    :: c, s
        end subroutine
        subroutine ddrotg(a, b, c, s)
            import :: float64x2
            type(float64x2), intent(inout) :: a, b
            type(float64x2), intent(out)   :: c, s
        end subroutine
        subroutine ddrotm(n, x, incx, y, incy, param)
            import :: float64x2
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(inout) :: x(*), y(*)
            type(float64x2), intent(in)    :: param(5)
        end subroutine
        subroutine ddrotmg(d1, d2, x1, y1, param)
            import :: float64x2
            type(float64x2), intent(inout) :: d1, d2, x1
            type(float64x2), intent(in)    :: y1
            type(float64x2), intent(out)   :: param(5)
        end subroutine

        ! ── zz-prefix (Z family) ─────────────────────────────────────
        function zzdotc(n, x, incx, y, incy) result(r)
            import :: complex64x2
            integer, intent(in) :: n, incx, incy
            type(complex64x2), intent(in) :: x(*), y(*)
            type(complex64x2) :: r
        end function
        function zzdotu(n, x, incx, y, incy) result(r)
            import :: complex64x2
            integer, intent(in) :: n, incx, incy
            type(complex64x2), intent(in) :: x(*), y(*)
            type(complex64x2) :: r
        end function
        function ddzzasum(n, x, incx) result(r)
            import :: float64x2, complex64x2
            integer, intent(in) :: n, incx
            type(complex64x2), intent(in) :: x(*)
            type(float64x2) :: r
        end function
        subroutine zzaxpy(n, alpha, x, incx, y, incy)
            import :: complex64x2
            integer, intent(in) :: n, incx, incy
            type(complex64x2), intent(in)    :: alpha, x(*)
            type(complex64x2), intent(inout) :: y(*)
        end subroutine
        subroutine zzscal(n, alpha, x, incx)
            import :: complex64x2
            integer, intent(in) :: n, incx
            type(complex64x2), intent(in)    :: alpha
            type(complex64x2), intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — real (DD) ──────────────────────────────────────
        subroutine ddgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: float64x2
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
            import :: float64x2
            character, intent(in) :: trans
            integer, intent(in) :: m, n, kl, ku, lda, incx, incy
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddger(m, n, alpha, x, incx, y, incy, A, lda)
            import :: float64x2
            integer, intent(in) :: m, n, incx, incy, lda
            type(float64x2), intent(in)    :: alpha, x(*), y(*)
            type(float64x2), intent(inout) :: A(lda,*)
        end subroutine
        subroutine ddsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: float64x2
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda, incx, incy
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: float64x2
            character, intent(in) :: uplo
            integer, intent(in) :: n, incx, incy
            type(float64x2), intent(in)    :: alpha, beta, ap(*), x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
            import :: float64x2
            character, intent(in) :: uplo
            integer, intent(in) :: n, k, lda, incx, incy
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), x(*)
            type(float64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ddtrmv(uplo, trans, diag, n, A, lda, x, incx)
            import :: float64x2
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, lda, incx
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: x(*)
        end subroutine
        subroutine ddtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
            import :: float64x2
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, k, lda, incx
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: x(*)
        end subroutine
        subroutine ddtpmv(uplo, trans, diag, n, ap, x, incx)
            import :: float64x2
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, incx
            type(float64x2), intent(in)    :: ap(*)
            type(float64x2), intent(inout) :: x(*)
        end subroutine
        subroutine ddtrsv(uplo, trans, diag, n, A, lda, x, incx)
            import :: float64x2
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, lda, incx
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex (ZZ) ───────────────────────────────────
        subroutine zzgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: complex64x2
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            type(complex64x2), intent(in)    :: alpha, beta
            type(complex64x2), intent(in)    :: A(lda,*), x(*)
            type(complex64x2), intent(inout) :: y(*)
        end subroutine
        subroutine zzhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: complex64x2
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda, incx, incy
            type(complex64x2), intent(in)    :: alpha, beta
            type(complex64x2), intent(in)    :: A(lda,*), x(*)
            type(complex64x2), intent(inout) :: y(*)
        end subroutine
        subroutine zzgerc(m, n, alpha, x, incx, y, incy, A, lda)
            import :: complex64x2
            integer, intent(in) :: m, n, incx, incy, lda
            type(complex64x2), intent(in)    :: alpha, x(*), y(*)
            type(complex64x2), intent(inout) :: A(lda,*)
        end subroutine

        ! ── Level 3 — real (DD) ──────────────────────────────────────
        subroutine ddgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: float64x2
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), B(ldb,*)
            type(float64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine ddsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: float64x2
            character, intent(in) :: side, uplo
            integer, intent(in) :: m, n, lda, ldb, ldc
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), B(ldb,*)
            type(float64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine ddsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: float64x2
            character, intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldc
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine ddsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: float64x2
            character, intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldb, ldc
            type(float64x2), intent(in)    :: alpha, beta
            type(float64x2), intent(in)    :: A(lda,*), B(ldb,*)
            type(float64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine ddtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: float64x2
            character, intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            type(float64x2), intent(in)    :: alpha
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: B(ldb,*)
        end subroutine
        subroutine ddtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: float64x2
            character, intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            type(float64x2), intent(in)    :: alpha
            type(float64x2), intent(in)    :: A(lda,*)
            type(float64x2), intent(inout) :: B(ldb,*)
        end subroutine

        ! ── Level 3 — complex (ZZ) ───────────────────────────────────
        subroutine zzgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: complex64x2
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            type(complex64x2), intent(in)    :: alpha, beta
            type(complex64x2), intent(in)    :: A(lda,*), B(ldb,*)
            type(complex64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine zzhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: complex64x2
            character, intent(in) :: side, uplo
            integer, intent(in) :: m, n, lda, ldb, ldc
            type(complex64x2), intent(in)    :: alpha, beta
            type(complex64x2), intent(in)    :: A(lda,*), B(ldb,*)
            type(complex64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine zzherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: float64x2, complex64x2
            character, intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldc
            type(float64x2),   intent(in)    :: alpha, beta
            type(complex64x2), intent(in)    :: A(lda,*)
            type(complex64x2), intent(inout) :: C(ldc,*)
        end subroutine
        subroutine zztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: complex64x2
            character, intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            type(complex64x2), intent(in)    :: alpha
            type(complex64x2), intent(in)    :: A(lda,*)
            type(complex64x2), intent(inout) :: B(ldb,*)
        end subroutine
    end interface

contains

    ! ── Conversion helpers ──────────────────────────────────────────
    ! Quad → double-double via the well-known split:
    !   hi = round-to-double(quad)
    !   lo = round-to-double(quad - hi)
    ! `lo` fits in a double because |lo| <= ulp(hi)/2 < 2^-52 * |hi|.

    elemental function q2dd(x) result(r)
        real(ep), intent(in) :: x
        type(float64x2) :: r
        real(dp) :: hi
        hi = real(x, dp)
        r%limbs(1) = hi
        r%limbs(2) = real(x - real(hi, ep), dp)
    end function

    elemental function dd2q(x) result(r)
        type(float64x2), intent(in) :: x
        real(ep) :: r
        r = real(x%limbs(1), ep) + real(x%limbs(2), ep)
    end function

    elemental function q2zz(z) result(r)
        complex(ep), intent(in) :: z
        type(complex64x2) :: r
        r%re = q2dd(real(z, ep))
        r%im = q2dd(aimag(z))
    end function

    elemental function zz2q(z) result(r)
        type(complex64x2), intent(in) :: z
        complex(ep) :: r
        r = cmplx(dd2q(z%re), dd2q(z%im), ep)
    end function

    ! ── Level 1 — real ───────────────────────────────────────────────
    function target_ddot(n, x, incx, y, incy) result(r)
        integer,  intent(in) :: n, incx, incy
        real(ep), intent(in) :: x(*), y(*)
        real(ep) :: r
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx))
        yt = q2dd(y(1:ny))
        r = dd2q(dddot(n, xt, incx, yt, incy))
    end function

    function target_dasum(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        type(float64x2), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        r = dd2q(ddasum(n, xt, incx))
    end function

    function target_dnrm2(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        type(float64x2), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        r = dd2q(ddnrm2(n, xt, incx))
    end function

    function target_idamax(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        integer :: r
        integer :: nx
        type(float64x2), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        r = iddamax(n, xt, incx)
    end function

    subroutine target_daxpy(n, alpha, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        type(float64x2) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        alpha_t = q2dd(alpha)
        call ddaxpy(n, alpha_t, xt, incx, yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dcopy(n, x, incx, y, incy)
        integer,  intent(in)  :: n, incx, incy
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx))
        call ddcopy(n, xt, incx, yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dscal(n, alpha, x, incx)
        integer,  intent(in)    :: n, incx
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        integer :: nx
        type(float64x2), allocatable :: xt(:)
        type(float64x2) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        alpha_t = q2dd(alpha)
        call ddscal(n, alpha_t, xt, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    subroutine target_dswap(n, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddswap(n, xt, incx, yt, incy)
        x(1:nx) = dd2q(xt)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_drot(n, x, incx, y, incy, c, s)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: c, s
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        type(float64x2) :: ct, st
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        ct = q2dd(c); st = q2dd(s)
        call ddrot(n, xt, incx, yt, incy, ct, st)
        x(1:nx) = dd2q(xt)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_drotg(a, b, c, s)
        real(ep), intent(inout) :: a, b
        real(ep), intent(out)   :: c, s
        type(float64x2) :: at, bt, ct, st
        at = q2dd(a); bt = q2dd(b)
        call ddrotg(at, bt, ct, st)
        a = dd2q(at); b = dd2q(bt)
        c = dd2q(ct); s = dd2q(st)
    end subroutine

    subroutine target_drotm(n, x, incx, y, incy, param)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: param(5)
        integer :: nx, ny
        type(float64x2), allocatable :: xt(:), yt(:)
        type(float64x2) :: pt(5)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        pt = q2dd(param)
        call ddrotm(n, xt, incx, yt, incy, pt)
        x(1:nx) = dd2q(xt)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_drotmg(d1, d2, x1, y1, param)
        real(ep), intent(inout) :: d1, d2, x1
        real(ep), intent(in)    :: y1
        real(ep), intent(out)   :: param(5)
        type(float64x2) :: d1t, d2t, x1t, y1t, pt(5)
        d1t = q2dd(d1); d2t = q2dd(d2); x1t = q2dd(x1); y1t = q2dd(y1)
        call ddrotmg(d1t, d2t, x1t, y1t, pt)
        d1 = dd2q(d1t); d2 = dd2q(d2t); x1 = dd2q(x1t)
        param = dd2q(pt)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    function target_zdotc(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        integer :: nx, ny
        type(complex64x2), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2zz(x(1:nx))
        yt = q2zz(y(1:ny))
        r = zz2q(zzdotc(n, xt, incx, yt, incy))
    end function

    function target_zdotu(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        integer :: nx, ny
        type(complex64x2), allocatable :: xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2zz(x(1:nx))
        yt = q2zz(y(1:ny))
        r = zz2q(zzdotu(n, xt, incx, yt, incy))
    end function

    function target_dzasum(n, x, incx) result(r)
        integer,     intent(in) :: n, incx
        complex(ep), intent(in) :: x(*)
        real(ep) :: r
        integer :: nx
        type(complex64x2), allocatable :: xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2zz(x(1:nx))
        r = dd2q(ddzzasum(n, xt, incx))
    end function

    subroutine target_zaxpy(n, alpha, x, incx, y, incy)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(complex64x2), allocatable :: xt(:), yt(:)
        type(complex64x2) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(xt(nx), yt(ny))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        alpha_t = q2zz(alpha)
        call zzaxpy(n, alpha_t, xt, incx, yt, incy)
        y(1:ny) = zz2q(yt)
    end subroutine

    subroutine target_zscal(n, alpha, x, incx)
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: x(*)
        integer :: nx
        type(complex64x2), allocatable :: xt(:)
        type(complex64x2) :: alpha_t
        nx = (n - 1) * abs(incx) + 1
        allocate(xt(nx))
        xt = q2zz(x(1:nx))
        alpha_t = q2zz(alpha)
        call zzscal(n, alpha_t, xt, incx)
        x(1:nx) = zz2q(xt)
    end subroutine

    ! ── Level 2 — real ───────────────────────────────────────────────
    subroutine target_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
        end if
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddgemv(trans, m, n, q2dd(alpha), At, lda, xt, incx, &
                    q2dd(beta), yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
        end if
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddgbmv(trans, m, n, kl, ku, q2dd(alpha), At, lda, &
                    xt, incx, q2dd(beta), yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dger(m, n, alpha, x, incx, y, incy, A, lda)
        integer,  intent(in)    :: m, n, incx, incy, lda
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(lda,*)
        integer :: nx, ny
        type(float64x2), allocatable :: At(:,:), xt(:), yt(:)
        nx = (m - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddger(m, n, q2dd(alpha), xt, incx, yt, incy, At, lda)
        A(1:lda, 1:n) = dd2q(At)
    end subroutine

    subroutine target_dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddsymv(uplo, n, q2dd(alpha), At, lda, xt, incx, &
                    q2dd(beta), yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, incx, incy
        real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: aps, nx, ny
        type(float64x2), allocatable :: apt(:), xt(:), yt(:)
        aps = n * (n + 1) / 2
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(apt(aps), xt(nx), yt(ny))
        apt = q2dd(ap(1:aps))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddspmv(uplo, n, q2dd(alpha), apt, xt, incx, &
                    q2dd(beta), yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, k, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: nx, ny
        type(float64x2), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ddsbmv(uplo, n, k, q2dd(alpha), At, lda, xt, incx, &
                    q2dd(beta), yt, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_dtrmv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        type(float64x2), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx))
        call ddtrmv(uplo, trans, diag, n, At, lda, xt, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    subroutine target_dtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, k, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        type(float64x2), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx))
        call ddtbmv(uplo, trans, diag, n, k, At, lda, xt, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    subroutine target_dtpmv(uplo, trans, diag, n, ap, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, incx
        real(ep),  intent(in)    :: ap(*)
        real(ep),  intent(inout) :: x(*)
        integer :: aps, nx
        type(float64x2), allocatable :: apt(:), xt(:)
        aps = n * (n + 1) / 2
        nx = (n - 1) * abs(incx) + 1
        allocate(apt(aps), xt(nx))
        apt = q2dd(ap(1:aps))
        xt = q2dd(x(1:nx))
        call ddtpmv(uplo, trans, diag, n, apt, xt, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    subroutine target_dtrsv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        integer :: nx
        type(float64x2), allocatable :: At(:,:), xt(:)
        nx = (n - 1) * abs(incx) + 1
        allocate(At(lda, n), xt(nx))
        At = q2dd(A(1:lda, 1:n))
        xt = q2dd(x(1:nx))
        call ddtrsv(uplo, trans, diag, n, At, lda, xt, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    ! ── Level 2 — complex ────────────────────────────────────────────
    subroutine target_zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(complex64x2), allocatable :: At(:,:), xt(:), yt(:)
        if (trans == 'N' .or. trans == 'n') then
            nx = (n - 1) * abs(incx) + 1
            ny = (m - 1) * abs(incy) + 1
        else
            nx = (m - 1) * abs(incx) + 1
            ny = (n - 1) * abs(incy) + 1
        end if
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2zz(A(1:lda, 1:n))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call zzgemv(trans, m, n, q2zz(alpha), At, lda, xt, incx, &
                    q2zz(beta), yt, incy)
        y(1:ny) = zz2q(yt)
    end subroutine

    subroutine target_zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(complex64x2), allocatable :: At(:,:), xt(:), yt(:)
        nx = (n - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2zz(A(1:lda, 1:n))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call zzhemv(uplo, n, q2zz(alpha), At, lda, xt, incx, &
                    q2zz(beta), yt, incy)
        y(1:ny) = zz2q(yt)
    end subroutine

    subroutine target_zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        integer,     intent(in)    :: m, n, incx, incy, lda
        complex(ep), intent(in)    :: alpha, x(*), y(*)
        complex(ep), intent(inout) :: A(lda,*)
        integer :: nx, ny
        type(complex64x2), allocatable :: At(:,:), xt(:), yt(:)
        nx = (m - 1) * abs(incx) + 1
        ny = (n - 1) * abs(incy) + 1
        allocate(At(lda, n), xt(nx), yt(ny))
        At = q2zz(A(1:lda, 1:n))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call zzgerc(m, n, q2zz(alpha), xt, incx, yt, incy, At, lda)
        A(1:lda, 1:n) = zz2q(At)
    end subroutine

    ! ── Level 3 — real ───────────────────────────────────────────────
    subroutine target_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an, bn
        type(float64x2), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (transa == 'N' .or. transa == 'n') then; an = k; else; an = m; end if
        if (transb == 'N' .or. transb == 'n') then; bn = n; else; bn = k; end if
        allocate(At(lda, an), Bt(ldb, bn), Ct(ldc, n))
        At = q2dd(A(1:lda, 1:an))
        Bt = q2dd(B(1:ldb, 1:bn))
        Ct = q2dd(C(1:ldc, 1:n))
        call ddgemm(transa, transb, m, n, k, q2dd(alpha), At, lda, &
                    Bt, ldb, q2dd(beta), Ct, ldc)
        C(1:ldc, 1:n) = dd2q(Ct)
    end subroutine

    subroutine target_dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        type(float64x2), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n), Ct(ldc, n))
        At = q2dd(A(1:lda, 1:an))
        Bt = q2dd(B(1:ldb, 1:n))
        Ct = q2dd(C(1:ldc, 1:n))
        call ddsymm(side, uplo, m, n, q2dd(alpha), At, lda, &
                    Bt, ldb, q2dd(beta), Ct, ldc)
        C(1:ldc, 1:n) = dd2q(Ct)
    end subroutine

    subroutine target_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        type(float64x2), allocatable :: At(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Ct(ldc, n))
        At = q2dd(A(1:lda, 1:an))
        Ct = q2dd(C(1:ldc, 1:n))
        call ddsyrk(uplo, trans, n, k, q2dd(alpha), At, lda, &
                    q2dd(beta), Ct, ldc)
        C(1:ldc, 1:n) = dd2q(Ct)
    end subroutine

    subroutine target_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        integer :: an
        type(float64x2), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, an), Ct(ldc, n))
        At = q2dd(A(1:lda, 1:an))
        Bt = q2dd(B(1:ldb, 1:an))
        Ct = q2dd(C(1:ldc, 1:n))
        call ddsyr2k(uplo, trans, n, k, q2dd(alpha), At, lda, &
                     Bt, ldb, q2dd(beta), Ct, ldc)
        C(1:ldc, 1:n) = dd2q(Ct)
    end subroutine

    subroutine target_dtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer :: an
        type(float64x2), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = q2dd(A(1:lda, 1:an))
        Bt = q2dd(B(1:ldb, 1:n))
        call ddtrmm(side, uplo, transa, diag, m, n, q2dd(alpha), &
                    At, lda, Bt, ldb)
        B(1:ldb, 1:n) = dd2q(Bt)
    end subroutine

    subroutine target_dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer :: an
        type(float64x2), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = q2dd(A(1:lda, 1:an))
        Bt = q2dd(B(1:ldb, 1:n))
        call ddtrsm(side, uplo, transa, diag, m, n, q2dd(alpha), &
                    At, lda, Bt, ldb)
        B(1:ldb, 1:n) = dd2q(Bt)
    end subroutine

    ! ── Level 3 — complex ────────────────────────────────────────────
    subroutine target_zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an, bn
        type(complex64x2), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (transa == 'N' .or. transa == 'n') then; an = k; else; an = m; end if
        if (transb == 'N' .or. transb == 'n') then; bn = n; else; bn = k; end if
        allocate(At(lda, an), Bt(ldb, bn), Ct(ldc, n))
        At = q2zz(A(1:lda, 1:an))
        Bt = q2zz(B(1:ldb, 1:bn))
        Ct = q2zz(C(1:ldc, 1:n))
        call zzgemm(transa, transb, m, n, k, q2zz(alpha), At, lda, &
                    Bt, ldb, q2zz(beta), Ct, ldc)
        C(1:ldc, 1:n) = zz2q(Ct)
    end subroutine

    subroutine target_zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: side, uplo
        integer,     intent(in)    :: m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an
        type(complex64x2), allocatable :: At(:,:), Bt(:,:), Ct(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n), Ct(ldc, n))
        At = q2zz(A(1:lda, 1:an))
        Bt = q2zz(B(1:ldb, 1:n))
        Ct = q2zz(C(1:ldc, 1:n))
        call zzhemm(side, uplo, m, n, q2zz(alpha), At, lda, &
                    Bt, ldb, q2zz(beta), Ct, ldc)
        C(1:ldc, 1:n) = zz2q(Ct)
    end subroutine

    subroutine target_zherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, lda, ldc
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: C(ldc,*)
        integer :: an
        type(complex64x2), allocatable :: At(:,:), Ct(:,:)
        if (trans == 'N' .or. trans == 'n') then; an = k; else; an = n; end if
        allocate(At(lda, an), Ct(ldc, n))
        At = q2zz(A(1:lda, 1:an))
        Ct = q2zz(C(1:ldc, 1:n))
        call zzherk(uplo, trans, n, k, q2dd(alpha), At, lda, &
                    q2dd(beta), Ct, ldc)
        C(1:ldc, 1:n) = zz2q(Ct)
    end subroutine

    subroutine target_ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character,   intent(in)    :: side, uplo, transa, diag
        integer,     intent(in)    :: m, n, lda, ldb
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer :: an
        type(complex64x2), allocatable :: At(:,:), Bt(:,:)
        if (side == 'L' .or. side == 'l') then; an = m; else; an = n; end if
        allocate(At(lda, an), Bt(ldb, n))
        At = q2zz(A(1:lda, 1:an))
        Bt = q2zz(B(1:ldb, 1:n))
        call zztrsm(side, uplo, transa, diag, m, n, q2zz(alpha), &
                    At, lda, Bt, ldb)
        B(1:ldb, 1:n) = zz2q(Bt)
    end subroutine

end module target_blas
