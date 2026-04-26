! Per-target wrapper for the multifloats (ddlapack) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers split each
! quad value into the high (double approx) + low (double remainder)
! representation that TYPE(real64x2) expects, call the migrated
! dd*/zz* LAPACK routine, then recombine the two limbs back into quad
! for the comparison. This preserves multifloats' full ~32-decimal-digit
! precision through the wrapper boundary.
!
! Migrated routine prefixes for multifloats:
!   D → DD   (tgesv, tgetrf, tgeqrf, …)
!   Z → ZZ   (vgesv, vheev, vgeqrf, …)

module target_lapack
    use prec_kinds, only: ep, dp
    use multifloats, only: real64x2, cmplx64x2
    implicit none
    private

    public :: target_name, target_eps
    public :: target_dgesv, target_dgetrf, target_dgetrs
    public :: target_dpotrf, target_dpotrs
    public :: target_zgesv
    public :: target_dgeqrf, target_dorgqr, target_zgeqrf
    public :: target_dsyev, target_zheev
    public :: target_dgesvd
    public :: target_dlange, target_dlacpy, target_dlaset

    character(len=*), parameter :: target_name = 'multifloats'
    real(ep),         parameter :: target_eps  = real(2.0_dp**(-104), ep)

    interface
        subroutine tgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: real64x2
            integer, intent(in) :: n, nrhs, lda, ldb
            type(real64x2), intent(inout) :: A(lda,*), B(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine tgetrf(m, n, A, lda, ipiv, info)
            import :: real64x2
            integer, intent(in) :: m, n, lda
            type(real64x2), intent(inout) :: A(lda,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine tgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: real64x2
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb
            type(real64x2), intent(in)    :: A(lda,*)
            integer,         intent(in)    :: ipiv(*)
            type(real64x2), intent(inout) :: B(ldb,*)
            integer, intent(out) :: info
        end subroutine
        subroutine tpotrf(uplo, n, A, lda, info)
            import :: real64x2
            character, intent(in) :: uplo
            integer,   intent(in) :: n, lda
            type(real64x2), intent(inout) :: A(lda,*)
            integer, intent(out) :: info
        end subroutine
        subroutine tpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: real64x2
            character, intent(in) :: uplo
            integer,   intent(in) :: n, nrhs, lda, ldb
            type(real64x2), intent(in)    :: A(lda,*)
            type(real64x2), intent(inout) :: B(ldb,*)
            integer, intent(out) :: info
        end subroutine
        subroutine vgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: cmplx64x2
            integer, intent(in) :: n, nrhs, lda, ldb
            type(cmplx64x2), intent(inout) :: A(lda,*), B(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine tgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: real64x2
            integer, intent(in) :: m, n, lda, lwork
            type(real64x2), intent(inout) :: A(lda,*)
            type(real64x2), intent(out)   :: tau(*), work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine torgqr(m, n, k, A, lda, tau, work, lwork, info)
            import :: real64x2
            integer, intent(in) :: m, n, k, lda, lwork
            type(real64x2), intent(inout) :: A(lda,*)
            type(real64x2), intent(in)    :: tau(*)
            type(real64x2), intent(out)   :: work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine vgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: cmplx64x2
            integer, intent(in) :: m, n, lda, lwork
            type(cmplx64x2), intent(inout) :: A(lda,*)
            type(cmplx64x2), intent(out)   :: tau(*), work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine tsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: real64x2
            character, intent(in) :: jobz, uplo
            integer,   intent(in) :: n, lda, lwork
            type(real64x2), intent(inout) :: A(lda,*)
            type(real64x2), intent(out)   :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine vheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: real64x2, cmplx64x2
            character, intent(in) :: jobz, uplo
            integer,   intent(in) :: n, lda, lwork
            type(cmplx64x2), intent(inout) :: A(lda,*)
            type(real64x2),   intent(out)   :: w(*), rwork(*)
            type(cmplx64x2), intent(out)   :: work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine tgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                           work, lwork, info)
            import :: real64x2
            character, intent(in) :: jobu, jobvt
            integer,   intent(in) :: m, n, lda, ldu, ldvt, lwork
            type(real64x2), intent(inout) :: A(lda,*)
            type(real64x2), intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer, intent(out) :: info
        end subroutine
        function tlange(norm, m, n, A, lda, work) result(r)
            import :: real64x2
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            type(real64x2), intent(in) :: A(lda,*)
            type(real64x2) :: work(*)
            type(real64x2) :: r
        end function
        subroutine tlacpy(uplo, m, n, A, lda, B, ldb)
            import :: real64x2
            character, intent(in) :: uplo
            integer,   intent(in) :: m, n, lda, ldb
            type(real64x2), intent(in)  :: A(lda,*)
            type(real64x2), intent(out) :: B(ldb,*)
        end subroutine
        subroutine tlaset(uplo, m, n, alpha, beta, A, lda)
            import :: real64x2
            character, intent(in) :: uplo
            integer,   intent(in) :: m, n, lda
            type(real64x2), intent(in)    :: alpha, beta
            type(real64x2), intent(inout) :: A(lda,*)
        end subroutine
    end interface

contains

    ! ── Conversion helpers ──────────────────────────────────────────
    elemental function q2dd(x) result(r)
        real(ep), intent(in) :: x
        type(real64x2) :: r
        real(dp) :: hi
        hi = real(x, dp)
        r%limbs(1) = hi
        r%limbs(2) = real(x - real(hi, ep), dp)
    end function

    elemental function dd2q(x) result(r)
        type(real64x2), intent(in) :: x
        real(ep) :: r
        r = real(x%limbs(1), ep) + real(x%limbs(2), ep)
    end function

    elemental function q2zz(z) result(r)
        complex(ep), intent(in) :: z
        type(cmplx64x2) :: r
        r%re = q2dd(real(z, ep))
        r%im = q2dd(aimag(z))
    end function

    elemental function zz2q(z) result(r)
        type(cmplx64x2), intent(in) :: z
        complex(ep) :: r
        r = cmplx(dd2q(z%re), dd2q(z%im), ep)
    end function

    ! ── Linear solve / factorization — real ──────────────────────────
    subroutine target_dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, lda, ldb
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,  intent(out)   :: ipiv(*), info
        type(real64x2), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = q2dd(A(1:lda, 1:n))
        Bt = q2dd(B(1:ldb, 1:nrhs))
        call tgesv(n, nrhs, At, lda, ipiv, Bt, ldb, info)
        A(1:lda, 1:n)    = dd2q(At)
        B(1:ldb, 1:nrhs) = dd2q(Bt)
    end subroutine

    subroutine target_dgetrf(m, n, A, lda, ipiv, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(out)   :: ipiv(*), info
        type(real64x2), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = q2dd(A(1:lda, 1:n))
        call tgetrf(m, n, At, lda, ipiv, info)
        A(1:lda, 1:n) = dd2q(At)
    end subroutine

    subroutine target_dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = q2dd(A(1:lda, 1:n))
        Bt = q2dd(B(1:ldb, 1:nrhs))
        call tgetrs(trans, n, nrhs, At, lda, ipiv, Bt, ldb, info)
        B(1:ldb, 1:nrhs) = dd2q(Bt)
    end subroutine

    subroutine target_dpotrf(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        type(real64x2), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = q2dd(A(1:lda, 1:n))
        call tpotrf(uplo, n, At, lda, info)
        A(1:lda, 1:n) = dd2q(At)
    end subroutine

    subroutine target_dpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = q2dd(A(1:lda, 1:n))
        Bt = q2dd(B(1:ldb, 1:nrhs))
        call tpotrs(uplo, n, nrhs, At, lda, Bt, ldb, info)
        B(1:ldb, 1:nrhs) = dd2q(Bt)
    end subroutine

    subroutine target_zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        type(cmplx64x2), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = q2zz(A(1:lda, 1:n))
        Bt = q2zz(B(1:ldb, 1:nrhs))
        call vgesv(n, nrhs, At, lda, ipiv, Bt, ldb, info)
        A(1:lda, 1:n)    = zz2q(At)
        B(1:ldb, 1:nrhs) = zz2q(Bt)
    end subroutine

    ! ── QR factorization ─────────────────────────────────────────────
    subroutine target_dgeqrf(m, n, A, lda, tau, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*)
        integer,  intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), taut(:), work(:)
        type(real64x2) :: wopt(1)
        integer :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), taut(kmn))
        At = q2dd(A(1:lda, 1:n))
        call tgeqrf(m, n, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(dd2q(wopt(1))))
        allocate(work(lwork))
        call tgeqrf(m, n, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = dd2q(At)
        tau(1:kmn)    = dd2q(taut)
    end subroutine

    subroutine target_dorgqr(m, n, k, A, lda, tau, info)
        integer,  intent(in)    :: m, n, k, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        integer,  intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), taut(:), work(:)
        type(real64x2) :: wopt(1)
        integer :: lwork
        allocate(At(lda, n), taut(k))
        At = q2dd(A(1:lda, 1:n))
        taut = q2dd(tau(1:k))
        call torgqr(m, n, k, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(dd2q(wopt(1))))
        allocate(work(lwork))
        call torgqr(m, n, k, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = dd2q(At)
    end subroutine

    subroutine target_zgeqrf(m, n, A, lda, tau, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*)
        integer,     intent(out)   :: info
        type(cmplx64x2), allocatable :: At(:,:), taut(:), work(:)
        type(cmplx64x2) :: wopt(1)
        integer :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), taut(kmn))
        At = q2zz(A(1:lda, 1:n))
        call vgeqrf(m, n, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(real(zz2q(wopt(1)), ep)))
        allocate(work(lwork))
        call vgeqrf(m, n, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = zz2q(At)
        tau(1:kmn)    = zz2q(taut)
    end subroutine

    ! ── Eigenvalue ───────────────────────────────────────────────────
    subroutine target_dsyev(jobz, uplo, n, A, lda, w, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*)
        integer,   intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), wt(:), work(:)
        type(real64x2) :: wopt(1)
        integer :: lwork
        allocate(At(lda, n), wt(n))
        At = q2dd(A(1:lda, 1:n))
        call tsyev(jobz, uplo, n, At, lda, wt, wopt, -1, info)
        lwork = max(1, int(dd2q(wopt(1))))
        allocate(work(lwork))
        call tsyev(jobz, uplo, n, At, lda, wt, work, lwork, info)
        A(1:lda, 1:n) = dd2q(At)
        w(1:n)        = dd2q(wt)
    end subroutine

    subroutine target_zheev(jobz, uplo, n, A, lda, w, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*)
        integer,     intent(out)   :: info
        type(cmplx64x2), allocatable :: At(:,:), work(:)
        type(real64x2),   allocatable :: wt(:), rwork(:)
        type(cmplx64x2) :: wopt(1)
        integer :: lwork, lrwork
        lrwork = max(1, 3*n - 2)
        allocate(At(lda, n), wt(n), rwork(lrwork))
        At = q2zz(A(1:lda, 1:n))
        call vheev(jobz, uplo, n, At, lda, wt, wopt, -1, rwork, info)
        lwork = max(1, int(real(zz2q(wopt(1)), ep)))
        allocate(work(lwork))
        call vheev(jobz, uplo, n, At, lda, wt, work, lwork, rwork, info)
        A(1:lda, 1:n) = zz2q(At)
        w(1:n)        = dd2q(wt)
    end subroutine

    ! ── SVD ──────────────────────────────────────────────────────────
    subroutine target_dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, lda, ldu, ldvt
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*)
        integer,   intent(out)   :: info
        type(real64x2), allocatable :: At(:,:), Ut(:,:), VTt(:,:), st(:), work(:)
        type(real64x2) :: wopt(1)
        integer :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), Ut(ldu, m), VTt(ldvt, n), st(kmn))
        At = q2dd(A(1:lda, 1:n))
        call tgesvd(jobu, jobvt, m, n, At, lda, st, Ut, ldu, VTt, ldvt, &
                     wopt, -1, info)
        lwork = max(1, int(dd2q(wopt(1))))
        allocate(work(lwork))
        call tgesvd(jobu, jobvt, m, n, At, lda, st, Ut, ldu, VTt, ldvt, &
                     work, lwork, info)
        A(1:lda, 1:n)    = dd2q(At)
        s(1:kmn)         = dd2q(st)
        U(1:ldu, 1:m)    = dd2q(Ut)
        VT(1:ldvt, 1:n)  = dd2q(VTt)
    end subroutine

    ! ── Auxiliary ────────────────────────────────────────────────────
    function target_dlange(norm, m, n, A, lda) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: r
        type(real64x2), allocatable :: At(:,:), work(:)
        allocate(At(lda, n), work(max(1, m)))
        At = q2dd(A(1:lda, 1:n))
        r = dd2q(tlange(norm, m, n, At, lda, work))
    end function

    subroutine target_dlacpy(uplo, m, n, A, lda, B, ldb)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        type(real64x2), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, n))
        At = q2dd(A(1:lda, 1:n))
        Bt = q2dd(B(1:ldb, 1:n))
        call tlacpy(uplo, m, n, At, lda, Bt, ldb)
        B(1:ldb, 1:n) = dd2q(Bt)
    end subroutine

    subroutine target_dlaset(uplo, m, n, alpha, beta, A, lda)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(lda,*)
        type(real64x2), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = q2dd(A(1:lda, 1:n))
        call tlaset(uplo, m, n, q2dd(alpha), q2dd(beta), At, lda)
        A(1:lda, 1:n) = dd2q(At)
    end subroutine

end module target_lapack
