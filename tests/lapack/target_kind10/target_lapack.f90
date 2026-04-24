! Per-target wrapper for the kind10 (elapack) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers cast quad
! inputs down to REAL(KIND=10) for the call into the migrated e/y
! prefix LAPACK routines, then cast results back up to quad.
!
! Migrated routine prefixes for kind10:
!   D → E   (egesv, egetrf, egeqrf, …)
!   Z → Y   (ygesv, yheev, ygeqrf, …)

module target_lapack
    use prec_kinds, only: ep
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

    integer,          parameter :: tk = 10
    character(len=*), parameter :: target_name = 'kind10'
    real(ep),         parameter :: target_eps  = real(epsilon(1.0_10), ep)

    interface
        subroutine egesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            integer,  intent(in)    :: n, nrhs, lda, ldb
            real(10), intent(inout) :: A(lda,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine
        subroutine egetrf(m, n, A, lda, ipiv, info)
            integer,  intent(in)    :: m, n, lda
            real(10), intent(inout) :: A(lda,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine
        subroutine egetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(10),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(10),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine epotrf(uplo, n, A, lda, info)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(10),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine epotrs(uplo, n, nrhs, A, lda, B, ldb, info)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(10),  intent(in)    :: A(lda,*)
            real(10),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine ygesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(10), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine
        subroutine egeqrf(m, n, A, lda, tau, work, lwork, info)
            integer,  intent(in)    :: m, n, lda, lwork
            real(10), intent(inout) :: A(lda,*)
            real(10), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine
        subroutine eorgqr(m, n, k, A, lda, tau, work, lwork, info)
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(10), intent(inout) :: A(lda,*)
            real(10), intent(in)    :: tau(*)
            real(10), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine
        subroutine ygeqrf(m, n, A, lda, tau, work, lwork, info)
            integer,     intent(in)    :: m, n, lda, lwork
            complex(10), intent(inout) :: A(lda,*)
            complex(10), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine
        subroutine esyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(10),  intent(inout) :: A(lda,*)
            real(10),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine yheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(10), intent(inout) :: A(lda,*)
            real(10),    intent(out)   :: w(*), rwork(*)
            complex(10), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine
        subroutine egesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, info)
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(10),  intent(inout) :: A(lda,*)
            real(10),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: info
        end subroutine
        function elange(norm, m, n, A, lda, work) result(r)
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            real(10),  intent(in) :: A(lda,*)
            real(10) :: work(*)
            real(10) :: r
        end function
        subroutine elacpy(uplo, m, n, A, lda, B, ldb)
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, lda, ldb
            real(10),  intent(in)  :: A(lda,*)
            real(10),  intent(out) :: B(ldb,*)
        end subroutine
        subroutine elaset(uplo, m, n, alpha, beta, A, lda)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, lda
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(inout) :: A(lda,*)
        end subroutine
    end interface

contains

    ! ── Linear solve / factorization — real ──────────────────────────
    subroutine target_dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, lda, ldb
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,  intent(out)   :: ipiv(*), info
        real(tk), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = real(A(1:lda, 1:n), tk)
        Bt = real(B(1:ldb, 1:nrhs), tk)
        call egesv(n, nrhs, At, lda, ipiv, Bt, ldb, info)
        A(1:lda, 1:n)    = real(At, ep)
        B(1:ldb, 1:nrhs) = real(Bt, ep)
    end subroutine

    subroutine target_dgetrf(m, n, A, lda, ipiv, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(out)   :: ipiv(*), info
        real(tk), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = real(A(1:lda, 1:n), tk)
        call egetrf(m, n, At, lda, ipiv, info)
        A(1:lda, 1:n) = real(At, ep)
    end subroutine

    subroutine target_dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        real(tk), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = real(A(1:lda, 1:n), tk)
        Bt = real(B(1:ldb, 1:nrhs), tk)
        call egetrs(trans, n, nrhs, At, lda, ipiv, Bt, ldb, info)
        B(1:ldb, 1:nrhs) = real(Bt, ep)
    end subroutine

    subroutine target_dpotrf(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        real(tk), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = real(A(1:lda, 1:n), tk)
        call epotrf(uplo, n, At, lda, info)
        A(1:lda, 1:n) = real(At, ep)
    end subroutine

    subroutine target_dpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        real(tk), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = real(A(1:lda, 1:n), tk)
        Bt = real(B(1:ldb, 1:nrhs), tk)
        call epotrs(uplo, n, nrhs, At, lda, Bt, ldb, info)
        B(1:ldb, 1:nrhs) = real(Bt, ep)
    end subroutine

    subroutine target_zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(tk), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, nrhs))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        Bt = cmplx(B(1:ldb, 1:nrhs), kind=tk)
        call ygesv(n, nrhs, At, lda, ipiv, Bt, ldb, info)
        A(1:lda, 1:n)    = cmplx(At, kind=ep)
        B(1:ldb, 1:nrhs) = cmplx(Bt, kind=ep)
    end subroutine

    ! ── QR factorization ─────────────────────────────────────────────
    subroutine target_dgeqrf(m, n, A, lda, tau, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*)
        integer,  intent(out)   :: info
        real(tk), allocatable :: At(:,:), taut(:), work(:)
        real(tk) :: wopt(1)
        integer  :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), taut(kmn))
        At = real(A(1:lda, 1:n), tk)
        call egeqrf(m, n, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call egeqrf(m, n, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = real(At, ep)
        tau(1:kmn)    = real(taut, ep)
    end subroutine

    subroutine target_dorgqr(m, n, k, A, lda, tau, info)
        integer,  intent(in)    :: m, n, k, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        integer,  intent(out)   :: info
        real(tk), allocatable :: At(:,:), taut(:), work(:)
        real(tk) :: wopt(1)
        integer  :: lwork
        allocate(At(lda, n), taut(k))
        At = real(A(1:lda, 1:n), tk)
        taut = real(tau(1:k), tk)
        call eorgqr(m, n, k, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call eorgqr(m, n, k, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = real(At, ep)
    end subroutine

    subroutine target_zgeqrf(m, n, A, lda, tau, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*)
        integer,     intent(out)   :: info
        complex(tk), allocatable :: At(:,:), taut(:), work(:)
        complex(tk) :: wopt(1)
        integer     :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), taut(kmn))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        call ygeqrf(m, n, At, lda, taut, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), tk)))
        allocate(work(lwork))
        call ygeqrf(m, n, At, lda, taut, work, lwork, info)
        A(1:lda, 1:n) = cmplx(At, kind=ep)
        tau(1:kmn)    = cmplx(taut, kind=ep)
    end subroutine

    ! ── Eigenvalue ───────────────────────────────────────────────────
    subroutine target_dsyev(jobz, uplo, n, A, lda, w, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*)
        integer,   intent(out)   :: info
        real(tk), allocatable :: At(:,:), wt(:), work(:)
        real(tk) :: wopt(1)
        integer  :: lwork
        allocate(At(lda, n), wt(n))
        At = real(A(1:lda, 1:n), tk)
        call esyev(jobz, uplo, n, At, lda, wt, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call esyev(jobz, uplo, n, At, lda, wt, work, lwork, info)
        A(1:lda, 1:n) = real(At, ep)
        w(1:n)        = real(wt, ep)
    end subroutine

    subroutine target_zheev(jobz, uplo, n, A, lda, w, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*)
        integer,     intent(out)   :: info
        complex(tk), allocatable :: At(:,:), work(:)
        real(tk),    allocatable :: wt(:), rwork(:)
        complex(tk) :: wopt(1)
        integer     :: lwork, lrwork
        lrwork = max(1, 3*n - 2)
        allocate(At(lda, n), wt(n), rwork(lrwork))
        At = cmplx(A(1:lda, 1:n), kind=tk)
        call yheev(jobz, uplo, n, At, lda, wt, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), tk)))
        allocate(work(lwork))
        call yheev(jobz, uplo, n, At, lda, wt, work, lwork, rwork, info)
        A(1:lda, 1:n) = cmplx(At, kind=ep)
        w(1:n)        = real(wt, ep)
    end subroutine

    ! ── SVD ──────────────────────────────────────────────────────────
    subroutine target_dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, lda, ldu, ldvt
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*)
        integer,   intent(out)   :: info
        real(tk), allocatable :: At(:,:), Ut(:,:), VTt(:,:), st(:), work(:)
        real(tk) :: wopt(1)
        integer  :: lwork, kmn
        kmn = min(m, n)
        allocate(At(lda, n), Ut(ldu, m), VTt(ldvt, n), st(kmn))
        At = real(A(1:lda, 1:n), tk)
        call egesvd(jobu, jobvt, m, n, At, lda, st, Ut, ldu, VTt, ldvt, &
                    wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call egesvd(jobu, jobvt, m, n, At, lda, st, Ut, ldu, VTt, ldvt, &
                    work, lwork, info)
        A(1:lda, 1:n)    = real(At, ep)
        s(1:kmn)         = real(st, ep)
        U(1:ldu, 1:m)    = real(Ut, ep)
        VT(1:ldvt, 1:n)  = real(VTt, ep)
    end subroutine

    ! ── Auxiliary ────────────────────────────────────────────────────
    function target_dlange(norm, m, n, A, lda) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: r
        real(tk), allocatable :: At(:,:), work(:)
        allocate(At(lda, n), work(max(1, m)))
        At = real(A(1:lda, 1:n), tk)
        r = real(elange(norm, m, n, At, lda, work), ep)
    end function

    subroutine target_dlacpy(uplo, m, n, A, lda, B, ldb)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        real(tk), allocatable :: At(:,:), Bt(:,:)
        allocate(At(lda, n), Bt(ldb, n))
        At = real(A(1:lda, 1:n), tk)
        Bt = real(B(1:ldb, 1:n), tk)
        call elacpy(uplo, m, n, At, lda, Bt, ldb)
        B(1:ldb, 1:n) = real(Bt, ep)
    end subroutine

    subroutine target_dlaset(uplo, m, n, alpha, beta, A, lda)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(lda,*)
        real(tk), allocatable :: At(:,:)
        allocate(At(lda, n))
        At = real(A(1:lda, 1:n), tk)
        call elaset(uplo, m, n, real(alpha, tk), real(beta, tk), At, lda)
        A(1:lda, 1:n) = real(At, ep)
    end subroutine

end module target_lapack
