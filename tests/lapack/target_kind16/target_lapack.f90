! Per-target wrapper for the kind16 (qlapack) build.
!
! kind16 uses REAL(KIND=16) directly, so wrappers are passthroughs to
! the migrated q-prefix (real) and x-prefix (complex) LAPACK routines.
! The wrapper interface lets test programs call target_<routine>
! without knowing the target's prefix or precision.
!
! Migrated routine prefixes for kind16:
!   D → Q   (qgesv, qgetrf, qgeqrf, …)
!   Z → X   (xgesv, xheev, xgeqrf, …)

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

    character(len=*), parameter :: target_name = 'kind16'
    real(ep),         parameter :: target_eps  = epsilon(1.0_ep)

    interface
        subroutine qgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, lda, ldb
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine
        subroutine qgetrf(m, n, A, lda, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine
        subroutine qgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine qpotrf(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine qpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine xgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine
        subroutine qgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine
        subroutine qorgqr(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine
        subroutine xgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine
        subroutine qsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine
        subroutine xheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine
        subroutine qgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: info
        end subroutine
        function qlange(norm, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function
        subroutine qlacpy(uplo, m, n, A, lda, B, ldb)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, lda, ldb
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: B(ldb,*)
        end subroutine
        subroutine qlaset(uplo, m, n, alpha, beta, A, lda)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, lda
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine
    end interface

contains

    ! ── Linear solve / factorization — real ──────────────────────────
    subroutine target_dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, lda, ldb
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,  intent(out)   :: ipiv(*), info
        call qgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine

    subroutine target_dgetrf(m, n, A, lda, ipiv, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(out)   :: ipiv(*), info
        call qgetrf(m, n, A, lda, ipiv, info)
    end subroutine

    subroutine target_dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call qgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine

    subroutine target_dpotrf(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call qpotrf(uplo, n, A, lda, info)
    end subroutine

    subroutine target_dpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call qpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine

    subroutine target_zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call xgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine

    ! ── QR factorization (workspace managed internally) ──────────────
    subroutine target_dgeqrf(m, n, A, lda, tau, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*)
        integer,  intent(out)   :: info
        real(ep), allocatable :: work(:)
        real(ep) :: wopt(1)
        integer  :: lwork
        call qgeqrf(m, n, A, lda, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call qgeqrf(m, n, A, lda, tau, work, lwork, info)
    end subroutine

    subroutine target_dorgqr(m, n, k, A, lda, tau, info)
        integer,  intent(in)    :: m, n, k, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        integer,  intent(out)   :: info
        real(ep), allocatable :: work(:)
        real(ep) :: wopt(1)
        integer  :: lwork
        call qorgqr(m, n, k, A, lda, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call qorgqr(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine

    subroutine target_zgeqrf(m, n, A, lda, tau, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*)
        integer,     intent(out)   :: info
        complex(ep), allocatable :: work(:)
        complex(ep) :: wopt(1)
        integer     :: lwork
        call xgeqrf(m, n, A, lda, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call xgeqrf(m, n, A, lda, tau, work, lwork, info)
    end subroutine

    ! ── Eigenvalue (workspace managed internally) ────────────────────
    subroutine target_dsyev(jobz, uplo, n, A, lda, w, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*)
        integer,   intent(out)   :: info
        real(ep), allocatable :: work(:)
        real(ep) :: wopt(1)
        integer  :: lwork
        call qsyev(jobz, uplo, n, A, lda, w, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call qsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
    end subroutine

    subroutine target_zheev(jobz, uplo, n, A, lda, w, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*)
        integer,     intent(out)   :: info
        complex(ep), allocatable :: work(:)
        real(ep),    allocatable :: rwork(:)
        complex(ep) :: wopt(1)
        integer     :: lwork, lrwork
        lrwork = max(1, 3*n - 2)
        allocate(rwork(lrwork))
        call xheev(jobz, uplo, n, A, lda, w, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call xheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
    end subroutine

    ! ── SVD (workspace managed internally) ───────────────────────────
    subroutine target_dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, lda, ldu, ldvt
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*)
        integer,   intent(out)   :: info
        real(ep), allocatable :: work(:)
        real(ep) :: wopt(1)
        integer  :: lwork
        call qgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                    wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call qgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                    work, lwork, info)
    end subroutine

    ! ── Auxiliary ────────────────────────────────────────────────────
    function target_dlange(norm, m, n, A, lda) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: r
        real(ep), allocatable :: work(:)
        allocate(work(max(1, m)))
        r = qlange(norm, m, n, A, lda, work)
    end function

    subroutine target_dlacpy(uplo, m, n, A, lda, B, ldb)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        call qlacpy(uplo, m, n, A, lda, B, ldb)
    end subroutine

    subroutine target_dlaset(uplo, m, n, alpha, beta, A, lda)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(lda,*)
        call qlaset(uplo, m, n, alpha, beta, A, lda)
    end subroutine

end module target_lapack
