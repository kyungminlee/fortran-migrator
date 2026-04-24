! Per-target ScaLAPACK wrapper for the kind16 (qscalapack) build.
!
! kind16 uses REAL(KIND=16) directly, so the wrappers are passthroughs
! to the migrated q-prefix (real) and x-prefix (complex) ScaLAPACK
! routines. Follows the same contract as tests/pblas/target_kind16 —
! test code calls target_pdXXX in quad precision, the wrappers hand
! the arrays straight through to the migrated pq*/px* entry points.
!
!   pdgesv  → pqgesv    pzgesv  → pxgesv
!   pdgetrf → pqgetrf   pdpotrf → pqpotrf
!   pdgeqrf → pqgeqrf   pzgeqrf → pxgeqrf
!   pdsyev  → pqsyev    pzheev  → pxheev
!   pdgesvd → pqgesvd   pdlange → pqlange
!   pdlacpy → pqlacpy   pdlaset → pqlaset

module target_scalapack
    use prec_kinds, only: ep
    implicit none
    private

    public :: target_name, target_eps
    public :: target_pdgesv, target_pdgetrf, target_pdgetrs
    public :: target_pdpotrf, target_pdpotrs
    public :: target_pdgeqrf
    public :: target_pdsyev
    public :: target_pdgesvd
    public :: target_pdlange
    public :: target_pdlacpy, target_pdlaset
    public :: target_pzgesv, target_pzgeqrf, target_pzheev

    character(len=*), parameter :: target_name = 'kind16'
    real(ep),         parameter :: target_eps  = epsilon(1.0_ep)

    interface
        ! ── Linear solve / LU — real ─────────────────────────────────
        subroutine pqgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                          B, ib, jb, descb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,  intent(in)    :: desca(9), descb(9)
            integer,  intent(out)   :: ipiv(*), info
            real(ep), intent(inout) :: A(*), B(*)
        end subroutine

        subroutine pqgetrf(m, n, A, ia, ja, desca, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, ia, ja
            integer,  intent(in)    :: desca(9)
            integer,  intent(out)   :: ipiv(*), info
            real(ep), intent(inout) :: A(*)
        end subroutine

        subroutine pqgetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                           B, ib, jb, descb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9), ipiv(*)
            integer,   intent(out)   :: info
            real(ep),  intent(in)    :: A(*)
            real(ep),  intent(inout) :: B(*)
        end subroutine

        ! ── Cholesky — real ──────────────────────────────────────────
        subroutine pqpotrf(uplo, n, A, ia, ja, desca, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ia, ja
            integer,   intent(in)    :: desca(9)
            integer,   intent(out)   :: info
            real(ep),  intent(inout) :: A(*)
        end subroutine

        subroutine pqpotrs(uplo, n, nrhs, A, ia, ja, desca, &
                           B, ib, jb, descb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            integer,   intent(out)   :: info
            real(ep),  intent(in)    :: A(*)
            real(ep),  intent(inout) :: B(*)
        end subroutine

        ! ── QR — real ────────────────────────────────────────────────
        subroutine pqgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, ia, ja, lwork
            integer,  intent(in)    :: desca(9)
            integer,  intent(out)   :: info
            real(ep), intent(inout) :: A(*)
            real(ep), intent(out)   :: tau(*), work(*)
        end subroutine

        ! ── Symmetric eigenvalue — real ──────────────────────────────
        subroutine pqsyev(jobz, uplo, n, A, ia, ja, desca, w, &
                          Z, iz, jz, descz, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ia, ja, iz, jz, lwork
            integer,   intent(in)    :: desca(9), descz(9)
            integer,   intent(out)   :: info
            real(ep),  intent(inout) :: A(*)
            real(ep),  intent(out)   :: w(*), Z(*), work(*)
        end subroutine

        ! ── SVD — real ───────────────────────────────────────────────
        subroutine pqgesvd(jobu, jobvt, m, n, A, ia, ja, desca, s, &
                           U, iu, ju, descu, VT, ivt, jvt, descvt, &
                           work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, ia, ja, iu, ju, ivt, jvt, lwork
            integer,   intent(in)    :: desca(9), descu(9), descvt(9)
            integer,   intent(out)   :: info
            real(ep),  intent(inout) :: A(*)
            real(ep),  intent(out)   :: s(*), U(*), VT(*), work(*)
        end subroutine

        ! ── Auxiliary — real ─────────────────────────────────────────
        function pqlange(norm, m, n, A, ia, ja, desca, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, ia, ja
            integer,   intent(in) :: desca(9)
            real(ep),  intent(in) :: A(*)
            real(ep)              :: work(*)
            real(ep) :: r
        end function

        subroutine pqlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, ia, ja, ib, jb
            integer,   intent(in)  :: desca(9), descb(9)
            real(ep),  intent(in)  :: A(*)
            real(ep),  intent(out) :: B(*)
        end subroutine

        subroutine pqlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, ia, ja
            integer,   intent(in)    :: desca(9)
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(inout) :: A(*)
        end subroutine

        ! ── Linear solve / QR — complex ──────────────────────────────
        subroutine pxgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                          B, ib, jb, descb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,     intent(in)    :: desca(9), descb(9)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(inout) :: A(*), B(*)
        end subroutine

        subroutine pxgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, ia, ja, lwork
            integer,     intent(in)    :: desca(9)
            integer,     intent(out)   :: info
            complex(ep), intent(inout) :: A(*)
            complex(ep), intent(out)   :: tau(*), work(*)
        end subroutine

        ! ── Hermitian eigenvalue — complex ───────────────────────────
        subroutine pxheev(jobz, uplo, n, A, ia, ja, desca, w, &
                          Z, iz, jz, descz, work, lwork, rwork, lrwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ia, ja, iz, jz, lwork, lrwork
            integer,     intent(in)    :: desca(9), descz(9)
            integer,     intent(out)   :: info
            complex(ep), intent(inout) :: A(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: Z(*), work(*)
        end subroutine
    end interface

contains

    subroutine target_pdgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,  intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,  intent(in)    :: desca(9), descb(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*), B(*)
        call pqgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                    B, ib, jb, descb, info)
    end subroutine

    subroutine target_pdgetrf(m, n, A, ia, ja, desca, ipiv, info)
        integer,  intent(in)    :: m, n, ia, ja
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*)
        call pqgetrf(m, n, A, ia, ja, desca, ipiv, info)
    end subroutine

    subroutine target_pdgetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                              B, ib, jb, descb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9), ipiv(*)
        integer,   intent(out)   :: info
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: B(*)
        call pqgetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                     B, ib, jb, descb, info)
    end subroutine

    subroutine target_pdpotrf(uplo, n, A, ia, ja, desca, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja
        integer,   intent(in)    :: desca(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        call pqpotrf(uplo, n, A, ia, ja, desca, info)
    end subroutine

    subroutine target_pdpotrs(uplo, n, nrhs, A, ia, ja, desca, &
                              B, ib, jb, descb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        integer,   intent(out)   :: info
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: B(*)
        call pqpotrs(uplo, n, nrhs, A, ia, ja, desca, &
                     B, ib, jb, descb, info)
    end subroutine

    subroutine target_pdgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, ia, ja, lwork
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: info
        real(ep), intent(inout) :: A(*)
        real(ep), intent(out)   :: tau(*), work(*)
        call pqgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
    end subroutine

    subroutine target_pdsyev(jobz, uplo, n, A, ia, ja, desca, w, &
                             Z, iz, jz, descz, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ia, ja, iz, jz, lwork
        integer,   intent(in)    :: desca(9), descz(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        real(ep),  intent(out)   :: w(*), Z(*), work(*)
        call pqsyev(jobz, uplo, n, A, ia, ja, desca, w, &
                    Z, iz, jz, descz, work, lwork, info)
    end subroutine

    subroutine target_pdgesvd(jobu, jobvt, m, n, A, ia, ja, desca, s, &
                              U, iu, ju, descu, VT, ivt, jvt, descvt, &
                              work, lwork, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, ia, ja, iu, ju, ivt, jvt, lwork
        integer,   intent(in)    :: desca(9), descu(9), descvt(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        real(ep),  intent(out)   :: s(*), U(*), VT(*), work(*)
        call pqgesvd(jobu, jobvt, m, n, A, ia, ja, desca, s, &
                     U, iu, ju, descu, VT, ivt, jvt, descvt, &
                     work, lwork, info)
    end subroutine

    function target_pdlange(norm, m, n, A, ia, ja, desca, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, ia, ja
        integer,   intent(in) :: desca(9)
        real(ep),  intent(in) :: A(*)
        real(ep)              :: work(*)
        real(ep) :: r
        r = pqlange(norm, m, n, A, ia, ja, desca, work)
    end function

    subroutine target_pdlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: m, n, ia, ja, ib, jb
        integer,   intent(in)  :: desca(9), descb(9)
        real(ep),  intent(in)  :: A(*)
        real(ep),  intent(out) :: B(*)
        call pqlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
    end subroutine

    subroutine target_pdlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, ia, ja
        integer,   intent(in)    :: desca(9)
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(*)
        call pqlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
    end subroutine

    subroutine target_pzgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,     intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,     intent(in)    :: desca(9), descb(9)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(inout) :: A(*), B(*)
        call pxgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                    B, ib, jb, descb, info)
    end subroutine

    subroutine target_pzgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, ia, ja, lwork
        integer,     intent(in)    :: desca(9)
        integer,     intent(out)   :: info
        complex(ep), intent(inout) :: A(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        call pxgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
    end subroutine

    subroutine target_pzheev(jobz, uplo, n, A, ia, ja, desca, w, &
                             Z, iz, jz, descz, work, lwork, rwork, lrwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, ia, ja, iz, jz, lwork, lrwork
        integer,     intent(in)    :: desca(9), descz(9)
        integer,     intent(out)   :: info
        complex(ep), intent(inout) :: A(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: Z(*), work(*)
        call pxheev(jobz, uplo, n, A, ia, ja, desca, w, &
                    Z, iz, jz, descz, work, lwork, rwork, lrwork, info)
    end subroutine

end module target_scalapack
