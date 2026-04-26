! Explicit interfaces to reflapack_quad — the vendored Netlib LAPACK
! compiled with gfortran's -freal-8-real-16 so REAL(KIND=8) and
! DOUBLE PRECISION entities are promoted to REAL(KIND=16) without
! routine renaming. The D/Z routines below operate at quad precision
! and serve as the universal reference for the differential precision
! tests.
!
! Calling these via implicit typing would be wrong: a caller without
! the explicit interface would default DGESV etc. to DOUBLE PRECISION
! (KIND=8), mismatching the promoted KIND=16 symbol. Always `use
! ref_quad_lapack`.

module ref_quad_lapack
    use prec_kinds, only: ep
    implicit none

    interface
        ! ── Linear solve / factorization — real ──────────────────────
        subroutine dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, lda, ldb
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgesv

        subroutine dgetrf(m, n, A, lda, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgetrf

        subroutine dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgetrs

        subroutine dpotrf(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotrf

        subroutine dpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpotrs

        ! ── Linear solve — complex ───────────────────────────────────
        subroutine zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgesv

        subroutine zgetrf(m, n, A, lda, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgetrf

        subroutine zgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgetrs

        subroutine zpotrf(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotrf

        subroutine zpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpotrs

        subroutine zungqr(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungqr

        subroutine zgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobu, jobvt
            integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: s(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgesvd

        function zlange(norm, m, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: m, n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlange

        function zlanhe(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlanhe

        subroutine zlacpy(uplo, m, n, A, lda, B, ldb)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: m, n, lda, ldb
            complex(ep), intent(in)  :: A(lda,*)
            complex(ep), intent(out) :: B(ldb,*)
        end subroutine zlacpy

        subroutine zlaset(uplo, m, n, alpha, beta, A, lda)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zlaset

        ! ── QR factorization — real ──────────────────────────────────
        subroutine dgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqrf

        subroutine dorgqr(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgqr

        ! ── QR factorization — complex ───────────────────────────────
        subroutine zgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrf

        ! ── Symmetric / Hermitian eigenvalue ─────────────────────────
        subroutine dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsyev

        subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zheev

        ! ── SVD ──────────────────────────────────────────────────────
        subroutine dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgesvd

        ! ── Auxiliary ────────────────────────────────────────────────
        function dlange(norm, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlange

        subroutine dlacpy(uplo, m, n, A, lda, B, ldb)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, lda, ldb
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: B(ldb,*)
        end subroutine dlacpy

        subroutine dlaset(uplo, m, n, alpha, beta, A, lda)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, lda
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine dlaset

        ! ── Phase 2 — real user-facing drivers ───────────────────────
        subroutine dposv(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dposv

        subroutine dsysv(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsysv

        subroutine dsytrf(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsytrf

        subroutine dsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsytrs

        subroutine dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgels

        subroutine dsyevd(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyevd

        subroutine dgeev(jobvl, jobvr, n, A, lda, wr, wi, &
                         vl, ldvl, vr, ldvr, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobvl, jobvr
            integer,   intent(in)    :: n, lda, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgeev

        subroutine dgesdd(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgesdd

        subroutine dpotri(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotri

        subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgetri

        ! ── Phase 3 — triangular + complex driver mirrors ────────────
        subroutine dtrtri(uplo, diag, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo, diag
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dtrtri

        subroutine dtrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtrtrs

        subroutine ztrtri(uplo, diag, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo, diag
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine ztrtri

        subroutine ztrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztrtrs

        subroutine zposv(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zposv

        subroutine zhesv(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhesv

        subroutine zheevd(jobz, uplo, n, A, lda, w, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zheevd

        subroutine zgeev(jobvl, jobvr, n, A, lda, w, vl, ldvl, vr, ldvr, &
                         work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgeev

        subroutine zgesdd(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, rwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobz
            integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: s(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zgesdd

        subroutine zgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgels
    end interface

end module ref_quad_lapack
