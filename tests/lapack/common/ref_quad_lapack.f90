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

        ! ── Phase 4 — selected/generalized eig + QR/LQ family ────────
        subroutine dsyevr(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dsyevr

        subroutine zheevr(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork, lrwork, liwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
        end subroutine zheevr

        subroutine dsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                         work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsygv

        subroutine zhegv(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                         work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhegv

        subroutine dgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, ilo, ihi, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgehrd

        subroutine dorghr(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, ilo, ihi, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorghr

        subroutine dgelqf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelqf

        subroutine dorglq(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorglq

        subroutine zgelqf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelqf

        subroutine zunglq(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunglq

        ! ── Phase 5 — complex mirrors + QR helpers ───────────────────
        subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetri

        subroutine zpotri(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotri

        subroutine zhetrf(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhetrf

        subroutine zhetrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhetrs

        subroutine zgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, ilo, ihi, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgehrd

        subroutine zunghr(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, ilo, ihi, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunghr

        subroutine dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(inout) :: jpvt(*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqp3

        subroutine zgeqp3(m, n, A, lda, jpvt, tau, work, lwork, rwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(inout) :: jpvt(*)
            complex(ep), intent(out)   :: tau(*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgeqp3

        subroutine dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormqr

        subroutine zunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmqr

        ! ── Phase 6 — banded + tridiagonal ───────────────────────────
        subroutine dgbtrf(m, n, kl, ku, AB, ldab, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, kl, ku, ldab
            real(ep), intent(inout) :: AB(ldab,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgbtrf

        subroutine dgbtrs(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgbtrs

        subroutine dgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgbsv

        subroutine dpbtrf(uplo, n, kd, AB, ldab, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, ldab
            real(ep),  intent(inout) :: AB(ldab,*)
            integer,   intent(out)   :: info
        end subroutine dpbtrf

        subroutine dpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpbtrs

        subroutine dpbsv(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpbsv

        subroutine dgttrf(n, dl, d, du, du2, ipiv, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: dl(*), d(*), du(*)
            real(ep), intent(out)   :: du2(*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgttrf

        subroutine dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: dl(*), d(*), du(*), du2(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgttrs

        subroutine dpttrf(n, d, e, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: d(*), e(*)
            integer,  intent(out)   :: info
        end subroutine dpttrf

        subroutine dpttrs(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(in)    :: d(*), e(*)
            real(ep), intent(inout) :: B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dpttrs
        ! ── Phase 7 — symmetric / Hermitian eigenvalue family ────────
        subroutine dsyevx(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, lwork, &
                          iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dsyevx

        subroutine zheevx(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, lwork, &
                          rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
        end subroutine zheevx

        subroutine dstev(jobz, n, d, e, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: d(*), e(*)
            real(ep),  intent(out)   :: z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dstev

        subroutine dstevd(jobz, n, d, e, z, ldz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: d(*), e(*)
            real(ep),  intent(out)   :: z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dstevd

        subroutine dstevx(jobz, range, n, d, e, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: d(*), e(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dstevx

        subroutine dstevr(jobz, range, n, d, e, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: d(*), e(*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dstevr

        subroutine dsbev(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbev

        subroutine dsbevd(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz, lwork, liwork
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsbevd

        subroutine dsbevx(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                          vl, vu, il, iu, abstol, m, w, z, ldz, &
                          work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: Q(ldq,*), w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsbevx

        subroutine zhbev(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                         work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhbev

        subroutine zhbevd(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhbevd

        subroutine zhbevx(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                          vl, vu, il, iu, abstol, m, w, z, ldz, &
                          work, rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: Q(ldq,*), z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhbevx

        subroutine dspev(jobz, uplo, n, AP, w, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dspev

        subroutine dspevd(jobz, uplo, n, AP, w, z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dspevd

        subroutine dspevx(jobz, range, uplo, n, AP, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dspevx

        subroutine zhpev(jobz, uplo, n, AP, w, z, ldz, work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ldz
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhpev

        subroutine zhpevd(jobz, uplo, n, AP, w, z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhpevd

        subroutine zhpevx(jobz, range, uplo, n, AP, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, rwork, &
                          iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, il, iu, ldz
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhpevx

        ! ── Phase 8 — QL / RQ family ─────────────────────────────────
        subroutine dgeqlf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqlf

        subroutine dorgql(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgql

        subroutine dormql(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormql

        subroutine zgeqlf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqlf

        subroutine zungql(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungql

        subroutine zunmql(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmql

        subroutine dgerqf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgerqf

        subroutine dorgrq(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgrq

        subroutine dormrq(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormrq

        subroutine zgerqf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgerqf

        subroutine zungrq(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungrq

        subroutine zunmrq(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmrq

        ! ── Phase 9 — real packed sym factor / solve / inverse ───────
        subroutine dsptrf(uplo, n, AP, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsptrf

        subroutine dsptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsptrs

        subroutine dspsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(inout) :: AP(*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dspsv

        subroutine dsptri(uplo, n, AP, ipiv, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsptri

        subroutine dpptrf(uplo, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dpptrf

        subroutine dpptrs(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpptrs

        subroutine dppsv(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(inout) :: AP(*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dppsv

        subroutine dpptri(uplo, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dpptri

        ! ── Phase 10 — complex packed factor / solve / inverse ───────
        subroutine zhptrf(uplo, n, AP, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhptrf

        subroutine zhptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhptrs

        subroutine zhpsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhpsv

        subroutine zhptri(uplo, n, AP, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhptri

        subroutine zsptrf(uplo, n, AP, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsptrf

        subroutine zsptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsptrs

        subroutine zspsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zspsv

        subroutine zsptri(uplo, n, AP, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsptri

        subroutine zpptrf(uplo, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine zpptrf

        subroutine zpptrs(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpptrs

        subroutine zppsv(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zppsv

        subroutine zpptri(uplo, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine zpptri

        ! ── Phase 11 — complex banded / tridiagonal ──────────────────
        subroutine zgbtrf(m, n, kl, ku, AB, ldab, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, kl, ku, ldab
            complex(ep), intent(inout) :: AB(ldab,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgbtrf

        subroutine zgbtrs(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgbtrs

        subroutine zgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgbsv

        subroutine zpbtrf(uplo, n, kd, AB, ldab, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, ldab
            complex(ep), intent(inout) :: AB(ldab,*)
            integer,     intent(out)   :: info
        end subroutine zpbtrf

        subroutine zpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpbtrs

        subroutine zpbsv(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpbsv

        subroutine zgttrf(n, dl, d, du, du2, ipiv, info)
            import :: ep
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: dl(*), d(*), du(*)
            complex(ep), intent(out)   :: du2(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgttrf

        subroutine zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: dl(*), d(*), du(*), du2(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgttrs

        subroutine zgtsv(n, nrhs, dl, d, du, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgtsv

        subroutine zpttrf(n, d, e, info)
            import :: ep
            integer,     intent(in)    :: n
            real(ep),    intent(inout) :: d(*)
            complex(ep), intent(inout) :: e(*)
            integer,     intent(out)   :: info
        end subroutine zpttrf

        subroutine zpttrs(uplo, n, nrhs, d, e, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            real(ep),    intent(in)    :: d(*)
            complex(ep), intent(in)    :: e(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpttrs

        subroutine zptsv(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            real(ep),    intent(inout) :: d(*)
            complex(ep), intent(inout) :: e(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zptsv

        ! ── Phase 12 — triangular packed / banded ────────────────────
        subroutine dtbtrs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtbtrs

        subroutine ztbtrs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztbtrs

        subroutine dtptrs(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtptrs

        subroutine ztptrs(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztptrs

        subroutine dtptri(uplo, diag, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo, diag
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dtptri

        subroutine ztptri(uplo, diag, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo, diag
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine ztptri

        ! ── Phase 13 — condition number estimators ───────────────────
        subroutine dgecon(norm, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgecon

        subroutine zgecon(norm, n, A, lda, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgecon

        subroutine dpocon(uplo, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dpocon

        subroutine zpocon(uplo, n, A, lda, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zpocon

        subroutine dgbcon(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, &
                          work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n, kl, ku, ldab
            real(ep),  intent(in)  :: AB(ldab,*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgbcon

        subroutine zgbcon(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, &
                          work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgbcon

        subroutine dgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, &
                          work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: dl(*), d(*), du(*), du2(*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgtcon

        subroutine zgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: dl(*), d(*), du(*), du2(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgtcon

        subroutine dsycon(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dsycon

        subroutine zhecon(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zhecon

        ! ── Phase 14 — condition estimators part 2 ───────────────────
        subroutine dpbcon(uplo, n, kd, AB, ldab, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dpbcon

        subroutine zpbcon(uplo, n, kd, AB, ldab, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zpbcon

        subroutine dppcon(uplo, n, AP, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dppcon

        subroutine zppcon(uplo, n, AP, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zppcon

        subroutine dptcon(n, d, e, anorm, rcond, work, info)
            import :: ep
            integer,  intent(in)  :: n
            real(ep), intent(in)  :: d(*), e(*), anorm
            real(ep), intent(out) :: rcond, work(*)
            integer,  intent(out) :: info
        end subroutine dptcon

        subroutine zptcon(n, d, e, anorm, rcond, rwork, info)
            import :: ep
            integer,     intent(in)  :: n
            real(ep),    intent(in)  :: d(*), anorm
            complex(ep), intent(in)  :: e(*)
            real(ep),    intent(out) :: rcond, rwork(*)
            integer,     intent(out) :: info
        end subroutine zptcon

        subroutine dspcon(uplo, n, AP, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dspcon

        subroutine zhpcon(uplo, n, AP, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zhpcon

        subroutine zspcon(uplo, n, AP, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zspcon

        subroutine dtrcon(norm, uplo, diag, n, A, lda, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrcon

        subroutine ztrcon(norm, uplo, diag, n, A, lda, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztrcon

        subroutine dtpcon(norm, uplo, diag, n, AP, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtpcon

        subroutine ztpcon(norm, uplo, diag, n, AP, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztpcon

        subroutine dtbcon(norm, uplo, diag, n, kd, AB, ldab, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtbcon

        subroutine ztbcon(norm, uplo, diag, n, kd, AB, ldab, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztbcon

        ! ── Phase 15 — iterative refinement ──────────────────────────
        subroutine dgerfs(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgerfs

        subroutine zgerfs(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgerfs

        subroutine dporfs(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dporfs

        subroutine zporfs(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zporfs

        ! ── Phase 16 — banded / tridiag refinement ───────────────────
        subroutine dgbrfs(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgbrfs

        subroutine zgbrfs(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgbrfs

        subroutine dpbrfs(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpbrfs

        subroutine zpbrfs(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zpbrfs

        subroutine dgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(in)    :: B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgtrfs

        subroutine zgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgtrfs

        subroutine dptrfs(n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          ferr, berr, work, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb, ldx
            real(ep), intent(in)    :: d(*), e(*), df(*), ef(*), B(ldb,*)
            real(ep), intent(inout) :: X(ldx,*)
            real(ep), intent(out)   :: ferr(*), berr(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dptrfs

        subroutine zptrfs(uplo, n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            real(ep),    intent(in)    :: d(*), df(*)
            complex(ep), intent(in)    :: e(*), ef(*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zptrfs

        ! ── Phase 17 — sym/Herm + packed refinement ──────────────────
        subroutine dsyrfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyrfs

        subroutine zsyrfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsyrfs

        subroutine zherfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zherfs

        subroutine dsprfs(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsprfs

        subroutine zsprfs(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsprfs

        subroutine dpprfs(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpprfs

        subroutine zpprfs(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zpprfs
    end interface

end module ref_quad_lapack
