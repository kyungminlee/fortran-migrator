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

        ! ── Phase 18 — triangular refinement ─────────────────────────
        subroutine dtrrfs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, nrhs, lda, ldb, ldx
            real(ep),  intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrrfs

        subroutine ztrrfs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, nrhs, lda, ldb, ldx
            complex(ep), intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztrrfs

        subroutine dtprfs(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, nrhs, ldb, ldx
            real(ep),  intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtprfs

        subroutine ztprfs(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, nrhs, ldb, ldx
            complex(ep), intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztprfs

        subroutine dtbrfs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
            real(ep),  intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtbrfs

        subroutine ztbrfs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
            complex(ep), intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztbrfs

        ! ── Equilibration — real ─────────────────────────────────────
        subroutine dgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgeequ

        subroutine dgbequ(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, kl, ku, ldab
            real(ep), intent(in)  :: AB(ldab,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgbequ

        subroutine dpoequ(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,  intent(in)  :: n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: S(*), scond, amax
            integer,  intent(out) :: info
        end subroutine dpoequ

        subroutine dppequ(uplo, n, AP, S, scond, amax, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*)
            real(ep),  intent(out) :: S(*), scond, amax
            integer,   intent(out) :: info
        end subroutine dppequ

        subroutine dpbequ(uplo, n, kd, AB, ldab, S, scond, amax, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*)
            real(ep),  intent(out) :: S(*), scond, amax
            integer,   intent(out) :: info
        end subroutine dpbequ

        ! ── Equilibration — complex ──────────────────────────────────
        subroutine zgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgeequ

        subroutine zgbequ(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgbequ

        subroutine zpoequ(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpoequ

        subroutine zppequ(uplo, n, AP, S, scond, amax, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zppequ

        subroutine zpbequ(uplo, n, kd, AB, ldab, S, scond, amax, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpbequ

        ! ── Improved equilibration (equb) ────────────────────────────
        subroutine dgeequb(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgeequb

        subroutine dgbequb(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, kl, ku, ldab
            real(ep), intent(in)  :: AB(ldab,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgbequb

        subroutine dpoequb(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,  intent(in)  :: n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: S(*), scond, amax
            integer,  intent(out) :: info
        end subroutine dpoequb

        subroutine dsyequb(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: S(*), scond, amax, work(*)
            integer,   intent(out) :: info
        end subroutine dsyequb

        subroutine zgeequb(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgeequb

        subroutine zgbequb(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgbequb

        subroutine zpoequb(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpoequb

        subroutine zsyequb(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zsyequb

        subroutine zheequb(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zheequb

        ! ── Tridiagonal reduction ────────────────────────────────────
        subroutine dsytrd(uplo, n, A, lda, D, E, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: D(*), E(*), tau(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrd

        subroutine zhetrd(uplo, n, A, lda, D, E, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrd

        subroutine dsptrd(uplo, n, AP, D, E, tau, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: D(*), E(*), tau(*)
            integer,   intent(out)   :: info
        end subroutine dsptrd

        subroutine zhptrd(uplo, n, AP, D, E, tau, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tau(*)
            integer,     intent(out)   :: info
        end subroutine zhptrd

        subroutine dsbtrd(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
            import :: ep
            character, intent(in)    :: vect, uplo
            integer,   intent(in)    :: n, kd, ldab, ldq
            real(ep),  intent(inout) :: AB(ldab,*), Q(ldq,*)
            real(ep),  intent(out)   :: D(*), E(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbtrd

        subroutine zhbtrd(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
            import :: ep
            character,   intent(in)    :: vect, uplo
            integer,     intent(in)    :: n, kd, ldab, ldq
            complex(ep), intent(inout) :: AB(ldab,*), Q(ldq,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhbtrd

        ! ── Orthogonal/unitary Q generation/application for *trd ─────
        subroutine dorgtr(uplo, n, A, lda, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: tau(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorgtr

        subroutine zungtr(uplo, n, A, lda, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungtr

        subroutine dormtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, uplo, trans
            integer,   intent(in)    :: m, n, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormtr

        subroutine zunmtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, uplo, trans
            integer,     intent(in)    :: m, n, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmtr

        subroutine dopgtr(uplo, n, AP, tau, Q, ldq, work, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, ldq
            real(ep),  intent(in)  :: AP(*), tau(*)
            real(ep),  intent(out) :: Q(ldq,*), work(*)
            integer,   intent(out) :: info
        end subroutine dopgtr

        subroutine zupgtr(uplo, n, AP, tau, Q, ldq, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, ldq
            complex(ep), intent(in)  :: AP(*), tau(*)
            complex(ep), intent(out) :: Q(ldq,*), work(*)
            integer,     intent(out) :: info
        end subroutine zupgtr

        subroutine dopmtr(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, uplo, trans
            integer,   intent(in)    :: m, n, ldc
            real(ep),  intent(in)    :: AP(*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dopmtr

        subroutine zupmtr(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, uplo, trans
            integer,     intent(in)    :: m, n, ldc
            complex(ep), intent(in)    :: AP(*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zupmtr

        ! ── Bidiagonal reduction ─────────────────────────────────────
        subroutine dgebrd(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: D(*), E(*), tauq(*), taup(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgebrd

        subroutine zgebrd(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgebrd

        subroutine dorgbr(vect, m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: vect
            integer,   intent(in)    :: m, n, k, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: tau(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorgbr

        subroutine zungbr(vect, m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungbr

        subroutine dormbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: vect, side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormbr

        subroutine zunmbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect, side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmbr

        ! ── Generalized eigenvalue (packed / banded) ─────────────────
        subroutine dspgv(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, ldz
            real(ep),  intent(inout) :: AP(*), BP(*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dspgv

        subroutine zhpgv(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, ldz
            complex(ep), intent(inout) :: AP(*), BP(*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhpgv

        subroutine dsbgv(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                         work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz
            real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbgv

        subroutine zhbgv(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                         work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz
            complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhbgv

        ! ── Generalized eigenvalue D&C variants ──────────────────────
        subroutine dsygvd(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, lda, ldb, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: W(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsygvd

        subroutine zhegvd(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, lda, ldb, lwork, lrwork, liwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhegvd

        subroutine dspgvd(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, ldz, lwork, liwork
            real(ep),  intent(inout) :: AP(*), BP(*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dspgvd

        subroutine zhpgvd(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AP(*), BP(*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhpgvd

        subroutine dsbgvd(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, liwork
            real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsbgvd

        subroutine zhbgvd(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, &
                                          lrwork, liwork
            complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhbgvd

        ! ── Tridiagonal eigenvalue solvers ───────────────────────────
        subroutine dsterf(n, D, E, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: D(*), E(*)
            integer,  intent(out)   :: info
        end subroutine dsterf

        subroutine dsteqr(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character, intent(in)    :: compz
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsteqr

        subroutine zsteqr(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character,   intent(in)    :: compz
            integer,     intent(in)    :: n, ldz
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: Z(ldz,*)
            real(ep),    intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsteqr

        subroutine dstedc(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: compz
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dstedc

        subroutine zstedc(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, &
                          iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: compz
            integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: Z(ldz,*)
            real(ep),    intent(out)   :: rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zstedc

        subroutine dstebz(range, order, n, vl, vu, il, iu, abstol, D, E, &
                          m, nsplit, W, iblock, isplit, work, iwork, info)
            import :: ep
            character, intent(in)  :: range, order
            integer,   intent(in)  :: n, il, iu
            real(ep),  intent(in)  :: vl, vu, abstol, D(*), E(*)
            integer,   intent(out) :: m, nsplit
            real(ep),  intent(out) :: W(*), work(*)
            integer,   intent(out) :: iblock(*), isplit(*), iwork(*), info
        end subroutine dstebz

        ! ── Bidiagonal SVD ───────────────────────────────────────────
        subroutine dbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, &
                          C, ldc, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
            real(ep),  intent(inout) :: D(*), E(*), VT(ldvt,*), U(ldu,*), C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dbdsqr

        subroutine zbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, &
                          C, ldc, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: VT(ldvt,*), U(ldu,*), C(ldc,*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zbdsqr

        subroutine dbdsdc(uplo, compq, n, D, E, U, ldu, VT, ldvt, Q, IQ, &
                          work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo, compq
            integer,   intent(in)    :: n, ldu, ldvt
            real(ep),  intent(inout) :: D(*), E(*)
            real(ep),  intent(out)   :: U(ldu,*), VT(ldvt,*), Q(*), work(*)
            integer,   intent(out)   :: IQ(*), iwork(*), info
        end subroutine dbdsdc

        ! ── ddisna + balancing/back-transform ────────────────────────
        subroutine ddisna(job, m, n, D, sep, info)
            import :: ep
            character, intent(in)  :: job
            integer,   intent(in)  :: m, n
            real(ep),  intent(in)  :: D(*)
            real(ep),  intent(out) :: sep(*)
            integer,   intent(out) :: info
        end subroutine ddisna

        subroutine dgebal(job, n, A, lda, ilo, ihi, scale, info)
            import :: ep
            character, intent(in)    :: job
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ilo, ihi
            real(ep),  intent(out)   :: scale(*)
            integer,   intent(out)   :: info
        end subroutine dgebal

        subroutine zgebal(job, n, A, lda, ilo, ihi, scale, info)
            import :: ep
            character,   intent(in)    :: job
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ilo, ihi
            real(ep),    intent(out)   :: scale(*)
            integer,     intent(out)   :: info
        end subroutine zgebal

        subroutine dgebak(job, side, n, ilo, ihi, scale, m, V, ldv, info)
            import :: ep
            character, intent(in)    :: job, side
            integer,   intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),  intent(in)    :: scale(*)
            real(ep),  intent(inout) :: V(ldv,*)
            integer,   intent(out)   :: info
        end subroutine dgebak

        subroutine zgebak(job, side, n, ilo, ihi, scale, m, V, ldv, info)
            import :: ep
            character,   intent(in)    :: job, side
            integer,     intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),    intent(in)    :: scale(*)
            complex(ep), intent(inout) :: V(ldv,*)
            integer,     intent(out)   :: info
        end subroutine zgebak

        ! ── Generalized balance/back-transform/Hessenberg ────────────
        subroutine dggbal(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, &
                          work, info)
            import :: ep
            character, intent(in)    :: job
            integer,   intent(in)    :: n, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ilo, ihi
            real(ep),  intent(out)   :: lscale(*), rscale(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dggbal

        subroutine zggbal(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, &
                          work, info)
            import :: ep
            character,   intent(in)    :: job
            integer,     intent(in)    :: n, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ilo, ihi
            real(ep),    intent(out)   :: lscale(*), rscale(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggbal

        subroutine dggbak(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
            import :: ep
            character, intent(in)    :: job, side
            integer,   intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),  intent(in)    :: lscale(*), rscale(*)
            real(ep),  intent(inout) :: V(ldv,*)
            integer,   intent(out)   :: info
        end subroutine dggbak

        subroutine zggbak(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
            import :: ep
            character,   intent(in)    :: job, side
            integer,     intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),    intent(in)    :: lscale(*), rscale(*)
            complex(ep), intent(inout) :: V(ldv,*)
            integer,     intent(out)   :: info
        end subroutine zggbak

        subroutine dgghrd(compq, compz, n, ilo, ihi, A, lda, B, ldb, &
                          Q, ldq, Z, ldz, info)
            import :: ep
            character, intent(in)    :: compq, compz
            integer,   intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            integer,   intent(out)   :: info
        end subroutine dgghrd

        subroutine zgghrd(compq, compz, n, ilo, ihi, A, lda, B, ldb, &
                          Q, ldq, Z, ldz, info)
            import :: ep
            character,   intent(in)    :: compq, compz
            integer,     intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            integer,     intent(out)   :: info
        end subroutine zgghrd

        ! ── Generalized QR / RQ factorizations ───────────────────────
        subroutine dggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, m, p, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggqrf

        subroutine zggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, m, p, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggqrf

        subroutine dggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, n, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggrqf

        subroutine zggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggrqf

        ! ── Sylvester equation ───────────────────────────────────────
        subroutine dtrsyl(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, &
                          scale, info)
            import :: ep
            character, intent(in)    :: trana, tranb
            integer,   intent(in)    :: isgn, m, n, lda, ldb, ldc
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: scale
            integer,   intent(out)   :: info
        end subroutine dtrsyl

        subroutine ztrsyl(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, &
                          scale, info)
            import :: ep
            character,   intent(in)    :: trana, tranb
            integer,     intent(in)    :: isgn, m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
            real(ep),    intent(out)   :: scale
            integer,     intent(out)   :: info
        end subroutine ztrsyl

        ! ── Hessenberg Schur factorization ───────────────────────────
        subroutine dhseqr(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: job, compz
            integer,   intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
            real(ep),  intent(inout) :: H(ldh,*), Z(ldz,*)
            real(ep),  intent(out)   :: WR(*), WI(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dhseqr

        subroutine zhseqr(job, compz, n, ilo, ihi, H, ldh, W, Z, ldz, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: job, compz
            integer,     intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
            complex(ep), intent(inout) :: H(ldh,*), Z(ldz,*)
            complex(ep), intent(out)   :: W(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhseqr

        ! ── Eigenvectors of (quasi-)triangular Schur form ────────────
        subroutine dtrevc(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          mm, m, work, info)
            import :: ep
            character, intent(in)    :: side, howmny
            integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm
            logical,   intent(inout) :: sel(*)
            real(ep),  intent(in)    :: T(ldt,*)
            real(ep),  intent(inout) :: VL(ldvl,*), VR(ldvr,*)
            integer,   intent(out)   :: m
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtrevc

        subroutine ztrevc(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          mm, m, work, rwork, info)
            import :: ep
            character,   intent(in)    :: side, howmny
            integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm
            logical,     intent(in)    :: sel(*)
            complex(ep), intent(inout) :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
            integer,     intent(out)   :: m
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine ztrevc

        ! ── Generalized non-symmetric eigenvalue ─────────────────────
        subroutine dggev(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, &
                         VL, ldvl, VR, ldvr, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobvl, jobvr
            integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dggev

        subroutine zggev(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, &
                         VL, ldvl, VR, ldvr, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zggev

        ! ── Modern (LAPACK 3.7+) blocked QR/LQ ───────────────────────
        subroutine dgeqr(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, tsize, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqr

        subroutine zgeqr(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, tsize, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqr

        subroutine dgelq(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, tsize, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelq

        subroutine zgelq(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, tsize, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelq

        subroutine dgemqr(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), T(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemqr

        subroutine zgemqr(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), T(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemqr

        subroutine dgemlq(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), T(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemlq

        subroutine zgemlq(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), T(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemlq

        ! ── Recursive blocked QR/LQ with explicit T (LAPACK 3.6+) ────
        subroutine dgeqrt(m, n, nb, A, lda, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, nb, lda, ldt
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqrt

        subroutine zgeqrt(m, n, nb, A, lda, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, nb, lda, ldt
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrt

        subroutine dgelqt(m, n, mb, A, lda, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb, lda, ldt
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelqt

        subroutine zgelqt(m, n, mb, A, lda, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb, lda, ldt
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelqt

        subroutine dgemqrt(side, trans, m, n, k, nb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, nb, ldv, ldt, ldc
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemqrt

        subroutine zgemqrt(side, trans, m, n, k, nb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, nb, ldv, ldt, ldc
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemqrt

        subroutine dgemlqt(side, trans, m, n, k, mb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, mb, ldv, ldt, ldc
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemlqt

        subroutine zgemlqt(side, trans, m, n, k, mb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, mb, ldv, ldt, ldc
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemlqt

        ! ── Tall-skinny least squares ────────────────────────────────
        subroutine dgetsls(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgetsls

        subroutine zgetsls(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetsls
    end interface

end module ref_quad_lapack
