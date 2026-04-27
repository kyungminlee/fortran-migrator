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

        subroutine dgeqr2(m, n, A, lda, tau, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqr2

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

        subroutine dormlq(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormlq

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

        subroutine dgels(trans, m, n, nrhs, A, lda, B, ldb, &
                         work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgels

        subroutine dposv(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dposv

        subroutine dpotri(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotri

        subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, lda, lwork
            integer,  intent(in)    :: ipiv(*)
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dgetri

        subroutine dgebrd(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: d(*), e(*), tauq(*), taup(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgebrd

        ! ── QR factorization — complex ───────────────────────────────
        subroutine zgeqrf(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrf

        subroutine zgeqr2(m, n, A, lda, tau, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqr2

        subroutine zungqr(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungqr

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

        subroutine zunmlq(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmlq

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

        subroutine zgels(trans, m, n, nrhs, A, lda, B, ldb, &
                         work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgels

        subroutine zposv(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zposv

        subroutine zpotrf(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotrf

        subroutine zgetrf(m, n, A, lda, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgetrf

        subroutine zpotri(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotri

        subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, lda, lwork
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetri

        subroutine zgebrd(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: d(*), e(*)
            complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgebrd

        ! ── Symmetric / Hermitian eigenvalue ─────────────────────────
        subroutine dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsyev

        subroutine dsyevd(jobz, uplo, n, A, lda, w, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyevd

        subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zheev

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
    end interface

end module ref_quad_lapack
