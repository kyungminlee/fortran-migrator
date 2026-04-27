! Quad-precision reference implementations of the PTZBLAS auxiliary
! kernels that have no direct BLAS analogue.  Each routine matches the
! upstream PTZBLAS semantics (lifted from /external/ comment headers,
! never reading the migrated source) and produces a REAL(KIND=ep) /
! COMPLEX(KIND=ep) result for differential comparison against the
! migrated quad-precision implementation.
module ptzblas_ref_quad_aux
    use prec_kinds, only: ep
    implicit none
    private

    public :: ref_dset, ref_zset
    public :: ref_dascal
    public :: ref_zhescal
    public :: ref_dmmadd, ref_dmmtadd
    public :: ref_zmmadd, ref_zmmtadd
    public :: ref_drshft, ref_dcshft
    public :: ref_zrshft, ref_zcshft
    public :: ref_dasqrtb
    public :: ref_dtzpad, ref_ztzpad
    public :: ref_dtzscal, ref_ztzscal
    public :: ref_ztzcnjg
    public :: ref_dagemv, ref_xagemv
    public :: ref_dasymv, ref_zahemv, ref_zasymv
    public :: ref_datrmv, ref_zatrmv

contains

    ! ── x := alpha (scalar fill) ────────────────────────────────────
    subroutine ref_dset(n, alpha, x)
        integer,  intent(in)    :: n
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(:)
        x(1:n) = alpha
    end subroutine ref_dset

    subroutine ref_zset(n, alpha, x)
        integer,     intent(in)    :: n
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: x(:)
        x(1:n) = alpha
    end subroutine ref_zset

    ! ── x := |alpha| * |x| ──────────────────────────────────────────
    subroutine ref_dascal(n, alpha, x)
        integer,  intent(in)    :: n
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(:)
        x(1:n) = abs(alpha) * abs(x(1:n))
    end subroutine ref_dascal

    ! ── Hermitian scaling: A := alpha*A; zero imag of diagonal ──────
    ! Mirrors upstream PTZBLAS/zhescal.f.  ALPHA==1: only zero imag part
    ! of the diagonal (rest of A untouched).  ALPHA==0: delegate to
    ! ztzpad(UPLO,'N',...,ZERO,ZERO).  Otherwise: scale trapezoid by
    ! ALPHA (real) and force imag(diagonal)=0.
    subroutine ref_zhescal(uplo, m, n, ioffd, alpha, A)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: m, n, ioffd
        real(ep),    intent(in)    :: alpha
        complex(ep), intent(inout) :: A(:,:)
        integer     :: j, jtmp, mn
        complex(ep) :: zalpha, zzero
        if (m <= 0 .or. n <= 0) return
        zzero = cmplx(0.0_ep, 0.0_ep, ep)
        if (alpha == 1.0_ep) then
            if (uplo == 'L' .or. uplo == 'l' .or. &
                uplo == 'U' .or. uplo == 'u' .or. &
                uplo == 'D' .or. uplo == 'd') then
                do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                    jtmp = j + ioffd
                    A(jtmp, j) = cmplx(real(A(jtmp, j), ep), 0.0_ep, ep)
                end do
            end if
            return
        else if (alpha == 0.0_ep) then
            call ref_ztzpad(uplo, 'N', m, n, ioffd, zzero, zzero, A)
            return
        end if
        zalpha = cmplx(alpha, 0.0_ep, ep)
        ! General path: scale the slice, with imag(diagonal) forced to 0.
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            do j = 1, min(mn, n)
                A(1:m, j) = zalpha * A(1:m, j)
            end do
            do j = mn + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                A(jtmp, j) = cmplx(alpha * real(A(jtmp, j), ep), 0.0_ep, ep)
                if (m > jtmp) A(jtmp+1:m, j) = zalpha * A(jtmp+1:m, j)
            end do
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            do j = max(0, -ioffd) + 1, mn
                jtmp = j + ioffd
                if (jtmp > 1) A(1:jtmp-1, j) = zalpha * A(1:jtmp-1, j)
                A(jtmp, j) = cmplx(alpha * real(A(jtmp, j), ep), 0.0_ep, ep)
            end do
            do j = max(0, mn) + 1, n
                A(1:m, j) = zalpha * A(1:m, j)
            end do
        else if (uplo == 'D' .or. uplo == 'd') then
            do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                A(jtmp, j) = cmplx(alpha * real(A(jtmp, j), ep), 0.0_ep, ep)
            end do
        else
            A(1:m, 1:n) = zalpha * A(1:m, 1:n)
        end if
    end subroutine ref_zhescal

    ! ── B := alpha*A + beta*B ──────────────────────────────────────
    subroutine ref_dmmadd(m, n, alpha, A, beta, B)
        integer,  intent(in)    :: m, n
        real(ep), intent(in)    :: alpha, beta
        real(ep), intent(in)    :: A(:,:)
        real(ep), intent(inout) :: B(:,:)
        B(1:m, 1:n) = alpha * A(1:m, 1:n) + beta * B(1:m, 1:n)
    end subroutine ref_dmmadd

    subroutine ref_zmmadd(m, n, alpha, A, beta, B)
        integer,     intent(in)    :: m, n
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:)
        complex(ep), intent(inout) :: B(:,:)
        B(1:m, 1:n) = alpha * A(1:m, 1:n) + beta * B(1:m, 1:n)
    end subroutine ref_zmmadd

    ! ── B := alpha*A^T + beta*B (A is M×N; B is N×M) ───────────────
    ! Matches upstream qmmtadd / xmmtadd semantics: M = rows(A) = cols(B);
    ! N = cols(A) = rows(B).
    subroutine ref_dmmtadd(m, n, alpha, A, beta, B)
        integer,  intent(in)    :: m, n
        real(ep), intent(in)    :: alpha, beta
        real(ep), intent(in)    :: A(:,:)
        real(ep), intent(inout) :: B(:,:)
        B(1:n, 1:m) = alpha * transpose(A(1:m, 1:n)) + beta * B(1:n, 1:m)
    end subroutine ref_dmmtadd

    subroutine ref_zmmtadd(m, n, alpha, A, beta, B)
        integer,     intent(in)    :: m, n
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:)
        complex(ep), intent(inout) :: B(:,:)
        B(1:n, 1:m) = alpha * transpose(A(1:m, 1:n)) + beta * B(1:n, 1:m)
    end subroutine ref_zmmtadd

    ! ── Row / column shift (NOT circular — upstream q[r|c]shft simply
    !    copies the leading m×n block by `offset` rows / columns and
    !    leaves the source positions unchanged). For offset>0 with
    !    rshft, the rows 1+offset..m+offset get the original 1..m;
    !    rows 1..offset are untouched. The buffer must be at least
    !    (m+|offset|) × n / m × (n+|offset|) — the test allocates that
    !    extra slack. ─────────────────────────────────────────────
    subroutine ref_drshft(m, n, offset, A)
        integer,  intent(in)    :: m, n, offset
        real(ep), intent(inout) :: A(:,:)
        integer :: i, j
        if (offset == 0 .or. m <= 0 .or. n <= 0) return
        if (offset > 0) then
            do j = 1, n
                do i = m, 1, -1
                    A(i + offset, j) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i, j) = A(i - offset, j)
                end do
            end do
        end if
    end subroutine ref_drshft

    subroutine ref_dcshft(m, n, offset, A)
        integer,  intent(in)    :: m, n, offset
        real(ep), intent(inout) :: A(:,:)
        integer :: i, j
        if (offset == 0 .or. m <= 0 .or. n <= 0) return
        if (offset > 0) then
            do j = n, 1, -1
                do i = 1, m
                    A(i, j + offset) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i, j) = A(i, j - offset)
                end do
            end do
        end if
    end subroutine ref_dcshft

    subroutine ref_zrshft(m, n, offset, A)
        integer,     intent(in)    :: m, n, offset
        complex(ep), intent(inout) :: A(:,:)
        integer :: i, j
        if (offset == 0 .or. m <= 0 .or. n <= 0) return
        if (offset > 0) then
            do j = 1, n
                do i = m, 1, -1
                    A(i + offset, j) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i, j) = A(i - offset, j)
                end do
            end do
        end if
    end subroutine ref_zrshft

    subroutine ref_zcshft(m, n, offset, A)
        integer,     intent(in)    :: m, n, offset
        complex(ep), intent(inout) :: A(:,:)
        integer :: i, j
        if (offset == 0 .or. m <= 0 .or. n <= 0) return
        if (offset > 0) then
            do j = n, 1, -1
                do i = 1, m
                    A(i, j + offset) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i, j) = A(i, j - offset)
                end do
            end do
        end if
    end subroutine ref_zcshft

    ! ── c := a * sqrt(b) ───────────────────────────────────────────
    subroutine ref_dasqrtb(a, b, c)
        real(ep), intent(in)  :: a, b
        real(ep), intent(out) :: c
        c = a * sqrt(b)
    end subroutine ref_dasqrtb

    ! ── Trapezoidal pad: diag := beta, offdiag := alpha ────────────
    ! Mirrors upstream PTZBLAS/dtzpad.f exactly: a "full-fill" pass for
    ! columns where the trapezoid covers the entire column, then a
    ! "diagonal + strict-trapezoid" pass for the remaining columns.
    ! Note: upstream dtzpad has no UPLO='D' branch (real type, no imag
    ! to zero); leave that case as a no-op to match.
    subroutine ref_dtzpad(uplo, herm, m, n, ioffd, alpha, beta, A)
        character, intent(in)    :: uplo, herm
        integer,   intent(in)    :: m, n, ioffd
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(:,:)
        integer  :: i, j, jtmp, mn
        logical  :: zero_imag
        if (m <= 0 .or. n <= 0) return
        zero_imag = (herm == 'Z' .or. herm == 'z')
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            ! Full-fill leading columns (entire trapezoid covers column).
            do j = 1, min(mn, n)
                do i = 1, m
                    A(i, j) = alpha
                end do
            end do
            ! Diagonal + strict-lower for remaining columns.
            if (zero_imag) then
                do j = mn + 1, min(m - ioffd, n)
                    jtmp = j + ioffd
                    do i = jtmp + 1, m
                        A(i, j) = alpha
                    end do
                end do
            else
                do j = mn + 1, min(m - ioffd, n)
                    jtmp = j + ioffd
                    A(jtmp, j) = beta
                    do i = jtmp + 1, m
                        A(i, j) = alpha
                    end do
                end do
            end if
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            if (zero_imag) then
                do j = max(0, -ioffd) + 1, mn
                    jtmp = j + ioffd
                    do i = 1, jtmp - 1
                        A(i, j) = alpha
                    end do
                end do
            else
                do j = max(0, -ioffd) + 1, mn
                    jtmp = j + ioffd
                    do i = 1, jtmp - 1
                        A(i, j) = alpha
                    end do
                    A(jtmp, j) = beta
                end do
            end if
            ! Trailing columns: entire column is in upper trapezoid.
            do j = max(0, mn) + 1, n
                do i = 1, m
                    A(i, j) = alpha
                end do
            end do
        end if
        ! UPLO='D' for real: upstream dtzpad branch sets diagonal to beta
        ! when HERM is not 'Z' (and is unsupported for HERM='Z' since real
        ! has no imag).  Mirror that.
        if (uplo == 'D' .or. uplo == 'd') then
            if (.not. zero_imag) then
                if (ioffd < m .and. ioffd > -n) then
                    do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                        A(j + ioffd, j) = beta
                    end do
                end if
            end if
        end if
    end subroutine ref_dtzpad

    ! Mirrors upstream PTZBLAS/ztzpad.f exactly, including UPLO='D'.
    subroutine ref_ztzpad(uplo, herm, m, n, ioffd, alpha, beta, A)
        character,   intent(in)    :: uplo, herm
        integer,     intent(in)    :: m, n, ioffd
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(inout) :: A(:,:)
        integer :: i, j, jtmp, mn
        logical :: zero_imag
        if (m <= 0 .or. n <= 0) return
        zero_imag = (herm == 'Z' .or. herm == 'z')
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            do j = 1, min(mn, n)
                do i = 1, m
                    A(i, j) = alpha
                end do
            end do
            if (zero_imag) then
                do j = mn + 1, min(m - ioffd, n)
                    jtmp = j + ioffd
                    A(jtmp, j) = cmplx(real(A(jtmp, j), ep), 0.0_ep, ep)
                    do i = jtmp + 1, m
                        A(i, j) = alpha
                    end do
                end do
            else
                do j = mn + 1, min(m - ioffd, n)
                    jtmp = j + ioffd
                    A(jtmp, j) = beta
                    do i = jtmp + 1, m
                        A(i, j) = alpha
                    end do
                end do
            end if
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            if (zero_imag) then
                do j = max(0, -ioffd) + 1, mn
                    jtmp = j + ioffd
                    do i = 1, jtmp - 1
                        A(i, j) = alpha
                    end do
                    A(jtmp, j) = cmplx(real(A(jtmp, j), ep), 0.0_ep, ep)
                end do
            else
                do j = max(0, -ioffd) + 1, mn
                    jtmp = j + ioffd
                    do i = 1, jtmp - 1
                        A(i, j) = alpha
                    end do
                    A(jtmp, j) = beta
                end do
            end if
            do j = max(0, mn) + 1, n
                do i = 1, m
                    A(i, j) = alpha
                end do
            end do
        else if (uplo == 'D' .or. uplo == 'd') then
            ! Diagonal-only set (or zero-imag of diagonal).
            if (ioffd < m .and. ioffd > -n) then
                if (zero_imag) then
                    do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                        jtmp = j + ioffd
                        A(jtmp, j) = cmplx(real(A(jtmp, j), ep), 0.0_ep, ep)
                    end do
                else
                    do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                        A(j + ioffd, j) = beta
                    end do
                end if
            end if
        end if
    end subroutine ref_ztzpad

    ! ── Trapezoidal scaling A := alpha * A on UPLO/IOFFD slice ─────
    ! Mirrors upstream PTZBLAS/dtzscal.f (and ztzscal.f).
    subroutine ref_dtzscal(uplo, m, n, ioffd, alpha, A)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, ioffd
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(inout) :: A(:,:)
        integer :: i, j, jtmp, mn
        if (m <= 0 .or. n <= 0) return
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            do j = 1, min(mn, n)
                do i = 1, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
            do j = mn + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                do i = jtmp, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            do j = max(0, -ioffd) + 1, mn
                jtmp = j + ioffd
                do i = 1, jtmp
                    A(i, j) = alpha * A(i, j)
                end do
            end do
            do j = max(0, mn) + 1, n
                do i = 1, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
        else if (uplo == 'D' .or. uplo == 'd') then
            do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                A(jtmp, j) = alpha * A(jtmp, j)
            end do
        else
            A(1:m, 1:n) = alpha * A(1:m, 1:n)
        end if
    end subroutine ref_dtzscal

    subroutine ref_ztzscal(uplo, m, n, ioffd, alpha, A)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: m, n, ioffd
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: A(:,:)
        integer :: i, j, jtmp, mn
        if (m <= 0 .or. n <= 0) return
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            do j = 1, min(mn, n)
                do i = 1, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
            do j = mn + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                do i = jtmp, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            do j = max(0, -ioffd) + 1, mn
                jtmp = j + ioffd
                do i = 1, jtmp
                    A(i, j) = alpha * A(i, j)
                end do
            end do
            do j = max(0, mn) + 1, n
                do i = 1, m
                    A(i, j) = alpha * A(i, j)
                end do
            end do
        else if (uplo == 'D' .or. uplo == 'd') then
            do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                A(jtmp, j) = alpha * A(jtmp, j)
            end do
        else
            A(1:m, 1:n) = alpha * A(1:m, 1:n)
        end if
    end subroutine ref_ztzscal

    ! ── Trapezoidal conjugate-and-scale A := alpha*conj(A) ─────────
    ! Mirrors upstream PTZBLAS/ztzcnjg.f.
    subroutine ref_ztzcnjg(uplo, m, n, ioffd, alpha, A)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: m, n, ioffd
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: A(:,:)
        integer :: i, j, jtmp, mn
        if (m <= 0 .or. n <= 0) return
        if (uplo == 'L' .or. uplo == 'l') then
            mn = max(0, -ioffd)
            do j = 1, min(mn, n)
                do i = 1, m
                    A(i, j) = alpha * conjg(A(i, j))
                end do
            end do
            do j = mn + 1, min(m - ioffd, n)
                do i = j + ioffd, m
                    A(i, j) = alpha * conjg(A(i, j))
                end do
            end do
        else if (uplo == 'U' .or. uplo == 'u') then
            mn = min(m - ioffd, n)
            do j = max(0, -ioffd) + 1, mn
                do i = 1, j + ioffd
                    A(i, j) = alpha * conjg(A(i, j))
                end do
            end do
            do j = max(0, mn) + 1, n
                do i = 1, m
                    A(i, j) = alpha * conjg(A(i, j))
                end do
            end do
        else if (uplo == 'D' .or. uplo == 'd') then
            do j = max(0, -ioffd) + 1, min(m - ioffd, n)
                jtmp = j + ioffd
                A(jtmp, j) = alpha * conjg(A(jtmp, j))
            end do
        else
            A(1:m, 1:n) = alpha * conjg(A(1:m, 1:n))
        end if
    end subroutine ref_ztzcnjg

    ! ── y := |alpha| * |A|*|x| + |beta*y|  (TRANS='N')
    !        y := |alpha| * |A^T|*|x| + |beta*y|  (TRANS='T' or 'C') ─
    subroutine ref_dagemv(trans, m, n, alpha, A, x, beta, y)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(:,:), x(:)
        real(ep),  intent(inout) :: y(:)
        if (trans == 'N' .or. trans == 'n') then
            y(1:m) = abs(alpha) * matmul(abs(A(1:m,1:n)), abs(x(1:n))) &
                     + abs(beta * y(1:m))
        else
            y(1:n) = abs(alpha) * matmul(transpose(abs(A(1:m,1:n))), abs(x(1:m))) &
                     + abs(beta * y(1:n))
        end if
    end subroutine ref_dagemv

    ! Complex variant: A is complex, x is complex, y is real.
    ! Upstream xagemv uses CABS1(z) = |Re(z)| + |Im(z)| as the complex
    ! "absolute value", NOT the Euclidean magnitude — match that.
    subroutine ref_xagemv(trans, m, n, alpha, A, x, beta, y)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:), x(:)
        real(ep),    intent(inout) :: y(:)
        real(ep), allocatable :: aA(:,:), aX(:)
        integer :: i, j
        if (trans == 'N' .or. trans == 'n') then
            allocate(aA(m, n), aX(n))
            do j = 1, n; aX(j) = abs(real(x(j), ep)) + abs(aimag(x(j))); end do
            do j = 1, n; do i = 1, m
                aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
            end do; end do
            y(1:m) = abs(alpha) * matmul(aA, aX) + abs(beta * y(1:m))
            deallocate(aA, aX)
        else
            allocate(aA(m, n), aX(m))
            do i = 1, m; aX(i) = abs(real(x(i), ep)) + abs(aimag(x(i))); end do
            do j = 1, n; do i = 1, m
                aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
            end do; end do
            y(1:n) = abs(alpha) * matmul(transpose(aA), aX) + abs(beta * y(1:n))
            deallocate(aA, aX)
        end if
    end subroutine ref_xagemv

    ! Symmetric abs-matvec — UPLO selects which triangle is referenced;
    ! the missing triangle is implicitly mirrored.
    subroutine ref_dasymv(uplo, n, alpha, A, x, beta, y)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(:,:), x(:)
        real(ep),  intent(inout) :: y(:)
        real(ep), allocatable :: As(:,:)
        integer :: i, j
        allocate(As(n, n))
        As = 0.0_ep
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    As(i, j) = A(i, j)
                    As(j, i) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    As(i, j) = A(i, j)
                    As(j, i) = A(i, j)
                end do
            end do
        end if
        y(1:n) = abs(alpha) * matmul(abs(As), abs(x(1:n))) + abs(beta * y(1:n))
        deallocate(As)
    end subroutine ref_dasymv

    ! Hermitian variant (complex A). Upstream xahemv uses CABS1 for
    ! off-diagonal A and |Re(A_jj)| for the diagonal (Hermitian => imag
    ! diagonal is zero by convention).
    subroutine ref_zahemv(uplo, n, alpha, A, x, beta, y)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:), x(:)
        real(ep),    intent(inout) :: y(:)
        real(ep), allocatable :: aA(:,:), aX(:)
        integer :: i, j
        allocate(aA(n, n), aX(n))
        do j = 1, n; aX(j) = abs(real(x(j), ep)) + abs(aimag(x(j))); end do
        ! Build symmetric |A| using CABS1 off-diagonal, |Re| on diagonal.
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j - 1
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                    aA(j, i) = aA(i, j)
                end do
                aA(j, j) = abs(real(A(j, j), ep))
            end do
        else
            do j = 1, n
                aA(j, j) = abs(real(A(j, j), ep))
                do i = j + 1, n
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                    aA(j, i) = aA(i, j)
                end do
            end do
        end if
        y(1:n) = abs(alpha) * matmul(aA, aX) + abs(beta * y(1:n))
        deallocate(aA, aX)
    end subroutine ref_zahemv

    ! Symmetric (not Hermitian) abs-matvec for complex A. Upstream
    ! xasymv uses CABS1 throughout (no special diagonal handling).
    subroutine ref_zasymv(uplo, n, alpha, A, x, beta, y)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:), x(:)
        real(ep),    intent(inout) :: y(:)
        real(ep), allocatable :: aA(:,:), aX(:)
        integer :: i, j
        allocate(aA(n, n), aX(n))
        do j = 1, n; aX(j) = abs(real(x(j), ep)) + abs(aimag(x(j))); end do
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                    aA(j, i) = aA(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                    aA(j, i) = aA(i, j)
                end do
            end do
        end if
        y(1:n) = abs(alpha) * matmul(aA, aX) + abs(beta * y(1:n))
        deallocate(aA, aX)
    end subroutine ref_zasymv

    ! Triangular abs-matvec — the upstream qatrmv is in fact a *gemv*
    ! restricted to the trapezoidal slice selected by UPLO/DIAG; in
    ! abs-arithmetic that means we mask off the unused triangle, treat
    ! DIAG='U' as setting the diagonal of |A| to 1, and apply the same
    ! formula as agemv.
    subroutine ref_datrmv(uplo, trans, diag, n, alpha, A, x, beta, y)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(:,:), x(:)
        real(ep),  intent(inout) :: y(:)
        real(ep), allocatable :: At(:,:)
        integer :: i, j
        allocate(At(n, n))
        At = 0.0_ep
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    At(i, j) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    At(i, j) = A(i, j)
                end do
            end do
        end if
        if (diag == 'U' .or. diag == 'u') then
            do i = 1, n
                At(i, i) = 1.0_ep
            end do
        end if
        if (trans == 'N' .or. trans == 'n') then
            y(1:n) = abs(alpha) * matmul(abs(At), abs(x(1:n))) + abs(beta * y(1:n))
        else
            y(1:n) = abs(alpha) * matmul(transpose(abs(At)), abs(x(1:n))) &
                     + abs(beta * y(1:n))
        end if
        deallocate(At)
    end subroutine ref_datrmv

    subroutine ref_zatrmv(uplo, trans, diag, n, alpha, A, x, beta, y)
        character,   intent(in)    :: uplo, trans, diag
        integer,     intent(in)    :: n
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(:,:), x(:)
        real(ep),    intent(inout) :: y(:)
        real(ep), allocatable :: aA(:,:), aX(:)
        integer :: i, j
        allocate(aA(n, n), aX(n))
        aA = 0.0_ep
        do j = 1, n; aX(j) = abs(real(x(j), ep)) + abs(aimag(x(j))); end do
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    aA(i, j) = abs(real(A(i, j), ep)) + abs(aimag(A(i, j)))
                end do
            end do
        end if
        if (diag == 'U' .or. diag == 'u') then
            do i = 1, n; aA(i, i) = 1.0_ep; end do
        end if
        if (trans == 'N' .or. trans == 'n') then
            y(1:n) = abs(alpha) * matmul(aA, aX) + abs(beta * y(1:n))
        else
            y(1:n) = abs(alpha) * matmul(transpose(aA), aX) + abs(beta * y(1:n))
        end if
        deallocate(aA, aX)
    end subroutine ref_zatrmv

end module ptzblas_ref_quad_aux
