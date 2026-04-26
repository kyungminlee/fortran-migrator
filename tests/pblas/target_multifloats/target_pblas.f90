! Per-target PBLAS wrapper for the multifloats (ddpblas) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers split each
! quad value into the high (double approx) + low (double remainder)
! representation that TYPE(real64x2) expects, call the migrated
! dd*/zz*-prefix PBLAS routine, then recombine the two limbs back to
! quad for the gather/compare step. Distribution descriptors pass
! through unchanged — precision conversion is purely elementwise.
!
! Migrated prefixes for multifloats:
!   pdgemm → ptgemm, pddot → ptdot, …, pzgemm → pvgemm, pzdotc → pvdotc

module target_pblas
    use prec_kinds,  only: ep, dp
    use multifloats, only: real64x2, cmplx64x2
    implicit none
    private

    public :: target_name, target_eps
    public :: target_pddot, target_pdnrm2, target_pdasum
    public :: target_pdscal, target_pdaxpy, target_pdcopy
    public :: target_pzdotc, target_pzaxpy
    public :: target_pdgemv, target_pdger, target_pdsymv, target_pdtrsv
    public :: target_pzgemv
    public :: target_pdgemm, target_pdsymm, target_pdsyrk
    public :: target_pdtrmm, target_pdtrsm
    public :: target_pzgemm, target_pzherk

    character(len=*), parameter :: target_name = 'multifloats'
    ! Double-double effective epsilon ≈ 2^-104 ≈ 4.93e-32.
    real(ep),         parameter :: target_eps  = real(2.0_dp**(-104), ep)

    interface
        ! ── Level 1 — real (DD) ──────────────────────────────────────
        subroutine ptdot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx, iy, jy, incy
            integer, intent(in) :: descx(9), descy(9)
            type(real64x2), intent(in)  :: x(*), y(*)
            type(real64x2), intent(out) :: dot
        end subroutine
        subroutine ptnrm2(n, norm2, x, ix, jx, descx, incx)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx
            integer, intent(in) :: descx(9)
            type(real64x2), intent(in)  :: x(*)
            type(real64x2), intent(out) :: norm2
        end subroutine
        subroutine ptasum(n, asum, x, ix, jx, descx, incx)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx
            integer, intent(in) :: descx(9)
            type(real64x2), intent(in)  :: x(*)
            type(real64x2), intent(out) :: asum
        end subroutine
        subroutine ptscal(n, alpha, x, ix, jx, descx, incx)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx
            integer, intent(in) :: descx(9)
            type(real64x2), intent(in)    :: alpha
            type(real64x2), intent(inout) :: x(*)
        end subroutine
        subroutine ptaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx, iy, jy, incy
            integer, intent(in) :: descx(9), descy(9)
            type(real64x2), intent(in)    :: alpha, x(*)
            type(real64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ptcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: real64x2
            integer, intent(in) :: n, ix, jx, incx, iy, jy, incy
            integer, intent(in) :: descx(9), descy(9)
            type(real64x2), intent(in)  :: x(*)
            type(real64x2), intent(out) :: y(*)
        end subroutine

        ! ── Level 1 — complex (ZZ) ───────────────────────────────────
        subroutine pvdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: cmplx64x2
            integer, intent(in) :: n, ix, jx, incx, iy, jy, incy
            integer, intent(in) :: descx(9), descy(9)
            type(cmplx64x2), intent(in)  :: x(*), y(*)
            type(cmplx64x2), intent(out) :: dot
        end subroutine
        subroutine pvaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            import :: cmplx64x2
            integer, intent(in) :: n, ix, jx, incx, iy, jy, incy
            integer, intent(in) :: descx(9), descy(9)
            type(cmplx64x2), intent(in)    :: alpha, x(*)
            type(cmplx64x2), intent(inout) :: y(*)
        end subroutine

        ! ── Level 2 — real (DD) ──────────────────────────────────────
        subroutine ptgemv(trans, m, n, alpha, A, ia, ja, desca, &
                           x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: real64x2
            character, intent(in) :: trans
            integer,   intent(in) :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in) :: desca(9), descx(9), descy(9)
            type(real64x2), intent(in)    :: alpha, beta, A(*), x(*)
            type(real64x2), intent(inout) :: y(*)
        end subroutine
        subroutine ptger(m, n, alpha, x, ix, jx, descx, incx, &
                          y, iy, jy, descy, incy, A, ia, ja, desca)
            import :: real64x2
            integer, intent(in) :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
            integer, intent(in) :: desca(9), descx(9), descy(9)
            type(real64x2), intent(in)    :: alpha, x(*), y(*)
            type(real64x2), intent(inout) :: A(*)
        end subroutine
        subroutine ptsymv(uplo, n, alpha, A, ia, ja, desca, &
                           x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: real64x2
            character, intent(in) :: uplo
            integer,   intent(in) :: n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in) :: desca(9), descx(9), descy(9)
            type(real64x2), intent(in)    :: alpha, beta, A(*), x(*)
            type(real64x2), intent(inout) :: y(*)
        end subroutine
        subroutine pttrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                           x, ix, jx, descx, incx)
            import :: real64x2
            character, intent(in) :: uplo, trans, diag
            integer,   intent(in) :: n, ia, ja, ix, jx, incx
            integer,   intent(in) :: desca(9), descx(9)
            type(real64x2), intent(in)    :: A(*)
            type(real64x2), intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex (ZZ) ───────────────────────────────────
        subroutine pvgemv(trans, m, n, alpha, A, ia, ja, desca, &
                           x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            import :: cmplx64x2
            character, intent(in) :: trans
            integer,   intent(in) :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in) :: desca(9), descx(9), descy(9)
            type(cmplx64x2), intent(in)    :: alpha, beta, A(*), x(*)
            type(cmplx64x2), intent(inout) :: y(*)
        end subroutine

        ! ── Level 3 — real (DD) ──────────────────────────────────────
        subroutine ptgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                           B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: real64x2
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,   intent(in) :: desca(9), descb(9), descc(9)
            type(real64x2), intent(in)    :: alpha, beta, A(*), B(*)
            type(real64x2), intent(inout) :: C(*)
        end subroutine
        subroutine ptsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                           B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: real64x2
            character, intent(in) :: side, uplo
            integer,   intent(in) :: m, n, ia, ja, ib, jb, ic, jc
            integer,   intent(in) :: desca(9), descb(9), descc(9)
            type(real64x2), intent(in)    :: alpha, beta, A(*), B(*)
            type(real64x2), intent(inout) :: C(*)
        end subroutine
        subroutine ptsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                           beta, C, ic, jc, descc)
            import :: real64x2
            character, intent(in) :: uplo, trans
            integer,   intent(in) :: n, k, ia, ja, ic, jc
            integer,   intent(in) :: desca(9), descc(9)
            type(real64x2), intent(in)    :: alpha, beta, A(*)
            type(real64x2), intent(inout) :: C(*)
        end subroutine
        subroutine pttrmm(side, uplo, trans, diag, m, n, alpha, &
                           A, ia, ja, desca, B, ib, jb, descb)
            import :: real64x2
            character, intent(in) :: side, uplo, trans, diag
            integer,   intent(in) :: m, n, ia, ja, ib, jb
            integer,   intent(in) :: desca(9), descb(9)
            type(real64x2), intent(in)    :: alpha, A(*)
            type(real64x2), intent(inout) :: B(*)
        end subroutine
        subroutine pttrsm(side, uplo, trans, diag, m, n, alpha, &
                           A, ia, ja, desca, B, ib, jb, descb)
            import :: real64x2
            character, intent(in) :: side, uplo, trans, diag
            integer,   intent(in) :: m, n, ia, ja, ib, jb
            integer,   intent(in) :: desca(9), descb(9)
            type(real64x2), intent(in)    :: alpha, A(*)
            type(real64x2), intent(inout) :: B(*)
        end subroutine

        ! ── Level 3 — complex (ZZ) ───────────────────────────────────
        subroutine pvgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                           B, ib, jb, descb, beta, C, ic, jc, descc)
            import :: cmplx64x2
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,   intent(in) :: desca(9), descb(9), descc(9)
            type(cmplx64x2), intent(in)    :: alpha, beta, A(*), B(*)
            type(cmplx64x2), intent(inout) :: C(*)
        end subroutine
        subroutine pvherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                           beta, C, ic, jc, descc)
            import :: real64x2, cmplx64x2
            character, intent(in) :: uplo, trans
            integer,   intent(in) :: n, k, ia, ja, ic, jc
            integer,   intent(in) :: desca(9), descc(9)
            type(real64x2),   intent(in)    :: alpha, beta
            type(cmplx64x2), intent(in)    :: A(*)
            type(cmplx64x2), intent(inout) :: C(*)
        end subroutine
    end interface

contains

    ! ── Quad ↔ double-double splitting (same pattern as tests/blas) ──
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

    ! ── Level 1 — real ───────────────────────────────────────────────
    subroutine target_pddot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*), y(*)
        real(ep), intent(out) :: dot
        integer :: nx, ny
        type(real64x2), allocatable :: xt(:), yt(:)
        type(real64x2) :: dot_t
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ptdot(n, dot_t, xt, ix, jx, descx, incx, yt, iy, jy, descy, incy)
        dot = dd2q(dot_t)
    end subroutine

    subroutine target_pdnrm2(n, norm2, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: norm2
        integer :: nx
        type(real64x2), allocatable :: xt(:)
        type(real64x2) :: norm_t
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        call ptnrm2(n, norm_t, xt, ix, jx, descx, incx)
        norm2 = dd2q(norm_t)
    end subroutine

    subroutine target_pdasum(n, asum, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: asum
        integer :: nx
        type(real64x2), allocatable :: xt(:)
        type(real64x2) :: asum_t
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        call ptasum(n, asum_t, xt, ix, jx, descx, incx)
        asum = dd2q(asum_t)
    end subroutine

    subroutine target_pdscal(n, alpha, x, ix, jx, descx, incx)
        integer,  intent(in)    :: n, ix, jx, incx
        integer,  intent(in)    :: descx(9)
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        integer :: nx
        type(real64x2), allocatable :: xt(:)
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = q2dd(x(1:nx))
        call ptscal(n, q2dd(alpha), xt, ix, jx, descx, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    subroutine target_pdaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)    :: descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(real64x2), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ptaxpy(n, q2dd(alpha), xt, ix, jx, descx, incx, &
                     yt, iy, jy, descy, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_pdcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        integer :: nx, ny
        type(real64x2), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = q2dd(x(1:nx))
        call ptcopy(n, xt, ix, jx, descx, incx, yt, iy, jy, descy, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    subroutine target_pzdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)  :: descx(9), descy(9)
        complex(ep), intent(in)  :: x(*), y(*)
        complex(ep), intent(out) :: dot
        integer :: nx, ny
        type(cmplx64x2), allocatable :: xt(:), yt(:)
        type(cmplx64x2) :: dot_t
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call pvdotc(n, dot_t, xt, ix, jx, descx, incx, yt, iy, jy, descy, incy)
        dot = zz2q(dot_t)
    end subroutine

    subroutine target_pzaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)    :: descx(9), descy(9)
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        type(cmplx64x2), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call pvaxpy(n, q2zz(alpha), xt, ix, jx, descx, incx, &
                     yt, iy, jy, descy, incy)
        y(1:ny) = zz2q(yt)
    end subroutine

    ! ── Level 2 — real ───────────────────────────────────────────────
    subroutine target_pdgemv(trans, m, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,   intent(in)    :: desca(9), descx(9), descy(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: na, nx, ny
        type(real64x2), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = q2dd(A(1:na))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ptgemv(trans, m, n, q2dd(alpha), At, ia, ja, desca, &
                     xt, ix, jx, descx, incx, q2dd(beta), &
                     yt, iy, jy, descy, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_pdger(m, n, alpha, x, ix, jx, descx, incx, &
                            y, iy, jy, descy, incy, A, ia, ja, desca)
        integer,  intent(in)    :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
        integer,  intent(in)    :: desca(9), descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(*)
        integer :: na, nx, ny
        type(real64x2), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = q2dd(A(1:na))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ptger(m, n, q2dd(alpha), xt, ix, jx, descx, incx, &
                    yt, iy, jy, descy, incy, At, ia, ja, desca)
        A(1:na) = dd2q(At)
    end subroutine

    subroutine target_pdsymv(uplo, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,   intent(in)    :: desca(9), descx(9), descy(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: na, nx, ny
        type(real64x2), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = q2dd(A(1:na))
        xt = q2dd(x(1:nx)); yt = q2dd(y(1:ny))
        call ptsymv(uplo, n, q2dd(alpha), At, ia, ja, desca, &
                     xt, ix, jx, descx, incx, q2dd(beta), &
                     yt, iy, jy, descy, incy)
        y(1:ny) = dd2q(yt)
    end subroutine

    subroutine target_pdtrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                             x, ix, jx, descx, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx
        integer,   intent(in)    :: desca(9), descx(9)
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: x(*)
        integer :: na, nx
        type(real64x2), allocatable :: At(:), xt(:)
        na = desca(9) * desca(4); nx = descx(9) * descx(4)
        allocate(At(na), xt(nx))
        At = q2dd(A(1:na)); xt = q2dd(x(1:nx))
        call pttrsv(uplo, trans, diag, n, At, ia, ja, desca, &
                     xt, ix, jx, descx, incx)
        x(1:nx) = dd2q(xt)
    end subroutine

    ! ── Level 2 — complex ────────────────────────────────────────────
    subroutine target_pzgemv(trans, m, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,     intent(in)    :: desca(9), descx(9), descy(9)
        complex(ep), intent(in)    :: alpha, beta, A(*), x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: na, nx, ny
        type(cmplx64x2), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4); ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = q2zz(A(1:na))
        xt = q2zz(x(1:nx)); yt = q2zz(y(1:ny))
        call pvgemv(trans, m, n, q2zz(alpha), At, ia, ja, desca, &
                     xt, ix, jx, descx, incx, q2zz(beta), &
                     yt, iy, jy, descy, incy)
        y(1:ny) = zz2q(yt)
    end subroutine

    ! ── Level 3 — real ───────────────────────────────────────────────
    subroutine target_pdgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
        integer,   intent(in)    :: desca(9), descb(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
        real(ep),  intent(inout) :: C(*)
        integer :: na, nb, nc
        type(real64x2), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4); nb = descb(9) * descb(4); nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb)); Ct = q2dd(C(1:nc))
        call ptgemm(transa, transb, m, n, k, q2dd(alpha), &
                     At, ia, ja, desca, Bt, ib, jb, descb, &
                     q2dd(beta), Ct, ic, jc, descc)
        C(1:nc) = dd2q(Ct)
    end subroutine

    subroutine target_pdsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, ia, ja, ib, jb, ic, jc
        integer,   intent(in)    :: desca(9), descb(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
        real(ep),  intent(inout) :: C(*)
        integer :: na, nb, nc
        type(real64x2), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4); nb = descb(9) * descb(4); nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb)); Ct = q2dd(C(1:nc))
        call ptsymm(side, uplo, m, n, q2dd(alpha), &
                     At, ia, ja, desca, Bt, ib, jb, descb, &
                     q2dd(beta), Ct, ic, jc, descc)
        C(1:nc) = dd2q(Ct)
    end subroutine

    subroutine target_pdsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                             beta, C, ic, jc, descc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, ia, ja, ic, jc
        integer,   intent(in)    :: desca(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*)
        real(ep),  intent(inout) :: C(*)
        integer :: na, nc
        type(real64x2), allocatable :: At(:), Ct(:)
        na = desca(9) * desca(4); nc = descc(9) * descc(4)
        allocate(At(na), Ct(nc))
        At = q2dd(A(1:na)); Ct = q2dd(C(1:nc))
        call ptsyrk(uplo, trans, n, k, q2dd(alpha), &
                     At, ia, ja, desca, q2dd(beta), Ct, ic, jc, descc)
        C(1:nc) = dd2q(Ct)
    end subroutine

    subroutine target_pdtrmm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        type(real64x2), allocatable :: At(:), Bt(:)
        na = desca(9) * desca(4); nb = descb(9) * descb(4)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb))
        call pttrmm(side, uplo, trans, diag, m, n, q2dd(alpha), &
                     At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = dd2q(Bt)
    end subroutine

    subroutine target_pdtrsm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        type(real64x2), allocatable :: At(:), Bt(:)
        na = desca(9) * desca(4); nb = descb(9) * descb(4)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb))
        call pttrsm(side, uplo, trans, diag, m, n, q2dd(alpha), &
                     At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = dd2q(Bt)
    end subroutine

    ! ── Level 3 — complex ────────────────────────────────────────────
    subroutine target_pzgemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
        integer,     intent(in)    :: desca(9), descb(9), descc(9)
        complex(ep), intent(in)    :: alpha, beta, A(*), B(*)
        complex(ep), intent(inout) :: C(*)
        integer :: na, nb, nc
        type(cmplx64x2), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4); nb = descb(9) * descb(4); nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = q2zz(A(1:na)); Bt = q2zz(B(1:nb)); Ct = q2zz(C(1:nc))
        call pvgemm(transa, transb, m, n, k, q2zz(alpha), &
                     At, ia, ja, desca, Bt, ib, jb, descb, &
                     q2zz(beta), Ct, ic, jc, descc)
        C(1:nc) = zz2q(Ct)
    end subroutine

    subroutine target_pzherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                             beta, C, ic, jc, descc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, ia, ja, ic, jc
        integer,     intent(in)    :: desca(9), descc(9)
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(*)
        complex(ep), intent(inout) :: C(*)
        integer :: na, nc
        type(cmplx64x2), allocatable :: At(:), Ct(:)
        na = desca(9) * desca(4); nc = descc(9) * descc(4)
        allocate(At(na), Ct(nc))
        At = q2zz(A(1:na)); Ct = q2zz(C(1:nc))
        call pvherk(uplo, trans, n, k, q2dd(alpha), &
                     At, ia, ja, desca, q2dd(beta), Ct, ic, jc, descc)
        C(1:nc) = zz2q(Ct)
    end subroutine

end module target_pblas
