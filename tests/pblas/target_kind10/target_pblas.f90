! Per-target PBLAS wrapper for the kind10 (epblas) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers cast local
! quad inputs down to REAL(KIND=10) (the x87 80-bit extended type) for
! the call into the migrated e/y-prefix PBLAS routines, then cast
! results back up to quad. Distribution bookkeeping (descriptors,
! IA/JA indices, INCX) passes through untouched — descriptors carry
! distribution info, not element precision.
!
! Migrated prefixes for kind10:
!   pdgemm → pegemm, pddot → pedot, ..., pzgemm → pygemm, pzdotc → pydotc

module target_pblas
    use prec_kinds, only: ep
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

    integer,          parameter :: tk = 10
    character(len=*), parameter :: target_name = 'kind10'
    real(ep),         parameter :: target_eps  = real(epsilon(1.0_10), ep)

    interface
        ! ── Level 1 — real (E) ───────────────────────────────────────
        subroutine pedot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)  :: descx(9), descy(9)
            real(10), intent(in)  :: x(*), y(*)
            real(10), intent(out) :: dot
        end subroutine
        subroutine penrm2(n, norm2, x, ix, jx, descx, incx)
            integer,  intent(in)  :: n, ix, jx, incx
            integer,  intent(in)  :: descx(9)
            real(10), intent(in)  :: x(*)
            real(10), intent(out) :: norm2
        end subroutine
        subroutine peasum(n, asum, x, ix, jx, descx, incx)
            integer,  intent(in)  :: n, ix, jx, incx
            integer,  intent(in)  :: descx(9)
            real(10), intent(in)  :: x(*)
            real(10), intent(out) :: asum
        end subroutine
        subroutine pescal(n, alpha, x, ix, jx, descx, incx)
            integer,  intent(in)    :: n, ix, jx, incx
            integer,  intent(in)    :: descx(9)
            real(10), intent(in)    :: alpha
            real(10), intent(inout) :: x(*)
        end subroutine
        subroutine peaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)    :: descx(9), descy(9)
            real(10), intent(in)    :: alpha, x(*)
            real(10), intent(inout) :: y(*)
        end subroutine
        subroutine pecopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
            integer,  intent(in)  :: descx(9), descy(9)
            real(10), intent(in)  :: x(*)
            real(10), intent(out) :: y(*)
        end subroutine

        ! ── Level 1 — complex (Y) ────────────────────────────────────
        subroutine pydotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            integer,     intent(in)  :: n, ix, jx, incx, iy, jy, incy
            integer,     intent(in)  :: descx(9), descy(9)
            complex(10), intent(in)  :: x(*), y(*)
            complex(10), intent(out) :: dot
        end subroutine
        subroutine pyaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
            integer,     intent(in)    :: n, ix, jx, incx, iy, jy, incy
            integer,     intent(in)    :: descx(9), descy(9)
            complex(10), intent(in)    :: alpha, x(*)
            complex(10), intent(inout) :: y(*)
        end subroutine

        ! ── Level 2 — real (E) ───────────────────────────────────────
        subroutine pegemv(trans, m, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in)    :: desca(9), descx(9), descy(9)
            real(10),  intent(in)    :: alpha, beta, A(*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine peger(m, n, alpha, x, ix, jx, descx, incx, &
                         y, iy, jy, descy, incy, A, ia, ja, desca)
            integer,  intent(in)    :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
            integer,  intent(in)    :: desca(9), descx(9), descy(9)
            real(10), intent(in)    :: alpha, x(*), y(*)
            real(10), intent(inout) :: A(*)
        end subroutine
        subroutine pesymv(uplo, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,   intent(in)    :: desca(9), descx(9), descy(9)
            real(10),  intent(in)    :: alpha, beta, A(*), x(*)
            real(10),  intent(inout) :: y(*)
        end subroutine
        subroutine petrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                          x, ix, jx, descx, incx)
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, ia, ja, ix, jx, incx
            integer,   intent(in)    :: desca(9), descx(9)
            real(10),  intent(in)    :: A(*)
            real(10),  intent(inout) :: x(*)
        end subroutine

        ! ── Level 2 — complex (Y) ────────────────────────────────────
        subroutine pygemv(trans, m, n, alpha, A, ia, ja, desca, &
                          x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, ia, ja, ix, jx, incx, iy, jy, incy
            integer,     intent(in)    :: desca(9), descx(9), descy(9)
            complex(10), intent(in)    :: alpha, beta, A(*), x(*)
            complex(10), intent(inout) :: y(*)
        end subroutine

        ! ── Level 3 — real (E) ───────────────────────────────────────
        subroutine pegemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,   intent(in)    :: desca(9), descb(9), descc(9)
            real(10),  intent(in)    :: alpha, beta, A(*), B(*)
            real(10),  intent(inout) :: C(*)
        end subroutine
        subroutine pesymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, ia, ja, ib, jb, ic, jc
            integer,   intent(in)    :: desca(9), descb(9), descc(9)
            real(10),  intent(in)    :: alpha, beta, A(*), B(*)
            real(10),  intent(inout) :: C(*)
        end subroutine
        subroutine pesyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                          beta, C, ic, jc, descc)
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, ia, ja, ic, jc
            integer,   intent(in)    :: desca(9), descc(9)
            real(10),  intent(in)    :: alpha, beta, A(*)
            real(10),  intent(inout) :: C(*)
        end subroutine
        subroutine petrmm(side, uplo, trans, diag, m, n, alpha, &
                          A, ia, ja, desca, B, ib, jb, descb)
            character, intent(in)    :: side, uplo, trans, diag
            integer,   intent(in)    :: m, n, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            real(10),  intent(in)    :: alpha, A(*)
            real(10),  intent(inout) :: B(*)
        end subroutine
        subroutine petrsm(side, uplo, trans, diag, m, n, alpha, &
                          A, ia, ja, desca, B, ib, jb, descb)
            character, intent(in)    :: side, uplo, trans, diag
            integer,   intent(in)    :: m, n, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            real(10),  intent(in)    :: alpha, A(*)
            real(10),  intent(inout) :: B(*)
        end subroutine

        ! ── Level 3 — complex (Y) ────────────────────────────────────
        subroutine pygemm(transa, transb, m, n, k, alpha, A, ia, ja, desca, &
                          B, ib, jb, descb, beta, C, ic, jc, descc)
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc
            integer,     intent(in)    :: desca(9), descb(9), descc(9)
            complex(10), intent(in)    :: alpha, beta, A(*), B(*)
            complex(10), intent(inout) :: C(*)
        end subroutine
        subroutine pyherk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                          beta, C, ic, jc, descc)
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, ia, ja, ic, jc
            integer,     intent(in)    :: desca(9), descc(9)
            real(10),    intent(in)    :: alpha, beta
            complex(10), intent(in)    :: A(*)
            complex(10), intent(inout) :: C(*)
        end subroutine
    end interface

contains

    ! ── Level 1 — real ───────────────────────────────────────────────
    subroutine target_pddot(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*), y(*)
        real(ep), intent(out) :: dot
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        real(tk) :: dot_t
        nx = descx(9) * descx(4)   ! lld * global_n (upper bound on array size)
        ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call pedot(n, dot_t, xt, ix, jx, descx, incx, &
                   yt, iy, jy, descy, incy)
        dot = real(dot_t, ep)
    end subroutine

    subroutine target_pdnrm2(n, norm2, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: norm2
        integer :: nx
        real(tk), allocatable :: xt(:)
        real(tk) :: norm_t
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        call penrm2(n, norm_t, xt, ix, jx, descx, incx)
        norm2 = real(norm_t, ep)
    end subroutine

    subroutine target_pdasum(n, asum, x, ix, jx, descx, incx)
        integer,  intent(in)  :: n, ix, jx, incx
        integer,  intent(in)  :: descx(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: asum
        integer :: nx
        real(tk), allocatable :: xt(:)
        real(tk) :: asum_t
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        call peasum(n, asum_t, xt, ix, jx, descx, incx)
        asum = real(asum_t, ep)
    end subroutine

    subroutine target_pdscal(n, alpha, x, ix, jx, descx, incx)
        integer,  intent(in)    :: n, ix, jx, incx
        integer,  intent(in)    :: descx(9)
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        integer :: nx
        real(tk), allocatable :: xt(:)
        nx = descx(9) * descx(4)
        allocate(xt(nx))
        xt = real(x(1:nx), tk)
        call pescal(n, real(alpha, tk), xt, ix, jx, descx, incx)
        x(1:nx) = real(xt, ep)
    end subroutine

    subroutine target_pdaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)    :: descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call peaxpy(n, real(alpha, tk), xt, ix, jx, descx, incx, &
                    yt, iy, jy, descy, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_pdcopy(n, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,  intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,  intent(in)  :: descx(9), descy(9)
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        integer :: nx, ny
        real(tk), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = real(x(1:nx), tk)
        call pecopy(n, xt, ix, jx, descx, incx, yt, iy, jy, descy, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    ! ── Level 1 — complex ────────────────────────────────────────────
    subroutine target_pzdotc(n, dot, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)  :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)  :: descx(9), descy(9)
        complex(ep), intent(in)  :: x(*), y(*)
        complex(ep), intent(out) :: dot
        integer :: nx, ny
        complex(tk), allocatable :: xt(:), yt(:)
        complex(tk) :: dot_t
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call pydotc(n, dot_t, xt, ix, jx, descx, incx, yt, iy, jy, descy, incy)
        dot = cmplx(dot_t, kind=ep)
    end subroutine

    subroutine target_pzaxpy(n, alpha, x, ix, jx, descx, incx, y, iy, jy, descy, incy)
        integer,     intent(in)    :: n, ix, jx, incx, iy, jy, incy
        integer,     intent(in)    :: descx(9), descy(9)
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: nx, ny
        complex(tk), allocatable :: xt(:), yt(:)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(xt(nx), yt(ny))
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call pyaxpy(n, cmplx(alpha, kind=tk), xt, ix, jx, descx, incx, &
                    yt, iy, jy, descy, incy)
        y(1:ny) = cmplx(yt, kind=ep)
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
        real(tk), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = real(A(1:na), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call pegemv(trans, m, n, real(alpha, tk), At, ia, ja, desca, &
                    xt, ix, jx, descx, incx, real(beta, tk), &
                    yt, iy, jy, descy, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_pdger(m, n, alpha, x, ix, jx, descx, incx, &
                            y, iy, jy, descy, incy, A, ia, ja, desca)
        integer,  intent(in)    :: m, n, ix, jx, incx, iy, jy, incy, ia, ja
        integer,  intent(in)    :: desca(9), descx(9), descy(9)
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(*)
        integer :: na, nx, ny
        real(tk), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = real(A(1:na), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call peger(m, n, real(alpha, tk), xt, ix, jx, descx, incx, &
                   yt, iy, jy, descy, incy, At, ia, ja, desca)
        A(1:na) = real(At, ep)
    end subroutine

    subroutine target_pdsymv(uplo, n, alpha, A, ia, ja, desca, &
                             x, ix, jx, descx, incx, beta, y, iy, jy, descy, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx, iy, jy, incy
        integer,   intent(in)    :: desca(9), descx(9), descy(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), x(*)
        real(ep),  intent(inout) :: y(*)
        integer :: na, nx, ny
        real(tk), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = real(A(1:na), tk)
        xt = real(x(1:nx), tk); yt = real(y(1:ny), tk)
        call pesymv(uplo, n, real(alpha, tk), At, ia, ja, desca, &
                    xt, ix, jx, descx, incx, real(beta, tk), &
                    yt, iy, jy, descy, incy)
        y(1:ny) = real(yt, ep)
    end subroutine

    subroutine target_pdtrsv(uplo, trans, diag, n, A, ia, ja, desca, &
                             x, ix, jx, descx, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, ia, ja, ix, jx, incx
        integer,   intent(in)    :: desca(9), descx(9)
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: x(*)
        integer :: na, nx
        real(tk), allocatable :: At(:), xt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4)
        allocate(At(na), xt(nx))
        At = real(A(1:na), tk)
        xt = real(x(1:nx), tk)
        call petrsv(uplo, trans, diag, n, At, ia, ja, desca, &
                    xt, ix, jx, descx, incx)
        x(1:nx) = real(xt, ep)
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
        complex(tk), allocatable :: At(:), xt(:), yt(:)
        na = desca(9) * desca(4)
        nx = descx(9) * descx(4)
        ny = descy(9) * descy(4)
        allocate(At(na), xt(nx), yt(ny))
        At = cmplx(A(1:na), kind=tk)
        xt = cmplx(x(1:nx), kind=tk); yt = cmplx(y(1:ny), kind=tk)
        call pygemv(trans, m, n, cmplx(alpha, kind=tk), At, ia, ja, desca, &
                    xt, ix, jx, descx, incx, cmplx(beta, kind=tk), &
                    yt, iy, jy, descy, incy)
        y(1:ny) = cmplx(yt, kind=ep)
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
        real(tk), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4)
        nb = descb(9) * descb(4)
        nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk); Ct = real(C(1:nc), tk)
        call pegemm(transa, transb, m, n, k, real(alpha, tk), &
                    At, ia, ja, desca, Bt, ib, jb, descb, &
                    real(beta, tk), Ct, ic, jc, descc)
        C(1:nc) = real(Ct, ep)
    end subroutine

    subroutine target_pdsymm(side, uplo, m, n, alpha, A, ia, ja, desca, &
                             B, ib, jb, descb, beta, C, ic, jc, descc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, ia, ja, ib, jb, ic, jc
        integer,   intent(in)    :: desca(9), descb(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*), B(*)
        real(ep),  intent(inout) :: C(*)
        integer :: na, nb, nc
        real(tk), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4)
        nb = descb(9) * descb(4)
        nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk); Ct = real(C(1:nc), tk)
        call pesymm(side, uplo, m, n, real(alpha, tk), &
                    At, ia, ja, desca, Bt, ib, jb, descb, &
                    real(beta, tk), Ct, ic, jc, descc)
        C(1:nc) = real(Ct, ep)
    end subroutine

    subroutine target_pdsyrk(uplo, trans, n, k, alpha, A, ia, ja, desca, &
                             beta, C, ic, jc, descc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, ia, ja, ic, jc
        integer,   intent(in)    :: desca(9), descc(9)
        real(ep),  intent(in)    :: alpha, beta, A(*)
        real(ep),  intent(inout) :: C(*)
        integer :: na, nc
        real(tk), allocatable :: At(:), Ct(:)
        na = desca(9) * desca(4)
        nc = descc(9) * descc(4)
        allocate(At(na), Ct(nc))
        At = real(A(1:na), tk); Ct = real(C(1:nc), tk)
        call pesyrk(uplo, trans, n, k, real(alpha, tk), &
                    At, ia, ja, desca, real(beta, tk), Ct, ic, jc, descc)
        C(1:nc) = real(Ct, ep)
    end subroutine

    subroutine target_pdtrmm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = desca(9) * desca(4)
        nb = descb(9) * descb(4)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk)
        call petrmm(side, uplo, trans, diag, m, n, real(alpha, tk), &
                    At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = real(Bt, ep)
    end subroutine

    subroutine target_pdtrsm(side, uplo, trans, diag, m, n, alpha, &
                             A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)    :: side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        real(ep),  intent(in)    :: alpha, A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = desca(9) * desca(4)
        nb = descb(9) * descb(4)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk)
        call petrsm(side, uplo, trans, diag, m, n, real(alpha, tk), &
                    At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = real(Bt, ep)
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
        complex(tk), allocatable :: At(:), Bt(:), Ct(:)
        na = desca(9) * desca(4)
        nb = descb(9) * descb(4)
        nc = descc(9) * descc(4)
        allocate(At(na), Bt(nb), Ct(nc))
        At = cmplx(A(1:na), kind=tk)
        Bt = cmplx(B(1:nb), kind=tk)
        Ct = cmplx(C(1:nc), kind=tk)
        call pygemm(transa, transb, m, n, k, cmplx(alpha, kind=tk), &
                    At, ia, ja, desca, Bt, ib, jb, descb, &
                    cmplx(beta, kind=tk), Ct, ic, jc, descc)
        C(1:nc) = cmplx(Ct, kind=ep)
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
        complex(tk), allocatable :: At(:), Ct(:)
        na = desca(9) * desca(4)
        nc = descc(9) * descc(4)
        allocate(At(na), Ct(nc))
        At = cmplx(A(1:na), kind=tk)
        Ct = cmplx(C(1:nc), kind=tk)
        call pyherk(uplo, trans, n, k, real(alpha, tk), &
                    At, ia, ja, desca, real(beta, tk), Ct, ic, jc, descc)
        C(1:nc) = cmplx(Ct, kind=ep)
    end subroutine

end module target_pblas
