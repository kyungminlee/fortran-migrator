! Per-target ScaLAPACK wrapper for the kind10 (escalapack) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers cast quad
! inputs down to REAL(KIND=10) (x87 80-bit extended), call the migrated
! e/y-prefix ScaLAPACK routines, then cast outputs back up to quad.
! Workspace arrays are allocated locally at REAL(KIND=10); for the
! lwork=-1 query, the optimal value returned in work_t(1) is cast back
! into the caller's quad-precision work(1).

module target_scalapack
    use prec_kinds, only: ep
    use pblas_grid, only: numroc_local
    implicit none
    private

    interface
        subroutine blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
            integer, intent(in)  :: icontxt
            integer, intent(out) :: nprow, npcol, myrow, mycol
        end subroutine
    end interface

    public :: target_name, target_eps
    public :: target_pdgesv, target_pdgetrf, target_pdgetrs
    public :: target_pdpotrf, target_pdpotrs
    public :: target_pdgeqrf
    public :: target_pdsyev
    public :: target_pdgesvd
    public :: target_pdlange
    public :: target_pdlacpy, target_pdlaset
    public :: target_pzgesv, target_pzgeqrf, target_pzheev

    integer,          parameter :: tk = 10
    character(len=*), parameter :: target_name = 'kind10'
    real(ep),         parameter :: target_eps  = real(epsilon(1.0_10), ep)

    interface
        subroutine pegesv(n, nrhs, A, ia, ja, desca, ipiv, &
                          B, ib, jb, descb, info)
            integer,  intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,  intent(in)    :: desca(9), descb(9)
            integer,  intent(out)   :: ipiv(*), info
            real(10), intent(inout) :: A(*), B(*)
        end subroutine
        subroutine pegetrf(m, n, A, ia, ja, desca, ipiv, info)
            integer,  intent(in)    :: m, n, ia, ja
            integer,  intent(in)    :: desca(9)
            integer,  intent(out)   :: ipiv(*), info
            real(10), intent(inout) :: A(*)
        end subroutine
        subroutine pegetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                           B, ib, jb, descb, info)
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9), ipiv(*)
            integer,   intent(out)   :: info
            real(10),  intent(in)    :: A(*)
            real(10),  intent(inout) :: B(*)
        end subroutine
        subroutine pepotrf(uplo, n, A, ia, ja, desca, info)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ia, ja
            integer,   intent(in)    :: desca(9)
            integer,   intent(out)   :: info
            real(10),  intent(inout) :: A(*)
        end subroutine
        subroutine pepotrs(uplo, n, nrhs, A, ia, ja, desca, &
                           B, ib, jb, descb, info)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            integer,   intent(out)   :: info
            real(10),  intent(in)    :: A(*)
            real(10),  intent(inout) :: B(*)
        end subroutine
        subroutine pegeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            integer,  intent(in)    :: m, n, ia, ja, lwork
            integer,  intent(in)    :: desca(9)
            integer,  intent(out)   :: info
            real(10), intent(inout) :: A(*)
            real(10), intent(out)   :: tau(*), work(*)
        end subroutine
        subroutine pesyev(jobz, uplo, n, A, ia, ja, desca, w, &
                          Z, iz, jz, descz, work, lwork, info)
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ia, ja, iz, jz, lwork
            integer,   intent(in)    :: desca(9), descz(9)
            integer,   intent(out)   :: info
            real(10),  intent(inout) :: A(*)
            real(10),  intent(out)   :: w(*), Z(*), work(*)
        end subroutine
        subroutine pegesvd(jobu, jobvt, m, n, A, ia, ja, desca, s, &
                           U, iu, ju, descu, VT, ivt, jvt, descvt, &
                           work, lwork, info)
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, ia, ja, iu, ju, ivt, jvt, lwork
            integer,   intent(in)    :: desca(9), descu(9), descvt(9)
            integer,   intent(out)   :: info
            real(10),  intent(inout) :: A(*)
            real(10),  intent(out)   :: s(*), U(*), VT(*), work(*)
        end subroutine
        function pelange(norm, m, n, A, ia, ja, desca, work) result(r)
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, ia, ja
            integer,   intent(in) :: desca(9)
            real(10),  intent(in) :: A(*)
            real(10)              :: work(*)
            real(10) :: r
        end function
        subroutine pelacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, ia, ja, ib, jb
            integer,   intent(in)  :: desca(9), descb(9)
            real(10),  intent(in)  :: A(*)
            real(10),  intent(out) :: B(*)
        end subroutine
        subroutine pelaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, ia, ja
            integer,   intent(in)    :: desca(9)
            real(10),  intent(in)    :: alpha, beta
            real(10),  intent(inout) :: A(*)
        end subroutine

        subroutine pygesv(n, nrhs, A, ia, ja, desca, ipiv, &
                          B, ib, jb, descb, info)
            integer,     intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,     intent(in)    :: desca(9), descb(9)
            integer,     intent(out)   :: ipiv(*), info
            complex(10), intent(inout) :: A(*), B(*)
        end subroutine
        subroutine pygeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            integer,     intent(in)    :: m, n, ia, ja, lwork
            integer,     intent(in)    :: desca(9)
            integer,     intent(out)   :: info
            complex(10), intent(inout) :: A(*)
            complex(10), intent(out)   :: tau(*), work(*)
        end subroutine
        subroutine pyheev(jobz, uplo, n, A, ia, ja, desca, w, &
                          Z, iz, jz, descz, work, lwork, rwork, lrwork, info)
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ia, ja, iz, jz, lwork, lrwork
            integer,     intent(in)    :: desca(9), descz(9)
            integer,     intent(out)   :: info
            complex(10), intent(inout) :: A(*)
            real(10),    intent(out)   :: w(*), rwork(*)
            complex(10), intent(out)   :: Z(*), work(*)
        end subroutine
    end interface

contains

    function nloc(desc) result(n)
        ! Local element count = lld * loccols. Test code allocates
        ! the local panel as A_loc(lld, locn), so this matches the
        ! caller's actual storage — using global N here would walk
        ! past the end of A_loc on every non-trivial grid.
        integer, intent(in) :: desc(9)
        integer :: n
        n = desc(9) * max(1, loccols(desc))
    end function

    ! Local row count = NUMROC(M, MB, myrow, RSRC, nprow).
    function locrows(desc) result(lr)
        integer, intent(in) :: desc(9)
        integer :: lr, nprow, npcol, myrow, mycol
        call blacs_gridinfo(desc(2), nprow, npcol, myrow, mycol)
        lr = numroc_local(desc(3), desc(5), myrow, desc(7), nprow)
    end function

    ! Local column count = NUMROC(N, NB, mycol, CSRC, npcol).
    function loccols(desc) result(lc)
        integer, intent(in) :: desc(9)
        integer :: lc, nprow, npcol, myrow, mycol
        call blacs_gridinfo(desc(2), nprow, npcol, myrow, mycol)
        lc = numroc_local(desc(4), desc(6), mycol, desc(8), npcol)
    end function

    subroutine target_pdgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,  intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,  intent(in)    :: desca(9), descb(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*), B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk)
        call pegesv(n, nrhs, At, ia, ja, desca, ipiv, &
                    Bt, ib, jb, descb, info)
        A(1:na) = real(At, ep); B(1:nb) = real(Bt, ep)
    end subroutine

    subroutine target_pdgetrf(m, n, A, ia, ja, desca, ipiv, info)
        integer,  intent(in)    :: m, n, ia, ja
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*)
        integer :: na
        real(tk), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = real(A(1:na), tk)
        call pegetrf(m, n, At, ia, ja, desca, ipiv, info)
        A(1:na) = real(At, ep)
    end subroutine

    subroutine target_pdgetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                              B, ib, jb, descb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9), ipiv(*)
        integer,   intent(out)   :: info
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk)
        call pegetrs(trans, n, nrhs, At, ia, ja, desca, ipiv, &
                     Bt, ib, jb, descb, info)
        B(1:nb) = real(Bt, ep)
    end subroutine

    subroutine target_pdpotrf(uplo, n, A, ia, ja, desca, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja
        integer,   intent(in)    :: desca(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        integer :: na
        real(tk), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = real(A(1:na), tk)
        call pepotrf(uplo, n, At, ia, ja, desca, info)
        A(1:na) = real(At, ep)
    end subroutine

    subroutine target_pdpotrs(uplo, n, nrhs, A, ia, ja, desca, &
                              B, ib, jb, descb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,   intent(in)    :: desca(9), descb(9)
        integer,   intent(out)   :: info
        real(ep),  intent(in)    :: A(*)
        real(ep),  intent(inout) :: B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk); Bt = real(B(1:nb), tk)
        call pepotrs(uplo, n, nrhs, At, ia, ja, desca, &
                     Bt, ib, jb, descb, info)
        B(1:nb) = real(Bt, ep)
    end subroutine

    subroutine target_pdgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, ia, ja, lwork
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: info
        real(ep), intent(inout) :: A(*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer :: na, ntau, lwt
        real(tk), allocatable :: At(:), tau_t(:), work_t(:)
        na = nloc(desca)
        ntau = max(1, loccols(desca))
        if (lwork == -1) then
            allocate(At(max(1,na)), tau_t(ntau), work_t(1))
            call pegeqrf(m, n, At, ia, ja, desca, tau_t, work_t, -1, info)
            work(1) = real(work_t(1), ep)
        else
            lwt = max(1, lwork)
            allocate(At(na), tau_t(ntau), work_t(lwt))
            At = real(A(1:na), tk)
            call pegeqrf(m, n, At, ia, ja, desca, tau_t, work_t, lwt, info)
            A(1:na) = real(At, ep)
            tau(1:ntau) = real(tau_t, ep)
        end if
    end subroutine

    subroutine target_pdsyev(jobz, uplo, n, A, ia, ja, desca, w, &
                             Z, iz, jz, descz, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ia, ja, iz, jz, lwork
        integer,   intent(in)    :: desca(9), descz(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        real(ep),  intent(out)   :: w(*), Z(*), work(*)
        integer :: na, nz, lwt
        real(tk), allocatable :: At(:), Zt(:), wt(:), work_t(:)
        na = nloc(desca); nz = nloc(descz)
        if (lwork == -1) then
            allocate(At(1), Zt(1), wt(1), work_t(1))
            call pesyev(jobz, uplo, n, At, ia, ja, desca, wt, &
                        Zt, iz, jz, descz, work_t, -1, info)
            work(1) = real(work_t(1), ep)
        else
            lwt = max(1, lwork)
            allocate(At(na), Zt(nz), wt(n), work_t(lwt))
            At = real(A(1:na), tk)
            call pesyev(jobz, uplo, n, At, ia, ja, desca, wt, &
                        Zt, iz, jz, descz, work_t, lwt, info)
            w(1:n) = real(wt, ep)
            Z(1:nz) = real(Zt, ep)
        end if
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
        integer :: na, nu, nvt, ns, lwt
        real(tk), allocatable :: At(:), Ut(:), VTt(:), st(:), work_t(:)
        na = nloc(desca); nu = nloc(descu); nvt = nloc(descvt)
        ns = min(m, n)
        if (lwork == -1) then
            allocate(At(1), Ut(1), VTt(1), st(1), work_t(1))
            call pegesvd(jobu, jobvt, m, n, At, ia, ja, desca, st, &
                         Ut, iu, ju, descu, VTt, ivt, jvt, descvt, &
                         work_t, -1, info)
            work(1) = real(work_t(1), ep)
        else
            lwt = max(1, lwork)
            allocate(At(na), Ut(nu), VTt(nvt), st(ns), work_t(lwt))
            At = real(A(1:na), tk)
            call pegesvd(jobu, jobvt, m, n, At, ia, ja, desca, st, &
                         Ut, iu, ju, descu, VTt, ivt, jvt, descvt, &
                         work_t, lwt, info)
            s(1:ns) = real(st, ep)
            U(1:nu) = real(Ut, ep)
            VT(1:nvt) = real(VTt, ep)
        end if
    end subroutine

    function target_pdlange(norm, m, n, A, ia, ja, desca, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, ia, ja
        integer,   intent(in) :: desca(9)
        real(ep),  intent(in) :: A(*)
        real(ep)              :: work(*)
        real(ep) :: r
        integer :: na, nw
        real(tk), allocatable :: At(:), work_t(:)
        real(tk) :: r_t
        na = nloc(desca)
        nw = max(m, n, 1)
        allocate(At(na), work_t(nw))
        At = real(A(1:na), tk)
        r_t = pelange(norm, m, n, At, ia, ja, desca, work_t)
        r = real(r_t, ep)
    end function

    subroutine target_pdlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: m, n, ia, ja, ib, jb
        integer,   intent(in)  :: desca(9), descb(9)
        real(ep),  intent(in)  :: A(*)
        real(ep),  intent(out) :: B(*)
        integer :: na, nb
        real(tk), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = real(A(1:na), tk)
        call pelacpy(uplo, m, n, At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = real(Bt, ep)
    end subroutine

    subroutine target_pdlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, ia, ja
        integer,   intent(in)    :: desca(9)
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(*)
        integer :: na
        real(tk), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = real(A(1:na), tk)
        call pelaset(uplo, m, n, real(alpha, tk), real(beta, tk), &
                     At, ia, ja, desca)
        A(1:na) = real(At, ep)
    end subroutine

    subroutine target_pzgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,     intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,     intent(in)    :: desca(9), descb(9)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(inout) :: A(*), B(*)
        integer :: na, nb
        complex(tk), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = cmplx(A(1:na), kind=tk); Bt = cmplx(B(1:nb), kind=tk)
        call pygesv(n, nrhs, At, ia, ja, desca, ipiv, &
                    Bt, ib, jb, descb, info)
        A(1:na) = cmplx(At, kind=ep); B(1:nb) = cmplx(Bt, kind=ep)
    end subroutine

    subroutine target_pzgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, ia, ja, lwork
        integer,     intent(in)    :: desca(9)
        integer,     intent(out)   :: info
        complex(ep), intent(inout) :: A(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer :: na, ntau, lwt
        complex(tk), allocatable :: At(:), tau_t(:), work_t(:)
        na = nloc(desca); ntau = max(1, loccols(desca))
        if (lwork == -1) then
            allocate(At(1), tau_t(1), work_t(1))
            call pygeqrf(m, n, At, ia, ja, desca, tau_t, work_t, -1, info)
            work(1) = cmplx(work_t(1), kind=ep)
        else
            lwt = max(1, lwork)
            allocate(At(na), tau_t(ntau), work_t(lwt))
            At = cmplx(A(1:na), kind=tk)
            call pygeqrf(m, n, At, ia, ja, desca, tau_t, work_t, lwt, info)
            A(1:na) = cmplx(At, kind=ep)
            tau(1:ntau) = cmplx(tau_t, kind=ep)
        end if
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
        integer :: na, nz, lwt, lrwt
        complex(tk), allocatable :: At(:), Zt(:), work_t(:)
        real(tk),    allocatable :: wt(:), rwork_t(:)
        na = nloc(desca); nz = nloc(descz)
        if (lwork == -1) then
            allocate(At(1), Zt(1), wt(1), work_t(1), rwork_t(1))
            call pyheev(jobz, uplo, n, At, ia, ja, desca, wt, &
                        Zt, iz, jz, descz, work_t, -1, rwork_t, -1, info)
            work(1) = cmplx(work_t(1), kind=ep)
            rwork(1) = real(rwork_t(1), ep)
        else
            lwt = max(1, lwork); lrwt = max(1, lrwork)
            allocate(At(na), Zt(nz), wt(n), work_t(lwt), rwork_t(lrwt))
            At = cmplx(A(1:na), kind=tk)
            call pyheev(jobz, uplo, n, At, ia, ja, desca, wt, &
                        Zt, iz, jz, descz, work_t, lwt, rwork_t, lrwt, info)
            w(1:n) = real(wt, ep)
            Z(1:nz) = cmplx(Zt, kind=ep)
        end if
    end subroutine

end module target_scalapack
