! Per-target ScaLAPACK wrapper for the multifloats (ddscalapack) build.
!
! Test code works in REAL(KIND=ep)=KIND=16. The wrappers split each
! quad value into a TYPE(real64x2) double-double, call the migrated
! dd/zz-prefix ScaLAPACK routines, then recombine back to quad. Same
! conversion pattern as tests/pblas/target_multifloats and the same
! NUMROC-based local sizing as tests/scalapack/target_kind10.

module target_scalapack
    use prec_kinds,  only: ep, dp
    use multifloats, only: real64x2, cmplx64x2
    use pblas_grid,  only: numroc_local
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

    character(len=*), parameter :: target_name = 'multifloats'
    real(ep),         parameter :: target_eps  = real(2.0_dp**(-104), ep)

    interface
        subroutine blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
            integer, intent(in)  :: icontxt
            integer, intent(out) :: nprow, npcol, myrow, mycol
        end subroutine
    end interface

    interface
        subroutine pddgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                           B, ib, jb, descb, info)
            import :: real64x2
            integer, intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer, intent(in)    :: desca(9), descb(9)
            integer, intent(out)   :: ipiv(*), info
            type(real64x2), intent(inout) :: A(*), B(*)
        end subroutine
        subroutine pddgetrf(m, n, A, ia, ja, desca, ipiv, info)
            import :: real64x2
            integer, intent(in)    :: m, n, ia, ja
            integer, intent(in)    :: desca(9)
            integer, intent(out)   :: ipiv(*), info
            type(real64x2), intent(inout) :: A(*)
        end subroutine
        subroutine pddgetrs(trans, n, nrhs, A, ia, ja, desca, ipiv, &
                            B, ib, jb, descb, info)
            import :: real64x2
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9), ipiv(*)
            integer,   intent(out)   :: info
            type(real64x2), intent(in)    :: A(*)
            type(real64x2), intent(inout) :: B(*)
        end subroutine
        subroutine pddpotrf(uplo, n, A, ia, ja, desca, info)
            import :: real64x2
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ia, ja
            integer,   intent(in)    :: desca(9)
            integer,   intent(out)   :: info
            type(real64x2), intent(inout) :: A(*)
        end subroutine
        subroutine pddpotrs(uplo, n, nrhs, A, ia, ja, desca, &
                            B, ib, jb, descb, info)
            import :: real64x2
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer,   intent(in)    :: desca(9), descb(9)
            integer,   intent(out)   :: info
            type(real64x2), intent(in)    :: A(*)
            type(real64x2), intent(inout) :: B(*)
        end subroutine
        subroutine pddgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            import :: real64x2
            integer, intent(in)    :: m, n, ia, ja, lwork
            integer, intent(in)    :: desca(9)
            integer, intent(out)   :: info
            type(real64x2), intent(inout) :: A(*)
            type(real64x2), intent(out)   :: tau(*), work(*)
        end subroutine
        subroutine pddsyev(jobz, uplo, n, A, ia, ja, desca, w, &
                           Z, iz, jz, descz, work, lwork, info)
            import :: real64x2
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ia, ja, iz, jz, lwork
            integer,   intent(in)    :: desca(9), descz(9)
            integer,   intent(out)   :: info
            type(real64x2), intent(inout) :: A(*)
            type(real64x2), intent(out)   :: w(*), Z(*), work(*)
        end subroutine
        subroutine pddgesvd(jobu, jobvt, m, n, A, ia, ja, desca, s, &
                            U, iu, ju, descu, VT, ivt, jvt, descvt, &
                            work, lwork, info)
            import :: real64x2
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, ia, ja, iu, ju, ivt, jvt, lwork
            integer,   intent(in)    :: desca(9), descu(9), descvt(9)
            integer,   intent(out)   :: info
            type(real64x2), intent(inout) :: A(*)
            type(real64x2), intent(out)   :: s(*), U(*), VT(*), work(*)
        end subroutine
        function pddlange(norm, m, n, A, ia, ja, desca, work) result(r)
            import :: real64x2
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, ia, ja
            integer,   intent(in) :: desca(9)
            type(real64x2), intent(in) :: A(*)
            type(real64x2)             :: work(*)
            type(real64x2) :: r
        end function
        subroutine pddlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
            import :: real64x2
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, ia, ja, ib, jb
            integer,   intent(in)  :: desca(9), descb(9)
            type(real64x2), intent(in)  :: A(*)
            type(real64x2), intent(out) :: B(*)
        end subroutine
        subroutine pddlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
            import :: real64x2
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, ia, ja
            integer,   intent(in)    :: desca(9)
            type(real64x2), intent(in)    :: alpha, beta
            type(real64x2), intent(inout) :: A(*)
        end subroutine

        subroutine pzzgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                           B, ib, jb, descb, info)
            import :: cmplx64x2
            integer, intent(in)    :: n, nrhs, ia, ja, ib, jb
            integer, intent(in)    :: desca(9), descb(9)
            integer, intent(out)   :: ipiv(*), info
            type(cmplx64x2), intent(inout) :: A(*), B(*)
        end subroutine
        subroutine pzzgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
            import :: cmplx64x2
            integer, intent(in)    :: m, n, ia, ja, lwork
            integer, intent(in)    :: desca(9)
            integer, intent(out)   :: info
            type(cmplx64x2), intent(inout) :: A(*)
            type(cmplx64x2), intent(out)   :: tau(*), work(*)
        end subroutine
        subroutine pzzheev(jobz, uplo, n, A, ia, ja, desca, w, &
                           Z, iz, jz, descz, work, lwork, rwork, lrwork, info)
            import :: real64x2, cmplx64x2
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ia, ja, iz, jz, lwork, lrwork
            integer,   intent(in)    :: desca(9), descz(9)
            integer,   intent(out)   :: info
            type(cmplx64x2), intent(inout) :: A(*)
            type(real64x2),   intent(out)   :: w(*), rwork(*)
            type(cmplx64x2), intent(out)   :: Z(*), work(*)
        end subroutine
    end interface

contains

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

    function locrows(desc) result(lr)
        integer, intent(in) :: desc(9)
        integer :: lr, nprow, npcol, myrow, mycol
        call blacs_gridinfo(desc(2), nprow, npcol, myrow, mycol)
        lr = numroc_local(desc(3), desc(5), myrow, desc(7), nprow)
    end function

    function loccols(desc) result(lc)
        integer, intent(in) :: desc(9)
        integer :: lc, nprow, npcol, myrow, mycol
        call blacs_gridinfo(desc(2), nprow, npcol, myrow, mycol)
        lc = numroc_local(desc(4), desc(6), mycol, desc(8), npcol)
    end function

    function nloc(desc) result(n)
        integer, intent(in) :: desc(9)
        integer :: n
        n = desc(9) * max(1, loccols(desc))
    end function

    subroutine target_pdgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,  intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,  intent(in)    :: desca(9), descb(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*), B(*)
        integer :: na, nb
        type(real64x2), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb))
        call pddgesv(n, nrhs, At, ia, ja, desca, ipiv, &
                     Bt, ib, jb, descb, info)
        A(1:na) = dd2q(At); B(1:nb) = dd2q(Bt)
    end subroutine

    subroutine target_pdgetrf(m, n, A, ia, ja, desca, ipiv, info)
        integer,  intent(in)    :: m, n, ia, ja
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: ipiv(*), info
        real(ep), intent(inout) :: A(*)
        integer :: na
        type(real64x2), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = q2dd(A(1:na))
        call pddgetrf(m, n, At, ia, ja, desca, ipiv, info)
        A(1:na) = dd2q(At)
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
        type(real64x2), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb))
        call pddgetrs(trans, n, nrhs, At, ia, ja, desca, ipiv, &
                      Bt, ib, jb, descb, info)
        B(1:nb) = dd2q(Bt)
    end subroutine

    subroutine target_pdpotrf(uplo, n, A, ia, ja, desca, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ia, ja
        integer,   intent(in)    :: desca(9)
        integer,   intent(out)   :: info
        real(ep),  intent(inout) :: A(*)
        integer :: na
        type(real64x2), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = q2dd(A(1:na))
        call pddpotrf(uplo, n, At, ia, ja, desca, info)
        A(1:na) = dd2q(At)
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
        type(real64x2), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na)); Bt = q2dd(B(1:nb))
        call pddpotrs(uplo, n, nrhs, At, ia, ja, desca, &
                      Bt, ib, jb, descb, info)
        B(1:nb) = dd2q(Bt)
    end subroutine

    subroutine target_pdgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, ia, ja, lwork
        integer,  intent(in)    :: desca(9)
        integer,  intent(out)   :: info
        real(ep), intent(inout) :: A(*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer :: na, ntau, lwt
        type(real64x2), allocatable :: At(:), tau_t(:), work_t(:)
        na = nloc(desca)
        ntau = max(1, loccols(desca))
        if (lwork == -1) then
            allocate(At(max(1,na)), tau_t(ntau), work_t(1))
            call pddgeqrf(m, n, At, ia, ja, desca, tau_t, work_t, -1, info)
            work(1) = dd2q(work_t(1))
        else
            lwt = max(1, lwork)
            allocate(At(na), tau_t(ntau), work_t(lwt))
            At = q2dd(A(1:na))
            call pddgeqrf(m, n, At, ia, ja, desca, tau_t, work_t, lwt, info)
            A(1:na) = dd2q(At)
            tau(1:ntau) = dd2q(tau_t)
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
        type(real64x2), allocatable :: At(:), Zt(:), wt(:), work_t(:)
        na = nloc(desca); nz = nloc(descz)
        if (lwork == -1) then
            allocate(At(1), Zt(1), wt(1), work_t(1))
            call pddsyev(jobz, uplo, n, At, ia, ja, desca, wt, &
                         Zt, iz, jz, descz, work_t, -1, info)
            work(1) = dd2q(work_t(1))
        else
            lwt = max(1, lwork)
            allocate(At(na), Zt(nz), wt(n), work_t(lwt))
            At = q2dd(A(1:na))
            call pddsyev(jobz, uplo, n, At, ia, ja, desca, wt, &
                         Zt, iz, jz, descz, work_t, lwt, info)
            w(1:n) = dd2q(wt)
            Z(1:nz) = dd2q(Zt)
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
        type(real64x2), allocatable :: At(:), Ut(:), VTt(:), st(:), work_t(:)
        na = nloc(desca); nu = nloc(descu); nvt = nloc(descvt)
        ns = min(m, n)
        if (lwork == -1) then
            allocate(At(1), Ut(1), VTt(1), st(1), work_t(1))
            call pddgesvd(jobu, jobvt, m, n, At, ia, ja, desca, st, &
                          Ut, iu, ju, descu, VTt, ivt, jvt, descvt, &
                          work_t, -1, info)
            work(1) = dd2q(work_t(1))
        else
            lwt = max(1, lwork)
            allocate(At(na), Ut(nu), VTt(nvt), st(ns), work_t(lwt))
            At = q2dd(A(1:na))
            call pddgesvd(jobu, jobvt, m, n, At, ia, ja, desca, st, &
                          Ut, iu, ju, descu, VTt, ivt, jvt, descvt, &
                          work_t, lwt, info)
            s(1:ns) = dd2q(st)
            U(1:nu) = dd2q(Ut)
            VT(1:nvt) = dd2q(VTt)
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
        type(real64x2), allocatable :: At(:), work_t(:)
        type(real64x2) :: r_t
        na = nloc(desca)
        nw = max(m, n, 1)
        allocate(At(na), work_t(nw))
        At = q2dd(A(1:na))
        r_t = pddlange(norm, m, n, At, ia, ja, desca, work_t)
        r = dd2q(r_t)
    end function

    subroutine target_pdlacpy(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: m, n, ia, ja, ib, jb
        integer,   intent(in)  :: desca(9), descb(9)
        real(ep),  intent(in)  :: A(*)
        real(ep),  intent(out) :: B(*)
        integer :: na, nb
        type(real64x2), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = q2dd(A(1:na))
        call pddlacpy(uplo, m, n, At, ia, ja, desca, Bt, ib, jb, descb)
        B(1:nb) = dd2q(Bt)
    end subroutine

    subroutine target_pdlaset(uplo, m, n, alpha, beta, A, ia, ja, desca)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, ia, ja
        integer,   intent(in)    :: desca(9)
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(*)
        integer :: na
        type(real64x2), allocatable :: At(:)
        na = nloc(desca)
        allocate(At(na))
        At = q2dd(A(1:na))
        call pddlaset(uplo, m, n, q2dd(alpha), q2dd(beta), &
                      At, ia, ja, desca)
        A(1:na) = dd2q(At)
    end subroutine

    subroutine target_pzgesv(n, nrhs, A, ia, ja, desca, ipiv, &
                             B, ib, jb, descb, info)
        integer,     intent(in)    :: n, nrhs, ia, ja, ib, jb
        integer,     intent(in)    :: desca(9), descb(9)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(inout) :: A(*), B(*)
        integer :: na, nb
        type(cmplx64x2), allocatable :: At(:), Bt(:)
        na = nloc(desca); nb = nloc(descb)
        allocate(At(na), Bt(nb))
        At = q2zz(A(1:na)); Bt = q2zz(B(1:nb))
        call pzzgesv(n, nrhs, At, ia, ja, desca, ipiv, &
                     Bt, ib, jb, descb, info)
        A(1:na) = zz2q(At); B(1:nb) = zz2q(Bt)
    end subroutine

    subroutine target_pzgeqrf(m, n, A, ia, ja, desca, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, ia, ja, lwork
        integer,     intent(in)    :: desca(9)
        integer,     intent(out)   :: info
        complex(ep), intent(inout) :: A(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer :: na, ntau, lwt
        type(cmplx64x2), allocatable :: At(:), tau_t(:), work_t(:)
        na = nloc(desca); ntau = max(1, loccols(desca))
        if (lwork == -1) then
            allocate(At(max(1,na)), tau_t(ntau), work_t(1))
            call pzzgeqrf(m, n, At, ia, ja, desca, tau_t, work_t, -1, info)
            work(1) = zz2q(work_t(1))
        else
            lwt = max(1, lwork)
            allocate(At(na), tau_t(ntau), work_t(lwt))
            At = q2zz(A(1:na))
            call pzzgeqrf(m, n, At, ia, ja, desca, tau_t, work_t, lwt, info)
            A(1:na) = zz2q(At)
            tau(1:ntau) = zz2q(tau_t)
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
        type(cmplx64x2), allocatable :: At(:), Zt(:), work_t(:)
        type(real64x2),   allocatable :: wt(:), rwork_t(:)
        na = nloc(desca); nz = nloc(descz)
        if (lwork == -1) then
            allocate(At(1), Zt(1), wt(1), work_t(1), rwork_t(1))
            call pzzheev(jobz, uplo, n, At, ia, ja, desca, wt, &
                         Zt, iz, jz, descz, work_t, -1, rwork_t, -1, info)
            work(1) = zz2q(work_t(1))
            rwork(1) = dd2q(rwork_t(1))
        else
            lwt = max(1, lwork); lrwt = max(1, lrwork)
            allocate(At(na), Zt(nz), wt(n), work_t(lwt), rwork_t(lrwt))
            At = q2zz(A(1:na))
            call pzzheev(jobz, uplo, n, At, ia, ja, desca, wt, &
                         Zt, iz, jz, descz, work_t, lwt, rwork_t, lrwt, info)
            w(1:n) = dd2q(wt)
            Z(1:nz) = zz2q(Zt)
        end if
    end subroutine

end module target_scalapack
