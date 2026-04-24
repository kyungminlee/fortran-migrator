! Link-time stubs for ZZDOTC / ZZDOTU referenced by the migrated
! xvvdotc.f / xvvdotu.f in qptzblas.
!
! The upstream zvvdotc.f / zvvdotu.f call ZZDOTC / ZZDOTU — a pair of
! ScaLAPACK-TOOLS Fortran wrappers that return the dot product through
! an output argument rather than a function result. The migrator's
! PTZBLAS recipe scans only PTZBLAS itself, so the reference to
! ZZDOTC inside xvvdotc doesn't get renamed to the kind16-appropriate
! symbol, and TOOLS/zzdotc.f itself isn't part of any migrated library.
! These stubs resolve the dangling reference by forwarding to the
! kind16 canonical complex dot products (xdotc / xdotu).
!
! The same gap exists for every target — see the sibling
! target_kind10 and target_multifloats wrappers for those variants.

subroutine zzdotc(n, dotc, x, incx, y, incy)
    implicit none
    integer,     intent(in)  :: n, incx, incy
    complex(16), intent(out) :: dotc
    complex(16), intent(in)  :: x(*), y(*)
    complex(16), external    :: xdotc
    dotc = xdotc(n, x, incx, y, incy)
end subroutine

subroutine zzdotu(n, dotu, x, incx, y, incy)
    implicit none
    integer,     intent(in)  :: n, incx, incy
    complex(16), intent(out) :: dotu
    complex(16), intent(in)  :: x(*), y(*)
    complex(16), external    :: xdotu
    dotu = xdotu(n, x, incx, y, incy)
end subroutine
