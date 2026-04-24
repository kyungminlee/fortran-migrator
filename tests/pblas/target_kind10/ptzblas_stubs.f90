! Link-time stubs for ZZDOTC / ZZDOTU in the kind10 (epblas) build.
! See target_kind16/ptzblas_stubs.f90 for the full rationale; briefly:
! the migrator leaves CALL ZZDOTC unchanged in yvvdotc.f / yvvdotu.f
! because TOOLS/zzdotc.f isn't in the PTZBLAS scan set, so we forward
! to ydotc / ydotu (kind10 complex dot products).

subroutine zzdotc(n, dotc, x, incx, y, incy)
    implicit none
    integer,     intent(in)  :: n, incx, incy
    complex(10), intent(out) :: dotc
    complex(10), intent(in)  :: x(*), y(*)
    complex(10), external    :: ydotc
    dotc = ydotc(n, x, incx, y, incy)
end subroutine

subroutine zzdotu(n, dotu, x, incx, y, incy)
    implicit none
    integer,     intent(in)  :: n, incx, incy
    complex(10), intent(out) :: dotu
    complex(10), intent(in)  :: x(*), y(*)
    complex(10), external    :: ydotu
    dotu = ydotu(n, x, incx, y, incy)
end subroutine
