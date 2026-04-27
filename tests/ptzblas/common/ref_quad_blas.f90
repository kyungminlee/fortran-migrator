! Explicit interfaces to refblas_quad — vendored Netlib BLAS compiled
! at REAL(KIND=16). PTZBLAS tests use the BLAS analogues for vasum /
! vvdot / dzvasum which match dasum / ddot / dzasum byte-for-byte.
!
! Module name is ptzblas_ref_quad_blas (not ref_quad_blas) to avoid
! .mod collisions with tests/blas and tests/pblas equivalents — all
! targets write into the same ${PROJECT_BINARY_DIR}/fmod directory.
module ptzblas_ref_quad_blas
    use prec_kinds, only: ep
    implicit none

    interface
        function ddot(n, x, incx, y, incy) result(r)
            import :: ep
            integer,  intent(in) :: n, incx, incy
            real(ep), intent(in) :: x(*), y(*)
            real(ep) :: r
        end function ddot

        function dasum(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dasum

        function zdotc(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotc

        function zdotu(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotu

        function dzasum(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dzasum
    end interface
end module ptzblas_ref_quad_blas
