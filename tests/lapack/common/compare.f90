module compare
    use prec_kinds, only: ep
    implicit none
    private
    public :: rel_err_scalar, rel_err_scalar_z
    public :: max_rel_err_vec, max_rel_err_vec_z
    public :: max_rel_err_mat, max_rel_err_mat_z
contains

    pure function rel_err_scalar(a, b) result(err)
        real(ep), intent(in) :: a, b
        real(ep) :: err, denom
        denom = abs(b)
        if (denom < tiny(1.0_ep)) then
            err = abs(a - b)
        else
            err = abs(a - b) / denom
        end if
    end function rel_err_scalar

    pure function rel_err_scalar_z(a, b) result(err)
        complex(ep), intent(in) :: a, b
        real(ep) :: err, denom
        denom = abs(b)
        if (denom < tiny(1.0_ep)) then
            err = abs(a - b)
        else
            err = abs(a - b) / denom
        end if
    end function rel_err_scalar_z

    pure function max_rel_err_vec(a, b) result(err)
        real(ep), intent(in) :: a(:), b(:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_vec

    pure function max_rel_err_vec_z(a, b) result(err)
        complex(ep), intent(in) :: a(:), b(:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_vec_z

    pure function max_rel_err_mat(a, b) result(err)
        real(ep), intent(in) :: a(:,:), b(:,:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_mat

    pure function max_rel_err_mat_z(a, b) result(err)
        complex(ep), intent(in) :: a(:,:), b(:,:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_mat_z

end module compare
