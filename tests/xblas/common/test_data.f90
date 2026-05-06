module test_data
    use prec_kinds, only: ep
    implicit none
    private
    public :: gen_vector_quad, gen_matrix_quad
    public :: gen_vector_complex, gen_matrix_complex

contains

    ! Seed Fortran's intrinsic RNG with a value derived from `seed` so
    ! every test gets reproducible inputs that differ across cases.
    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: s(:)
        call random_seed(size=n)
        allocate(s(n))
        do i = 1, n
            s(i) = seed + i * 1000003
        end do
        call random_seed(put=s)
    end subroutine seed_rng

    subroutine gen_vector_quad(n, x, seed)
        integer,   intent(in)  :: n, seed
        real(ep),  intent(out), allocatable :: x(:)
        real(8),   allocatable :: r(:)
        integer :: i

        call seed_rng(seed)
        allocate(r(n), x(n))
        call random_number(r)
        do i = 1, n
            x(i) = real(r(i), ep) * 2.0_ep - 1.0_ep
        end do
    end subroutine gen_vector_quad

    subroutine gen_matrix_quad(m, n, A, seed)
        integer,  intent(in)  :: m, n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(8),  allocatable :: R(:,:)
        integer :: i, j

        call seed_rng(seed)
        allocate(R(m, n), A(m, n))
        call random_number(R)
        do j = 1, n
            do i = 1, m
                A(i, j) = real(R(i, j), ep) * 2.0_ep - 1.0_ep
            end do
        end do
    end subroutine gen_matrix_quad

    subroutine gen_vector_complex(n, x, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: x(:)
        real(ep),    allocatable :: re(:), im(:)

        call gen_vector_quad(n, re, seed)
        call gen_vector_quad(n, im, seed + 1)
        allocate(x(n))
        x = cmplx(re, im, ep)
    end subroutine gen_vector_complex

    subroutine gen_matrix_complex(m, n, A, seed)
        integer,     intent(in)  :: m, n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        real(ep),    allocatable :: R(:,:), I_(:,:)

        call gen_matrix_quad(m, n, R,  seed)
        call gen_matrix_quad(m, n, I_, seed + 1)
        allocate(A(m, n))
        A = cmplx(R, I_, ep)
    end subroutine gen_matrix_complex

end module test_data
