module test_data
    use prec_kinds, only: ep
    implicit none
    private
    public :: gen_vector_quad, gen_matrix_quad
    public :: gen_vector_complex, gen_matrix_complex
    public :: gen_spd_matrix_quad, gen_hermitian_matrix_quad
    public :: gen_symmetric_matrix_quad

contains

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

    ! Symmetric positive-definite matrix: A = X*X^T + n*I. The diagonal
    ! shift guarantees SPD even for rank-deficient X. Output is n x n.
    subroutine gen_spd_matrix_quad(n, A, seed)
        integer,  intent(in)  :: n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(ep), allocatable :: X(:,:)
        integer :: i

        call gen_matrix_quad(n, n, X, seed)
        allocate(A(n, n))
        A = matmul(X, transpose(X))
        do i = 1, n
            A(i, i) = A(i, i) + real(n, ep)
        end do
    end subroutine gen_spd_matrix_quad

    ! Real symmetric matrix: A = (X + X^T)/2.
    subroutine gen_symmetric_matrix_quad(n, A, seed)
        integer,  intent(in)  :: n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(ep), allocatable :: X(:,:)

        call gen_matrix_quad(n, n, X, seed)
        allocate(A(n, n))
        A = 0.5_ep * (X + transpose(X))
    end subroutine gen_symmetric_matrix_quad

    ! Hermitian matrix: A = (X + X^H)/2.
    subroutine gen_hermitian_matrix_quad(n, A, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        complex(ep), allocatable :: X(:,:)

        call gen_matrix_complex(n, n, X, seed)
        allocate(A(n, n))
        A = 0.5_ep * (X + transpose(conjg(X)))
    end subroutine gen_hermitian_matrix_quad

end module test_data
