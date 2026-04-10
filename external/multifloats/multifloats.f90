











module multifloats
    implicit none
    private

    integer, parameter :: qp = 16  ! Quad precision
    integer, parameter :: dp = 8   ! Double precision
    integer, parameter :: sp = 4   ! Single precision

    ! SEQUENCE attribute is required so derived-type variables of these
    ! types may appear in EQUIVALENCE statements (used by some LAPACK
    ! routines such as DLALN2/wlaln2 to overlay 2x2 work arrays with
    ! linearized vectors). The DTIO procedures below use ``type(...)``
    ! (allowed by F2018 for non-extensible types) instead of the
    ! customary ``class(...)``, since SEQUENCE types may not be
    ! polymorphic.
    type, public :: float64x2
        sequence
        real(dp) :: limbs(2)
    end type float64x2

    type, public :: complex128x2
        sequence
        type(float64x2) :: re, im
    end type complex128x2

    ! ================================================================
    ! Constructors
    ! ================================================================

    interface float64x2
        module procedure mf_from_dp
        module procedure mf_from_sp
        module procedure mf_from_int
        module procedure mf_from_char
        module procedure mf_from_cx
        module procedure mf_from_cdp
        module procedure mf_from_csp
    end interface float64x2

    interface complex128x2
        module procedure cx_from_mf_2
        module procedure cx_from_dp_2
        module procedure cx_from_mf_1
        module procedure cx_from_dp_1
    end interface complex128x2

    ! ================================================================
    ! Arithmetic operators (+, -, *, /)
    ! ================================================================

    interface operator(+)
        module procedure mf_add_mf
        module procedure mf_add_dp
        module procedure mf_add_sp
        module procedure mf_add_int
        module procedure mf_add_cx
        module procedure mf_add_cdp
        module procedure mf_add_csp
        module procedure dp_add_mf
        module procedure dp_add_cx
        module procedure sp_add_mf
        module procedure sp_add_cx
        module procedure int_add_mf
        module procedure int_add_cx
        module procedure cx_add_mf
        module procedure cx_add_dp
        module procedure cx_add_sp
        module procedure cx_add_int
        module procedure cx_add_cx
        module procedure cx_add_cdp
        module procedure cx_add_csp
        module procedure cdp_add_mf
        module procedure cdp_add_cx
        module procedure csp_add_mf
        module procedure csp_add_cx
    end interface
    public :: operator(+)

    interface operator(-)
        module procedure mf_sub_mf
        module procedure mf_sub_dp
        module procedure mf_sub_sp
        module procedure mf_sub_int
        module procedure mf_sub_cx
        module procedure mf_sub_cdp
        module procedure mf_sub_csp
        module procedure dp_sub_mf
        module procedure dp_sub_cx
        module procedure sp_sub_mf
        module procedure sp_sub_cx
        module procedure int_sub_mf
        module procedure int_sub_cx
        module procedure cx_sub_mf
        module procedure cx_sub_dp
        module procedure cx_sub_sp
        module procedure cx_sub_int
        module procedure cx_sub_cx
        module procedure cx_sub_cdp
        module procedure cx_sub_csp
        module procedure cdp_sub_mf
        module procedure cdp_sub_cx
        module procedure csp_sub_mf
        module procedure csp_sub_cx
        module procedure mf_neg
        module procedure cx_neg
    end interface
    public :: operator(-)

    interface operator(*)
        module procedure mf_mul_mf
        module procedure mf_mul_dp
        module procedure mf_mul_sp
        module procedure mf_mul_int
        module procedure mf_mul_cx
        module procedure mf_mul_cdp
        module procedure mf_mul_csp
        module procedure dp_mul_mf
        module procedure dp_mul_cx
        module procedure sp_mul_mf
        module procedure sp_mul_cx
        module procedure int_mul_mf
        module procedure int_mul_cx
        module procedure cx_mul_mf
        module procedure cx_mul_dp
        module procedure cx_mul_sp
        module procedure cx_mul_int
        module procedure cx_mul_cx
        module procedure cx_mul_cdp
        module procedure cx_mul_csp
        module procedure cdp_mul_mf
        module procedure cdp_mul_cx
        module procedure csp_mul_mf
        module procedure csp_mul_cx
    end interface
    public :: operator(*)

    interface operator(/)
        module procedure mf_div_mf
        module procedure mf_div_dp
        module procedure mf_div_sp
        module procedure mf_div_int
        module procedure mf_div_cx
        module procedure mf_div_cdp
        module procedure mf_div_csp
        module procedure dp_div_mf
        module procedure dp_div_cx
        module procedure sp_div_mf
        module procedure sp_div_cx
        module procedure int_div_mf
        module procedure int_div_cx
        module procedure cx_div_mf
        module procedure cx_div_dp
        module procedure cx_div_sp
        module procedure cx_div_int
        module procedure cx_div_cx
        module procedure cx_div_cdp
        module procedure cx_div_csp
        module procedure cdp_div_mf
        module procedure cdp_div_cx
        module procedure csp_div_mf
        module procedure csp_div_cx
    end interface
    public :: operator(/)

    ! ================================================================
    ! Power operator (**)
    ! ================================================================

    interface operator(**)
        module procedure mf_pow_mf
        module procedure mf_pow_dp
        module procedure mf_pow_sp
        module procedure mf_pow_int
        module procedure mf_pow_cx
        module procedure mf_pow_cdp
        module procedure mf_pow_csp
        module procedure dp_pow_mf
        module procedure dp_pow_cx
        module procedure sp_pow_mf
        module procedure sp_pow_cx
        module procedure int_pow_mf
        module procedure int_pow_cx
        module procedure cx_pow_mf
        module procedure cx_pow_dp
        module procedure cx_pow_sp
        module procedure cx_pow_int
        module procedure cx_pow_cx
        module procedure cx_pow_cdp
        module procedure cx_pow_csp
        module procedure cdp_pow_mf
        module procedure cdp_pow_cx
        module procedure csp_pow_mf
        module procedure csp_pow_cx
    end interface
    public :: operator(**)

    ! ================================================================
    ! Equality comparison (==, /=) — all type pairs
    ! ================================================================

    interface operator(==)
        module procedure mf_eq_mf
        module procedure mf_eq_dp
        module procedure mf_eq_sp
        module procedure mf_eq_int
        module procedure mf_eq_cx
        module procedure mf_eq_cdp
        module procedure mf_eq_csp
        module procedure dp_eq_mf
        module procedure dp_eq_cx
        module procedure sp_eq_mf
        module procedure sp_eq_cx
        module procedure int_eq_mf
        module procedure int_eq_cx
        module procedure cx_eq_mf
        module procedure cx_eq_dp
        module procedure cx_eq_sp
        module procedure cx_eq_int
        module procedure cx_eq_cx
        module procedure cx_eq_cdp
        module procedure cx_eq_csp
        module procedure cdp_eq_mf
        module procedure cdp_eq_cx
        module procedure csp_eq_mf
        module procedure csp_eq_cx
    end interface
    public :: operator(==)

    interface operator(/=)
        module procedure mf_ne_mf
        module procedure mf_ne_dp
        module procedure mf_ne_sp
        module procedure mf_ne_int
        module procedure mf_ne_cx
        module procedure mf_ne_cdp
        module procedure mf_ne_csp
        module procedure dp_ne_mf
        module procedure dp_ne_cx
        module procedure sp_ne_mf
        module procedure sp_ne_cx
        module procedure int_ne_mf
        module procedure int_ne_cx
        module procedure cx_ne_mf
        module procedure cx_ne_dp
        module procedure cx_ne_sp
        module procedure cx_ne_int
        module procedure cx_ne_cx
        module procedure cx_ne_cdp
        module procedure cx_ne_csp
        module procedure cdp_ne_mf
        module procedure cdp_ne_cx
        module procedure csp_ne_mf
        module procedure csp_ne_cx
    end interface
    public :: operator(/=)

    ! ================================================================
    ! Ordered comparison (<, >, <=, >=) — real types only
    ! ================================================================

    interface operator(<)
        module procedure mf_lt_mf
        module procedure mf_lt_dp
        module procedure mf_lt_sp
        module procedure mf_lt_int
        module procedure dp_lt_mf
        module procedure sp_lt_mf
        module procedure int_lt_mf
    end interface
    public :: operator(<)

    interface operator(>)
        module procedure mf_gt_mf
        module procedure mf_gt_dp
        module procedure mf_gt_sp
        module procedure mf_gt_int
        module procedure dp_gt_mf
        module procedure sp_gt_mf
        module procedure int_gt_mf
    end interface
    public :: operator(>)

    interface operator(<=)
        module procedure mf_le_mf
        module procedure mf_le_dp
        module procedure mf_le_sp
        module procedure mf_le_int
        module procedure dp_le_mf
        module procedure sp_le_mf
        module procedure int_le_mf
    end interface
    public :: operator(<=)

    interface operator(>=)
        module procedure mf_ge_mf
        module procedure mf_ge_dp
        module procedure mf_ge_sp
        module procedure mf_ge_int
        module procedure dp_ge_mf
        module procedure sp_ge_mf
        module procedure int_ge_mf
    end interface
    public :: operator(>=)

    ! ================================================================
    ! Unary math generics
    ! ================================================================

    interface abs
        module procedure mf_abs
        module procedure cx_abs
    end interface
    public :: abs

    interface acos
        module procedure mf_acos
        module procedure cx_acos
    end interface
    public :: acos

    interface acosh
        module procedure mf_acosh
        module procedure cx_acosh
    end interface
    public :: acosh

    interface aimag
        module procedure cx_aimag
    end interface
    public :: aimag

    interface aint
        module procedure mf_aint
    end interface
    public :: aint

    interface anint
        module procedure mf_anint
    end interface
    public :: anint

    interface asin
        module procedure mf_asin
        module procedure cx_asin
    end interface
    public :: asin

    interface asinh
        module procedure mf_asinh
        module procedure cx_asinh
    end interface
    public :: asinh

    interface atan
        module procedure mf_atan
        module procedure cx_atan
    end interface
    public :: atan

    interface atanh
        module procedure mf_atanh
        module procedure cx_atanh
    end interface
    public :: atanh

    interface bessel_j0
        module procedure mf_bessel_j0
    end interface
    public :: bessel_j0

    interface bessel_j1
        module procedure mf_bessel_j1
    end interface
    public :: bessel_j1

    interface bessel_y0
        module procedure mf_bessel_y0
    end interface
    public :: bessel_y0

    interface bessel_y1
        module procedure mf_bessel_y1
    end interface
    public :: bessel_y1

    interface conjg
        module procedure cx_conjg
    end interface
    public :: conjg

    interface cos
        module procedure mf_cos
        module procedure cx_cos
    end interface
    public :: cos

    interface cosh
        module procedure mf_cosh
        module procedure cx_cosh
    end interface
    public :: cosh

    interface epsilon
        module procedure mf_epsilon
    end interface
    public :: epsilon

    interface erf
        module procedure mf_erf
    end interface
    public :: erf

    interface erfc
        module procedure mf_erfc
    end interface
    public :: erfc

    interface erfc_scaled
        module procedure mf_erfc_scaled
    end interface
    public :: erfc_scaled

    interface exp
        module procedure mf_exp
        module procedure cx_exp
    end interface
    public :: exp

    interface fraction
        module procedure mf_fraction
    end interface
    public :: fraction

    interface gamma
        module procedure mf_gamma
    end interface
    public :: gamma

    interface huge
        module procedure mf_huge
    end interface
    public :: huge

    interface log
        module procedure mf_log
        module procedure cx_log
    end interface
    public :: log

    interface log10
        module procedure mf_log10
    end interface
    public :: log10

    interface log_gamma
        module procedure mf_log_gamma
    end interface
    public :: log_gamma

    interface rrspacing
        module procedure mf_rrspacing
    end interface
    public :: rrspacing

    interface sin
        module procedure mf_sin
        module procedure cx_sin
    end interface
    public :: sin

    interface sinh
        module procedure mf_sinh
        module procedure cx_sinh
    end interface
    public :: sinh

    interface spacing
        module procedure mf_spacing
    end interface
    public :: spacing

    interface sqrt
        module procedure mf_sqrt
        module procedure cx_sqrt
    end interface
    public :: sqrt

    interface tan
        module procedure mf_tan
        module procedure cx_tan
    end interface
    public :: tan

    interface tanh
        module procedure mf_tanh
        module procedure cx_tanh
    end interface
    public :: tanh

    interface tiny
        module procedure mf_tiny
    end interface
    public :: tiny

    ! ================================================================
    ! Binary real functions (sign, mod, atan2, dim, modulo, hypot, nearest)
    ! ================================================================

    interface sign
        module procedure mf_sign_mf
        module procedure mf_sign_dp
        module procedure mf_sign_sp
        module procedure mf_sign_int
        module procedure dp_sign_mf
        module procedure sp_sign_mf
        module procedure int_sign_mf
    end interface
    public :: sign

    interface mod
        module procedure mf_mod_mf
        module procedure mf_mod_dp
        module procedure mf_mod_sp
        module procedure mf_mod_int
        module procedure dp_mod_mf
        module procedure sp_mod_mf
        module procedure int_mod_mf
    end interface
    public :: mod

    interface atan2
        module procedure mf_atan2_mf
        module procedure mf_atan2_dp
        module procedure mf_atan2_sp
        module procedure mf_atan2_int
        module procedure dp_atan2_mf
        module procedure sp_atan2_mf
        module procedure int_atan2_mf
    end interface
    public :: atan2

    interface dim
        module procedure mf_dim_mf
        module procedure mf_dim_dp
        module procedure mf_dim_sp
        module procedure mf_dim_int
        module procedure dp_dim_mf
        module procedure sp_dim_mf
        module procedure int_dim_mf
    end interface
    public :: dim

    interface modulo
        module procedure mf_modulo_mf
        module procedure mf_modulo_dp
        module procedure mf_modulo_sp
        module procedure mf_modulo_int
        module procedure dp_modulo_mf
        module procedure sp_modulo_mf
        module procedure int_modulo_mf
    end interface
    public :: modulo

    interface hypot
        module procedure mf_hypot_mf
        module procedure mf_hypot_dp
        module procedure mf_hypot_sp
        module procedure mf_hypot_int
        module procedure dp_hypot_mf
        module procedure sp_hypot_mf
        module procedure int_hypot_mf
    end interface
    public :: hypot

    interface nearest
        module procedure mf_nearest_mf
        module procedure mf_nearest_dp
        module procedure mf_nearest_sp
        module procedure mf_nearest_int
        module procedure dp_nearest_mf
        module procedure sp_nearest_mf
        module procedure int_nearest_mf
    end interface
    public :: nearest

    interface min
        module procedure mf_min_mf
        module procedure mf_min_dp
        module procedure mf_min_sp
        module procedure mf_min_int
        module procedure dp_min_mf
        module procedure sp_min_mf
        module procedure int_min_mf
        module procedure mf_min3
        module procedure mf_min4
        module procedure mf_min5
        module procedure mf_min6
        module procedure mf_min7
        module procedure mf_min8
    end interface
    public :: min

    interface max
        module procedure mf_max_mf
        module procedure mf_max_dp
        module procedure mf_max_sp
        module procedure mf_max_int
        module procedure dp_max_mf
        module procedure sp_max_mf
        module procedure int_max_mf
        module procedure mf_max3
        module procedure mf_max4
        module procedure mf_max5
        module procedure mf_max6
        module procedure mf_max7
        module procedure mf_max8
    end interface
    public :: max

    ! ================================================================
    ! Array reductions: value (maxval, minval)
    ! ================================================================

    interface maxval
        module procedure mf_maxval_1d
        module procedure mf_maxval_1d_dim
        module procedure mf_maxval_2d
        module procedure mf_maxval_2d_dim
        module procedure mf_maxval_3d
        module procedure mf_maxval_3d_dim
        module procedure mf_maxval_4d
        module procedure mf_maxval_4d_dim
        module procedure mf_maxval_5d
        module procedure mf_maxval_5d_dim
        module procedure mf_maxval_6d
        module procedure mf_maxval_6d_dim
        module procedure mf_maxval_7d
        module procedure mf_maxval_7d_dim
    end interface
    public :: maxval

    interface minval
        module procedure mf_minval_1d
        module procedure mf_minval_1d_dim
        module procedure mf_minval_2d
        module procedure mf_minval_2d_dim
        module procedure mf_minval_3d
        module procedure mf_minval_3d_dim
        module procedure mf_minval_4d
        module procedure mf_minval_4d_dim
        module procedure mf_minval_5d
        module procedure mf_minval_5d_dim
        module procedure mf_minval_6d
        module procedure mf_minval_6d_dim
        module procedure mf_minval_7d
        module procedure mf_minval_7d_dim
    end interface
    public :: minval

    ! ================================================================
    ! Array reductions: location (maxloc, minloc)
    ! ================================================================

    interface maxloc
        module procedure mf_maxloc_1d
        module procedure mf_maxloc_1d_dim
        module procedure mf_maxloc_2d
        module procedure mf_maxloc_2d_dim
        module procedure mf_maxloc_3d
        module procedure mf_maxloc_3d_dim
        module procedure mf_maxloc_4d
        module procedure mf_maxloc_4d_dim
        module procedure mf_maxloc_5d
        module procedure mf_maxloc_5d_dim
        module procedure mf_maxloc_6d
        module procedure mf_maxloc_6d_dim
        module procedure mf_maxloc_7d
        module procedure mf_maxloc_7d_dim
    end interface
    public :: maxloc

    interface minloc
        module procedure mf_minloc_1d
        module procedure mf_minloc_1d_dim
        module procedure mf_minloc_2d
        module procedure mf_minloc_2d_dim
        module procedure mf_minloc_3d
        module procedure mf_minloc_3d_dim
        module procedure mf_minloc_4d
        module procedure mf_minloc_4d_dim
        module procedure mf_minloc_5d
        module procedure mf_minloc_5d_dim
        module procedure mf_minloc_6d
        module procedure mf_minloc_6d_dim
        module procedure mf_minloc_7d
        module procedure mf_minloc_7d_dim
    end interface
    public :: minloc

    ! ================================================================
    ! Array reductions: sum, product (mf + cx)
    ! ================================================================

    interface sum
        module procedure mf_sum_1d
        module procedure mf_sum_1d_dim
        module procedure cx_sum_1d
        module procedure cx_sum_1d_dim
        module procedure mf_sum_2d
        module procedure mf_sum_2d_dim
        module procedure cx_sum_2d
        module procedure cx_sum_2d_dim
        module procedure mf_sum_3d
        module procedure mf_sum_3d_dim
        module procedure cx_sum_3d
        module procedure cx_sum_3d_dim
        module procedure mf_sum_4d
        module procedure mf_sum_4d_dim
        module procedure cx_sum_4d
        module procedure cx_sum_4d_dim
        module procedure mf_sum_5d
        module procedure mf_sum_5d_dim
        module procedure cx_sum_5d
        module procedure cx_sum_5d_dim
        module procedure mf_sum_6d
        module procedure mf_sum_6d_dim
        module procedure cx_sum_6d
        module procedure cx_sum_6d_dim
        module procedure mf_sum_7d
        module procedure mf_sum_7d_dim
        module procedure cx_sum_7d
        module procedure cx_sum_7d_dim
    end interface
    public :: sum

    interface product
        module procedure mf_product_1d
        module procedure mf_product_1d_dim
        module procedure cx_product_1d
        module procedure cx_product_1d_dim
        module procedure mf_product_2d
        module procedure mf_product_2d_dim
        module procedure cx_product_2d
        module procedure cx_product_2d_dim
        module procedure mf_product_3d
        module procedure mf_product_3d_dim
        module procedure cx_product_3d
        module procedure cx_product_3d_dim
        module procedure mf_product_4d
        module procedure mf_product_4d_dim
        module procedure cx_product_4d
        module procedure cx_product_4d_dim
        module procedure mf_product_5d
        module procedure mf_product_5d_dim
        module procedure cx_product_5d
        module procedure cx_product_5d_dim
        module procedure mf_product_6d
        module procedure mf_product_6d_dim
        module procedure cx_product_6d
        module procedure cx_product_6d_dim
        module procedure mf_product_7d
        module procedure mf_product_7d_dim
        module procedure cx_product_7d
        module procedure cx_product_7d_dim
    end interface
    public :: product

    ! ================================================================
    ! dot_product, norm2, findloc, matmul
    ! ================================================================

    interface dot_product
        module procedure mf_dot_product
        module procedure cx_dot_product
    end interface
    public :: dot_product

    interface norm2
        module procedure mf_norm2_1d
        module procedure mf_norm2_2d
        module procedure mf_norm2_2d_dim
        module procedure mf_norm2_3d
        module procedure mf_norm2_3d_dim
        module procedure mf_norm2_4d
        module procedure mf_norm2_4d_dim
        module procedure mf_norm2_5d
        module procedure mf_norm2_5d_dim
        module procedure mf_norm2_6d
        module procedure mf_norm2_6d_dim
        module procedure mf_norm2_7d
        module procedure mf_norm2_7d_dim
    end interface
    public :: norm2

    interface findloc
        module procedure mf_findloc_1d
        module procedure mf_findloc_1d_dim
        module procedure cx_findloc_1d
        module procedure cx_findloc_1d_dim
        module procedure mf_findloc_2d
        module procedure mf_findloc_2d_dim
        module procedure cx_findloc_2d
        module procedure cx_findloc_2d_dim
        module procedure mf_findloc_3d
        module procedure mf_findloc_3d_dim
        module procedure cx_findloc_3d
        module procedure cx_findloc_3d_dim
        module procedure mf_findloc_4d
        module procedure mf_findloc_4d_dim
        module procedure cx_findloc_4d
        module procedure cx_findloc_4d_dim
        module procedure mf_findloc_5d
        module procedure mf_findloc_5d_dim
        module procedure cx_findloc_5d
        module procedure cx_findloc_5d_dim
        module procedure mf_findloc_6d
        module procedure mf_findloc_6d_dim
        module procedure cx_findloc_6d
        module procedure cx_findloc_6d_dim
        module procedure mf_findloc_7d
        module procedure mf_findloc_7d_dim
        module procedure cx_findloc_7d
        module procedure cx_findloc_7d_dim
    end interface
    public :: findloc

    interface matmul
        module procedure mf_matmul_mm
        module procedure mf_matmul_mv
        module procedure mf_matmul_vm
        module procedure cx_matmul_mm
        module procedure cx_matmul_mv
        module procedure cx_matmul_vm
    end interface
    public :: matmul

    ! ================================================================
    ! Bessel functions with integer order
    ! ================================================================

    interface bessel_jn
        module procedure mf_bessel_jn_elem
        module procedure mf_bessel_jn_range
    end interface
    public :: bessel_jn

    interface bessel_yn
        module procedure mf_bessel_yn_elem
        module procedure mf_bessel_yn_range
    end interface
    public :: bessel_yn

    ! ================================================================
    ! Conversions and numeric functions
    ! ================================================================

    interface ceiling
        module procedure mf_ceiling
    end interface
    public :: ceiling

    interface floor
        module procedure mf_floor
    end interface
    public :: floor

    interface exponent
        module procedure mf_exponent
    end interface
    public :: exponent

    interface scale
        module procedure mf_scale
    end interface
    public :: scale

    interface set_exponent
        module procedure mf_set_exponent
    end interface
    public :: set_exponent

    ! ================================================================
    ! Numeric inquiry functions
    ! ================================================================

    interface digits
        module procedure mf_digits
    end interface
    public :: digits

    interface maxexponent
        module procedure mf_maxexponent
    end interface
    public :: maxexponent

    interface minexponent
        module procedure mf_minexponent
    end interface
    public :: minexponent

    interface radix
        module procedure mf_radix
    end interface
    public :: radix

    interface precision
        module procedure mf_precision
        module procedure cx_precision
    end interface
    public :: precision

    interface range
        module procedure mf_range
        module procedure cx_range
    end interface
    public :: range

    interface storage_size
        module procedure mf_storage_size
        module procedure cx_storage_size
    end interface
    public :: storage_size

    ! ================================================================
    ! Complex-specific: cmplx, real, mf_real
    ! ================================================================

    interface cmplx
        module procedure mf_cmplx_2
        module procedure mf_cmplx_1
    end interface
    public :: cmplx

    interface real
        module procedure cx_real
        module procedure mf_real_mf
    end interface
    public :: real

    interface mf_real
        module procedure mf_real_mf
        module procedure cx_real_mf
    end interface
    public :: mf_real

    ! ================================================================
    ! Conversions: dble, int, nint, mf_to_double
    ! ================================================================

    interface dble
        module procedure mf_dble
    end interface
    public :: dble

    interface int
        module procedure mf_int
    end interface
    public :: int

    interface nint
        module procedure mf_nint
    end interface
    public :: nint

    public :: mf_to_double

    ! ================================================================
    ! random_number (impure subroutine)
    ! ================================================================

    interface random_number
        module procedure mf_random_number_0d
        module procedure mf_random_number_1d
        module procedure mf_random_number_2d
        module procedure mf_random_number_3d
        module procedure mf_random_number_4d
        module procedure mf_random_number_5d
        module procedure mf_random_number_6d
        module procedure mf_random_number_7d
    end interface
    public :: random_number

    ! ================================================================
    ! Assignment (=)
    ! ================================================================

    interface assignment(=)
        module procedure mf_assign_dp
        module procedure mf_assign_sp
        module procedure mf_assign_int
        module procedure mf_assign_cx
        module procedure mf_assign_cdp
        module procedure mf_assign_csp
        module procedure cx_assign_mf
        module procedure cx_assign_dp
        module procedure cx_assign_sp
        module procedure cx_assign_int
        module procedure cx_assign_cdp
        module procedure cx_assign_csp
    end interface
    public :: assignment(=)

    ! ================================================================
    ! Defined I/O — intentionally OMITTED.
    ! ================================================================
    !
    ! float64x2 and complex128x2 must carry the SEQUENCE attribute so
    ! that LAPACK routines such as DLALN2 (wlaln2) can EQUIVALENCE
    ! their work arrays. gfortran (as of 15.2) treats any module that
    ! defines DTIO procedures for a derived type as displacing the
    ! SEQUENCE attribute on use-association — the use-associated type
    ! shows up with HAS-DTIO-PROCS instead, and any EQUIVALENCE of
    ! such a variable is then rejected. We therefore drop the DTIO
    ! generics entirely. Tests that need to print a value should fall
    ! back to component-wise I/O (e.g. ``write(*, *) x%limbs``).

    ! ================================================================
    ! Named constants
    ! ================================================================

    type(float64x2), parameter, public :: MF_ZERO  = float64x2([0.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_ONE   = float64x2([1.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_TWO   = float64x2([2.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_HALF  = float64x2([0.5_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_EIGHT = float64x2([8.0_dp, 0.0_dp])

    type(float64x2), parameter, public :: MF_SAFMIN = float64x2([tiny(0.0_dp),  0.0_dp])
    type(float64x2), parameter, public :: MF_SAFMAX = float64x2([huge(0.0_dp),  0.0_dp])
    type(float64x2), parameter, public :: MF_TSML   = float64x2([1.0e-100_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_TBIG   = float64x2([1.0e+100_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_SSML   = float64x2([1.0e-50_dp,    0.0_dp])
    type(float64x2), parameter, public :: MF_SBIG   = float64x2([1.0e+50_dp,    0.0_dp])
    type(float64x2), parameter, public :: MF_RTMIN  = float64x2([1.0e-150_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_RTMAX  = float64x2([1.0e+150_dp,   0.0_dp])

contains

    ! ================================================================
    ! Internal conversion helpers (elemental)
    ! ================================================================

    elemental function to_qp(mf) result(res)
        type(float64x2), intent(in) :: mf
        real(qp) :: res
        res = real(mf%limbs(1), qp) + real(mf%limbs(2), qp)
    end function

    elemental function from_qp(v) result(res)
        real(qp), intent(in) :: v
        type(float64x2) :: res
        res%limbs(1) = real(v, dp)
        res%limbs(2) = real(v - real(res%limbs(1), qp), dp)
    end function

    elemental function to_cqp(cx) result(res)
        type(complex128x2), intent(in) :: cx
        complex(qp) :: res
        res = cmplx(to_qp(cx%re), to_qp(cx%im), kind=qp)
    end function

    elemental function from_cqp(v) result(res)
        complex(qp), intent(in) :: v
        type(complex128x2) :: res
        res%re = from_qp(real(v, qp))
        res%im = from_qp(aimag(v))
    end function

    ! ================================================================
    ! Constructors
    ! ================================================================

    elemental function mf_from_dp(v) result(res)
        real(dp), intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    elemental function mf_from_sp(v) result(res)
        real(sp), intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    elemental function mf_from_int(v) result(res)
        integer, intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    function mf_from_char(v) result(res)
        character(len=*), intent(in) :: v
        type(float64x2) :: res
        real(qp) :: val_qp
        read(v, *) val_qp
        res = from_qp(val_qp)
    end function

    elemental function mf_from_cx(cx) result(res)
        type(complex128x2), intent(in) :: cx
        type(float64x2) :: res
        res = cx%re
    end function

    elemental function mf_from_cdp(x) result(res)
        complex(dp), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(x, qp))
    end function

    elemental function mf_from_csp(x) result(res)
        complex(sp), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(x, qp))
    end function

    elemental function cx_from_mf_2(r, i) result(res)
        type(float64x2), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = r
        res%im = i
    end function

    elemental function cx_from_dp_2(r, i) result(res)
        real(dp), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = from_qp(real(r, qp))
        res%im = from_qp(real(i, qp))
    end function

    elemental function cx_from_mf_1(r) result(res)
        type(float64x2), intent(in) :: r
        type(complex128x2) :: res
        res%re = r
        res%im = from_qp(0.0_qp)
    end function

    elemental function cx_from_dp_1(r) result(res)
        real(dp), intent(in) :: r
        type(complex128x2) :: res
        res%re = from_qp(real(r, qp))
        res%im = from_qp(0.0_qp)
    end function

    ! ================================================================
    ! Assignment (elemental)
    ! ================================================================

    elemental subroutine mf_assign_dp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    elemental subroutine mf_assign_sp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(sp), intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    elemental subroutine mf_assign_int(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    elemental subroutine mf_assign_cx(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        type(complex128x2), intent(in) :: rhs
        lhs = from_qp(real(to_cqp(rhs), qp))
    end subroutine

    elemental subroutine mf_assign_cdp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        complex(dp), intent(in) :: rhs
        lhs = from_qp(real(cmplx(rhs, kind=qp), qp))
    end subroutine

    elemental subroutine mf_assign_csp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        complex(sp), intent(in) :: rhs
        lhs = from_qp(real(cmplx(rhs, kind=qp), qp))
    end subroutine

    elemental subroutine cx_assign_mf(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        type(float64x2), intent(in) :: rhs
        lhs = from_cqp(cmplx(to_qp(rhs), 0.0_qp, qp))
    end subroutine

    elemental subroutine cx_assign_dp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    elemental subroutine cx_assign_sp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        real(sp), intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    elemental subroutine cx_assign_int(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    elemental subroutine cx_assign_cdp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        complex(dp), intent(in) :: rhs
        lhs = from_cqp(cmplx(rhs, kind=qp))
    end subroutine

    elemental subroutine cx_assign_csp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        complex(sp), intent(in) :: rhs
        lhs = from_cqp(cmplx(rhs, kind=qp))
    end subroutine

    ! ================================================================
    ! Conversions (elemental)
    ! ================================================================

    elemental function mf_dble(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(to_qp(x), dp)
    end function

    elemental function mf_int(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = int(to_qp(x))
    end function

    elemental function mf_nint(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = nint(to_qp(x))
    end function

    elemental function mf_to_double(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(to_qp(x), dp)
    end function

    elemental function mf_ceiling(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = ceiling(to_qp(x))
    end function

    elemental function mf_floor(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = floor(to_qp(x))
    end function

    elemental function mf_exponent(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = exponent(to_qp(x))
    end function

    ! ================================================================
    ! Unary operators (elemental)
    ! ================================================================

    elemental function mf_neg(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(-to_qp(x))
    end function

    elemental function cx_neg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(-to_cqp(x))
    end function

    ! ================================================================
    ! Unary math: mf → mf (elemental)
    ! ================================================================

    elemental function mf_abs(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(abs(to_qp(x)))
    end function

    elemental function mf_sqrt(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(sqrt(to_qp(x)))
    end function

    elemental function mf_sin(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(sin(to_qp(x)))
    end function

    elemental function mf_cos(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(cos(to_qp(x)))
    end function

    elemental function mf_tan(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(tan(to_qp(x)))
    end function

    elemental function mf_exp(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(exp(to_qp(x)))
    end function

    elemental function mf_log(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(log(to_qp(x)))
    end function

    elemental function mf_log10(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(log10(to_qp(x)))
    end function

    elemental function mf_atan(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(atan(to_qp(x)))
    end function

    elemental function mf_asin(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(asin(to_qp(x)))
    end function

    elemental function mf_acos(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(acos(to_qp(x)))
    end function

    elemental function mf_aint(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(aint(to_qp(x)))
    end function

    elemental function mf_anint(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(anint(to_qp(x)))
    end function

    elemental function mf_sinh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(sinh(to_qp(x)))
    end function

    elemental function mf_cosh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(cosh(to_qp(x)))
    end function

    elemental function mf_tanh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(tanh(to_qp(x)))
    end function

    elemental function mf_asinh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(asinh(to_qp(x)))
    end function

    elemental function mf_acosh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(acosh(to_qp(x)))
    end function

    elemental function mf_atanh(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(atanh(to_qp(x)))
    end function

    elemental function mf_erf(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(erf(to_qp(x)))
    end function

    elemental function mf_erfc(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(erfc(to_qp(x)))
    end function

    elemental function mf_erfc_scaled(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(erfc_scaled(to_qp(x)))
    end function

    elemental function mf_gamma(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(gamma(to_qp(x)))
    end function

    elemental function mf_log_gamma(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(log_gamma(to_qp(x)))
    end function

    elemental function mf_bessel_j0(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_j0(to_qp(x)))
    end function

    elemental function mf_bessel_j1(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_j1(to_qp(x)))
    end function

    elemental function mf_bessel_y0(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_y0(to_qp(x)))
    end function

    elemental function mf_bessel_y1(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_y1(to_qp(x)))
    end function

    elemental function mf_fraction(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(fraction(to_qp(x)))
    end function

    elemental function mf_rrspacing(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(rrspacing(to_qp(x)))
    end function

    elemental function mf_spacing(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(spacing(to_qp(x)))
    end function

    elemental function mf_epsilon(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(epsilon(to_qp(x)))
    end function

    elemental function mf_huge(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(huge(to_qp(x)))
    end function

    elemental function mf_tiny(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(tiny(to_qp(x)))
    end function

    ! ================================================================
    ! Unary math: cx → cx (elemental)
    ! ================================================================

    elemental function cx_sqrt(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sqrt(to_cqp(x)))
    end function

    elemental function cx_sin(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sin(to_cqp(x)))
    end function

    elemental function cx_cos(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(cos(to_cqp(x)))
    end function

    elemental function cx_tan(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(tan(to_cqp(x)))
    end function

    elemental function cx_exp(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(exp(to_cqp(x)))
    end function

    elemental function cx_log(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(log(to_cqp(x)))
    end function

    elemental function cx_atan(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(atan(to_cqp(x)))
    end function

    elemental function cx_asin(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(asin(to_cqp(x)))
    end function

    elemental function cx_acos(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(acos(to_cqp(x)))
    end function

    elemental function cx_conjg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(conjg(to_cqp(x)))
    end function

    elemental function cx_sinh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sinh(to_cqp(x)))
    end function

    elemental function cx_cosh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(cosh(to_cqp(x)))
    end function

    elemental function cx_tanh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(tanh(to_cqp(x)))
    end function

    elemental function cx_asinh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(asinh(to_cqp(x)))
    end function

    elemental function cx_acosh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(acosh(to_cqp(x)))
    end function

    elemental function cx_atanh(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(atanh(to_cqp(x)))
    end function

    ! ================================================================
    ! Unary math: cx → mf (elemental)
    ! ================================================================

    elemental function cx_abs(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(abs(to_cqp(x)))
    end function

    elemental function cx_aimag(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(aimag(to_cqp(x)))
    end function

    ! ================================================================
    ! Complex-specific: cmplx, real, mf_real (elemental)
    ! ================================================================

    elemental function mf_cmplx_2(re, im) result(res)
        type(float64x2), intent(in) :: re, im
        type(complex128x2) :: res
        res%re = re
        res%im = im
    end function

    elemental function mf_cmplx_1(re) result(res)
        type(float64x2), intent(in) :: re
        type(complex128x2) :: res
        res%re = re
        res%im = from_qp(0.0_qp)
    end function

    elemental function cx_real(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(to_cqp(x), qp))
    end function

    elemental function mf_real_mf(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = x
    end function

    elemental function cx_real_mf(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = x%re
    end function

    ! ================================================================
    ! Binary real functions (elemental)
    ! ================================================================

    elemental function mf_sign_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), to_qp(b)))
    end function

    elemental function mf_sign_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    elemental function mf_sign_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    elemental function mf_sign_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    elemental function dp_sign_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    elemental function sp_sign_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    elemental function int_sign_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    elemental function mf_mod_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), to_qp(b)))
    end function

    elemental function mf_mod_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    elemental function mf_mod_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    elemental function mf_mod_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    elemental function dp_mod_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    elemental function sp_mod_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    elemental function int_mod_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    elemental function mf_atan2_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), to_qp(b)))
    end function

    elemental function mf_atan2_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    elemental function mf_atan2_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    elemental function mf_atan2_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    elemental function dp_atan2_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    elemental function sp_atan2_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    elemental function int_atan2_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    elemental function mf_dim_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(to_qp(a), to_qp(b)))
    end function

    elemental function mf_dim_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(to_qp(a), real(b, qp)))
    end function

    elemental function mf_dim_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(to_qp(a), real(b, qp)))
    end function

    elemental function mf_dim_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(to_qp(a), real(b, qp)))
    end function

    elemental function dp_dim_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(real(a, qp), to_qp(b)))
    end function

    elemental function sp_dim_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(real(a, qp), to_qp(b)))
    end function

    elemental function int_dim_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(dim(real(a, qp), to_qp(b)))
    end function

    elemental function mf_modulo_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(to_qp(a), to_qp(b)))
    end function

    elemental function mf_modulo_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(to_qp(a), real(b, qp)))
    end function

    elemental function mf_modulo_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(to_qp(a), real(b, qp)))
    end function

    elemental function mf_modulo_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(to_qp(a), real(b, qp)))
    end function

    elemental function dp_modulo_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(real(a, qp), to_qp(b)))
    end function

    elemental function sp_modulo_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(real(a, qp), to_qp(b)))
    end function

    elemental function int_modulo_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(modulo(real(a, qp), to_qp(b)))
    end function

    elemental function mf_hypot_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(to_qp(a), to_qp(b)))
    end function

    elemental function mf_hypot_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(to_qp(a), real(b, qp)))
    end function

    elemental function mf_hypot_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(to_qp(a), real(b, qp)))
    end function

    elemental function mf_hypot_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(to_qp(a), real(b, qp)))
    end function

    elemental function dp_hypot_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(real(a, qp), to_qp(b)))
    end function

    elemental function sp_hypot_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(real(a, qp), to_qp(b)))
    end function

    elemental function int_hypot_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(hypot(real(a, qp), to_qp(b)))
    end function

    elemental function mf_nearest_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(to_qp(a), to_qp(b)))
    end function

    elemental function mf_nearest_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(to_qp(a), real(b, qp)))
    end function

    elemental function mf_nearest_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(to_qp(a), real(b, qp)))
    end function

    elemental function mf_nearest_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(to_qp(a), real(b, qp)))
    end function

    elemental function dp_nearest_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(real(a, qp), to_qp(b)))
    end function

    elemental function sp_nearest_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(real(a, qp), to_qp(b)))
    end function

    elemental function int_nearest_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(nearest(real(a, qp), to_qp(b)))
    end function

    elemental function mf_min_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), to_qp(b)))
    end function

    elemental function mf_min_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    elemental function mf_min_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    elemental function mf_min_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    elemental function dp_min_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    elemental function sp_min_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    elemental function int_min_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    elemental function mf_max_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), to_qp(b)))
    end function

    elemental function mf_max_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    elemental function mf_max_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    elemental function mf_max_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    elemental function dp_max_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    elemental function sp_max_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    elemental function int_max_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    ! min/max with N arguments (N = 3..8). LAPACK has call sites with
    ! up to 5 arguments; the upper bound here is generous so future
    ! call sites do not regress.
    elemental function mf_min3(a1, a2, a3) result(res)
        type(float64x2), intent(in) :: a1, a2, a3
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3)))
    end function

    elemental function mf_min4(a1, a2, a3, a4) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4)))
    end function

    elemental function mf_min5(a1, a2, a3, a4, a5) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5)))
    end function

    elemental function mf_min6(a1, a2, a3, a4, a5, a6) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6)))
    end function

    elemental function mf_min7(a1, a2, a3, a4, a5, a6, a7) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6, a7
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6), to_qp(a7)))
    end function

    elemental function mf_min8(a1, a2, a3, a4, a5, a6, a7, a8) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8
        type(float64x2) :: res
        res = from_qp(min(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6), to_qp(a7), to_qp(a8)))
    end function

    elemental function mf_max3(a1, a2, a3) result(res)
        type(float64x2), intent(in) :: a1, a2, a3
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3)))
    end function

    elemental function mf_max4(a1, a2, a3, a4) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4)))
    end function

    elemental function mf_max5(a1, a2, a3, a4, a5) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5)))
    end function

    elemental function mf_max6(a1, a2, a3, a4, a5, a6) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6)))
    end function

    elemental function mf_max7(a1, a2, a3, a4, a5, a6, a7) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6, a7
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6), to_qp(a7)))
    end function

    elemental function mf_max8(a1, a2, a3, a4, a5, a6, a7, a8) result(res)
        type(float64x2), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8
        type(float64x2) :: res
        res = from_qp(max(to_qp(a1), to_qp(a2), to_qp(a3), to_qp(a4), to_qp(a5), to_qp(a6), to_qp(a7), to_qp(a8)))
    end function

    ! ================================================================
    ! Bessel functions with integer order (elemental + transformational)
    ! ================================================================

    elemental function mf_bessel_jn_elem(n, x) result(res)
        integer, intent(in) :: n
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_jn(n, to_qp(x)))
    end function

    pure function mf_bessel_jn_range(n1, n2, x) result(res)
        integer, intent(in) :: n1, n2
        type(float64x2), intent(in) :: x
        type(float64x2) :: res(max(n2 - n1 + 1, 0))
        real(qp) :: tmp(max(n2 - n1 + 1, 0))
        integer :: i
        tmp = bessel_jn(n1, n2, to_qp(x))
        do i = 1, size(tmp)
            res(i) = from_qp(tmp(i))
        end do
    end function

    elemental function mf_bessel_yn_elem(n, x) result(res)
        integer, intent(in) :: n
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(bessel_yn(n, to_qp(x)))
    end function

    pure function mf_bessel_yn_range(n1, n2, x) result(res)
        integer, intent(in) :: n1, n2
        type(float64x2), intent(in) :: x
        type(float64x2) :: res(max(n2 - n1 + 1, 0))
        real(qp) :: tmp(max(n2 - n1 + 1, 0))
        integer :: i
        tmp = bessel_yn(n1, n2, to_qp(x))
        do i = 1, size(tmp)
            res(i) = from_qp(tmp(i))
        end do
    end function

    ! ================================================================
    ! Functions: (mf, integer) → mf (elemental)
    ! ================================================================

    elemental function mf_scale(x, i) result(res)
        type(float64x2), intent(in) :: x
        integer, intent(in) :: i
        type(float64x2) :: res
        res = from_qp(scale(to_qp(x), i))
    end function

    elemental function mf_set_exponent(x, i) result(res)
        type(float64x2), intent(in) :: x
        integer, intent(in) :: i
        type(float64x2) :: res
        res = from_qp(set_exponent(to_qp(x), i))
    end function

    ! ================================================================
    ! Numeric inquiry functions (elemental)
    ! ================================================================

    elemental function mf_digits(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = digits(0.0_qp)
    end function

    elemental function mf_maxexponent(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = maxexponent(0.0_qp)
    end function

    elemental function mf_minexponent(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = minexponent(0.0_qp)
    end function

    elemental function mf_radix(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = radix(0.0_qp)
    end function

    elemental function mf_precision(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = precision(0.0_qp)
    end function

    elemental function cx_precision(x) result(res)
        type(complex128x2), intent(in) :: x
        integer :: res
        res = precision(0.0_qp)
    end function

    elemental function mf_range(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = range(0.0_qp)
    end function

    elemental function cx_range(x) result(res)
        type(complex128x2), intent(in) :: x
        integer :: res
        res = range(0.0_qp)
    end function

    elemental function mf_storage_size(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = storage_size(0.0_qp)
    end function

    elemental function cx_storage_size(x) result(res)
        type(complex128x2), intent(in) :: x
        integer :: res
        res = storage_size(cmplx(0.0_qp, 0.0_qp, qp))
    end function

    ! ================================================================
    ! Array reductions: value (maxval, minval) — rank 1-7
    ! ================================================================

    pure function mf_maxval_1d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_1d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, dim=dim, mask=mask))
        else
            res = from_qp(maxval(tmp, dim=dim))
        end if
    end function

    pure function mf_maxval_2d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_2d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_qp(tmp_res(i1))
        end do
    end function

    pure function mf_maxval_3d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_3d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_qp(tmp_res(i1, i2))
        end do
        end do
    end function

    pure function mf_maxval_4d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_4d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_qp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    pure function mf_maxval_5d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_5d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_qp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    pure function mf_maxval_6d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_6d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_qp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    pure function mf_maxval_7d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(maxval(tmp, mask=mask))
        else
            res = from_qp(maxval(tmp))
        end if
    end function

    pure function mf_maxval_7d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = maxval(tmp, dim=dim, mask=mask)
        else
            tmp_res = maxval(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_qp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    pure function mf_minval_1d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_1d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, dim=dim, mask=mask))
        else
            res = from_qp(minval(tmp, dim=dim))
        end if
    end function

    pure function mf_minval_2d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_2d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_qp(tmp_res(i1))
        end do
    end function

    pure function mf_minval_3d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_3d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_qp(tmp_res(i1, i2))
        end do
        end do
    end function

    pure function mf_minval_4d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_4d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_qp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    pure function mf_minval_5d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_5d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_qp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    pure function mf_minval_6d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_6d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_qp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    pure function mf_minval_7d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(minval(tmp, mask=mask))
        else
            res = from_qp(minval(tmp))
        end if
    end function

    pure function mf_minval_7d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = minval(tmp, dim=dim, mask=mask)
        else
            tmp_res = minval(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_qp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! ================================================================
    ! Array reductions: location (maxloc, minloc) — rank 1-7
    ! ================================================================

    pure function mf_maxloc_1d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res(1)
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_1d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_2d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_2d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_3d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(3)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_3d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_4d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(4)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_4d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_5d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(5)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_5d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_6d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(6)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_6d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_maxloc_7d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(7)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, back=back)
        else
            res = maxloc(tmp)
        end if
    end function

    pure function mf_maxloc_7d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = maxloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = maxloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = maxloc(tmp, dim=dim, back=back)
        else
            res = maxloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_1d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res(1)
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_1d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_2d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_2d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_3d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(3)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_3d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_4d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(4)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_4d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_5d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(5)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_5d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_6d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(6)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_6d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    pure function mf_minloc_7d(array, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(7)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, back=back)
        else
            res = minloc(tmp)
        end if
    end function

    pure function mf_minloc_7d_dim(array, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = minloc(tmp, dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = minloc(tmp, dim=dim, mask=mask)
        else if (present(back)) then
            res = minloc(tmp, dim=dim, back=back)
        else
            res = minloc(tmp, dim=dim)
        end if
    end function

    ! ================================================================
    ! Array reductions: sum, product — rank 1-7 (mf + cx)
    ! ================================================================

    ! mf sum rank 1
    pure function mf_sum_1d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_1d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, dim=dim, mask=mask))
        else
            res = from_qp(sum(tmp, dim=dim))
        end if
    end function

    ! cx sum rank 1
    pure function cx_sum_1d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_1d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, dim=dim, mask=mask))
        else
            res = from_cqp(sum(tmp, dim=dim))
        end if
    end function

    ! mf sum rank 2
    pure function mf_sum_2d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_2d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_qp(tmp_res(i1))
        end do
    end function

    ! cx sum rank 2
    pure function cx_sum_2d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_2d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        complex(qp) :: tmp(size(array,1), size(array,2))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_cqp(tmp_res(i1))
        end do
    end function

    ! mf sum rank 3
    pure function mf_sum_3d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_3d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_qp(tmp_res(i1, i2))
        end do
        end do
    end function

    ! cx sum rank 3
    pure function cx_sum_3d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_3d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_cqp(tmp_res(i1, i2))
        end do
        end do
    end function

    ! mf sum rank 4
    pure function mf_sum_4d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_4d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_qp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    ! cx sum rank 4
    pure function cx_sum_4d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_4d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_cqp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    ! mf sum rank 5
    pure function mf_sum_5d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_5d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_qp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    ! cx sum rank 5
    pure function cx_sum_5d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_5d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_cqp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    ! mf sum rank 6
    pure function mf_sum_6d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_6d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_qp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    ! cx sum rank 6
    pure function cx_sum_6d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_6d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_cqp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    ! mf sum rank 7
    pure function mf_sum_7d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(sum(tmp, mask=mask))
        else
            res = from_qp(sum(tmp))
        end if
    end function

    pure function mf_sum_7d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_qp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! cx sum rank 7
    pure function cx_sum_7d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(sum(tmp, mask=mask))
        else
            res = from_cqp(sum(tmp))
        end if
    end function

    pure function cx_sum_7d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = sum(tmp, dim=dim, mask=mask)
        else
            tmp_res = sum(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_cqp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! mf product rank 1
    pure function mf_product_1d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_1d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, dim=dim, mask=mask))
        else
            res = from_qp(product(tmp, dim=dim))
        end if
    end function

    ! cx product rank 1
    pure function cx_product_1d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:)
        logical, intent(in), optional :: mask(:)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_1d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, dim=dim, mask=mask))
        else
            res = from_cqp(product(tmp, dim=dim))
        end if
    end function

    ! mf product rank 2
    pure function mf_product_2d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_2d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_qp(tmp_res(i1))
        end do
    end function

    ! cx product rank 2
    pure function cx_product_2d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        logical, intent(in), optional :: mask(:, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_2d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        complex(qp) :: tmp(size(array,1), size(array,2))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_cqp(tmp_res(i1))
        end do
    end function

    ! mf product rank 3
    pure function mf_product_3d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_3d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_qp(tmp_res(i1, i2))
        end do
        end do
    end function

    ! cx product rank 3
    pure function cx_product_3d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        logical, intent(in), optional :: mask(:, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_3d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_cqp(tmp_res(i1, i2))
        end do
        end do
    end function

    ! mf product rank 4
    pure function mf_product_4d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_4d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_qp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    ! cx product rank 4
    pure function cx_product_4d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_4d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_cqp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    ! mf product rank 5
    pure function mf_product_5d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_5d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_qp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    ! cx product rank 5
    pure function cx_product_5d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_5d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_cqp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    ! mf product rank 6
    pure function mf_product_6d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_6d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_qp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    ! cx product rank 6
    pure function cx_product_6d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_6d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_cqp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    ! mf product rank 7
    pure function mf_product_7d(array, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_qp(product(tmp, mask=mask))
        else
            res = from_qp(product(tmp))
        end if
    end function

    pure function mf_product_7d_dim(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_qp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! cx product rank 7
    pure function cx_product_7d(array, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(complex128x2) :: res
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            res = from_cqp(product(tmp, mask=mask))
        else
            res = from_cqp(product(tmp))
        end if
    end function

    pure function cx_product_7d_dim(array, dim, mask) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        type(complex128x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        complex(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask)) then
            tmp_res = product(tmp, dim=dim, mask=mask)
        else
            tmp_res = product(tmp, dim=dim)
        end if
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_cqp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! ================================================================
    ! dot_product (pure)
    ! ================================================================

    pure function mf_dot_product(a, b) result(res)
        type(float64x2), intent(in) :: a(:), b(:)
        type(float64x2) :: res
        real(qp) :: qa(size(a)), qb(size(b))
        integer :: i1
        do i1 = 1, size(a)
            qa(i1) = to_qp(a(i1))
        end do
        do i1 = 1, size(b)
            qb(i1) = to_qp(b(i1))
        end do
        res = from_qp(dot_product(qa, qb))
    end function

    pure function cx_dot_product(a, b) result(res)
        type(complex128x2), intent(in) :: a(:), b(:)
        type(complex128x2) :: res
        complex(qp) :: qa(size(a)), qb(size(b))
        integer :: i1
        do i1 = 1, size(a)
            qa(i1) = to_cqp(a(i1))
        end do
        do i1 = 1, size(b)
            qb(i1) = to_cqp(b(i1))
        end do
        res = from_cqp(dot_product(qa, qb))
    end function

    ! ================================================================
    ! norm2 (pure) — rank 1-7, NO mask per standard
    ! ================================================================

    pure function mf_norm2_1d(array) result(res)
        type(float64x2), intent(in) :: array(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        res = from_qp(norm2(tmp))
    end function


    pure function mf_norm2_2d(array) result(res)
        type(float64x2), intent(in) :: array(:, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_2d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i1 = 1, size(tmp_res, 1)
            res(i1) = from_qp(tmp_res(i1))
        end do
    end function

    pure function mf_norm2_3d(array) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_3d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2) = from_qp(tmp_res(i1, i2))
        end do
        end do
    end function

    pure function mf_norm2_4d(array) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_4d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3) = from_qp(tmp_res(i1, i2, i3))
        end do
        end do
        end do
    end function

    pure function mf_norm2_5d(array) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_5d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4) = from_qp(tmp_res(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end function

    pure function mf_norm2_6d(array) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_6d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5) = from_qp(tmp_res(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end function

    pure function mf_norm2_7d(array) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        res = from_qp(norm2(tmp))
    end function

    pure function mf_norm2_7d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        real(qp) :: tmp_res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        tmp_res = norm2(tmp, dim=dim)
        do i6 = 1, size(tmp_res, 6)
        do i5 = 1, size(tmp_res, 5)
        do i4 = 1, size(tmp_res, 4)
        do i3 = 1, size(tmp_res, 3)
        do i2 = 1, size(tmp_res, 2)
        do i1 = 1, size(tmp_res, 1)
            res(i1, i2, i3, i4, i5, i6) = from_qp(tmp_res(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end function

    ! ================================================================
    ! findloc (pure) — rank 1-7 (mf + cx)
    ! ================================================================

    ! mf findloc rank 1
    pure function mf_findloc_1d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res(1)
        real(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_1d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_qp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 1
    pure function cx_findloc_1d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res(1)
        complex(qp) :: tmp(size(array,1))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_1d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:)
        logical, intent(in), optional :: back
        integer :: res
        complex(qp) :: tmp(size(array))
        integer :: i1
        do i1 = 1, size(array, 1)
            tmp(i1) = to_cqp(array(i1))
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 2
    pure function mf_findloc_2d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_2d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_qp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 2
    pure function cx_findloc_2d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(2)
        complex(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_2d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1))
        complex(qp) :: tmp(size(array,1), size(array,2))
        integer :: i1, i2
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2) = to_cqp(array(i1, i2))
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 3
    pure function mf_findloc_3d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(3)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_3d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_qp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 3
    pure function cx_findloc_3d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(3)
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_3d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3))
        integer :: i1, i2, i3
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3) = to_cqp(array(i1, i2, i3))
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 4
    pure function mf_findloc_4d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(4)
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_4d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        real(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_qp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 4
    pure function cx_findloc_4d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(4)
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_4d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3))
        complex(qp) :: tmp(size(array,1), size(array,2), size(array,3), size(array,4))
        integer :: i1, i2, i3, i4
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4) = to_cqp(array(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 5
    pure function mf_findloc_5d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(5)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_5d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_qp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 5
    pure function cx_findloc_5d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(5)
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_5d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5))
        integer :: i1, i2, i3, i4, i5
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5) = to_cqp(array(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 6
    pure function mf_findloc_6d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(6)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_6d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_qp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 6
    pure function cx_findloc_6d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(6)
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_6d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6))
        integer :: i1, i2, i3, i4, i5, i6
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6) = to_cqp(array(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! mf findloc rank 7
    pure function mf_findloc_7d(array, value, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        type(float64x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(7)
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), back=back)
        else
            res = findloc(tmp, to_qp(value))
        end if
    end function

    pure function mf_findloc_7d_dim(array, value, dim, mask, back) result(res)
        type(float64x2), intent(in) :: array(:, :, :, :, :, :, :)
        type(float64x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        real(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_qp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_qp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_qp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_qp(value), dim=dim)
        end if
    end function

    ! cx findloc rank 7
    pure function cx_findloc_7d(array, value, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        type(complex128x2), intent(in) :: value
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(7)
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), back=back)
        else
            res = findloc(tmp, to_cqp(value))
        end if
    end function

    pure function cx_findloc_7d_dim(array, value, dim, mask, back) result(res)
        type(complex128x2), intent(in) :: array(:, :, :, :, :, :, :)
        type(complex128x2), intent(in) :: value
        integer, intent(in) :: dim
        logical, intent(in), optional :: mask(:, :, :, :, :, :, :)
        logical, intent(in), optional :: back
        integer :: res(merge(size(array,2),size(array,1),dim<=1), &
            merge(size(array,3),size(array,2),dim<=2), &
            merge(size(array,4),size(array,3),dim<=3), &
            merge(size(array,5),size(array,4),dim<=4), &
            merge(size(array,6),size(array,5),dim<=5), &
            merge(size(array,7),size(array,6),dim<=6))
        complex(qp) :: tmp(size(array,1), &
            size(array,2), &
            size(array,3), &
            size(array,4), &
            size(array,5), &
            size(array,6), &
            size(array,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        do i7 = 1, size(array, 7)
        do i6 = 1, size(array, 6)
        do i5 = 1, size(array, 5)
        do i4 = 1, size(array, 4)
        do i3 = 1, size(array, 3)
        do i2 = 1, size(array, 2)
        do i1 = 1, size(array, 1)
            tmp(i1, i2, i3, i4, i5, i6, i7) = to_cqp(array(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        if (present(mask) .and. present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask, back=back)
        else if (present(mask)) then
            res = findloc(tmp, to_cqp(value), dim=dim, mask=mask)
        else if (present(back)) then
            res = findloc(tmp, to_cqp(value), dim=dim, back=back)
        else
            res = findloc(tmp, to_cqp(value), dim=dim)
        end if
    end function

    ! ================================================================
    ! matmul (pure) — mf + cx
    ! ================================================================

    ! mf matmul: 2D × 2D → 2D
    pure function mf_matmul_mm(a, b) result(res)
        type(float64x2), intent(in) :: a(:,:), b(:,:)
        type(float64x2) :: res(size(a,1), size(b,2))
        real(qp) :: qa(size(a,1), size(a,2)), qb(size(b,1), size(b,2))
        real(qp) :: qr(size(a,1), size(b,2))
        integer :: i1, i2
        do i2 = 1, size(a,2)
            do i1 = 1, size(a,1)
                qa(i1,i2) = to_qp(a(i1,i2))
            end do
        end do
        do i2 = 1, size(b,2)
            do i1 = 1, size(b,1)
                qb(i1,i2) = to_qp(b(i1,i2))
            end do
        end do
        qr = matmul(qa, qb)
        do i2 = 1, size(res,2)
            do i1 = 1, size(res,1)
                res(i1,i2) = from_qp(qr(i1,i2))
            end do
        end do
    end function

    ! mf matmul: 2D × 1D → 1D
    pure function mf_matmul_mv(a, b) result(res)
        type(float64x2), intent(in) :: a(:,:), b(:)
        type(float64x2) :: res(size(a,1))
        real(qp) :: qa(size(a,1), size(a,2)), qb(size(b))
        real(qp) :: qr(size(a,1))
        integer :: i1, i2
        do i2 = 1, size(a,2)
            do i1 = 1, size(a,1)
                qa(i1,i2) = to_qp(a(i1,i2))
            end do
        end do
        do i1 = 1, size(b)
            qb(i1) = to_qp(b(i1))
        end do
        qr = matmul(qa, qb)
        do i1 = 1, size(res)
            res(i1) = from_qp(qr(i1))
        end do
    end function

    ! mf matmul: 1D × 2D → 1D
    pure function mf_matmul_vm(a, b) result(res)
        type(float64x2), intent(in) :: a(:), b(:,:)
        type(float64x2) :: res(size(b,2))
        real(qp) :: qa(size(a)), qb(size(b,1), size(b,2))
        real(qp) :: qr(size(b,2))
        integer :: i1, i2
        do i1 = 1, size(a)
            qa(i1) = to_qp(a(i1))
        end do
        do i2 = 1, size(b,2)
            do i1 = 1, size(b,1)
                qb(i1,i2) = to_qp(b(i1,i2))
            end do
        end do
        qr = matmul(qa, qb)
        do i1 = 1, size(res)
            res(i1) = from_qp(qr(i1))
        end do
    end function

    ! cx matmul: 2D × 2D → 2D
    pure function cx_matmul_mm(a, b) result(res)
        type(complex128x2), intent(in) :: a(:,:), b(:,:)
        type(complex128x2) :: res(size(a,1), size(b,2))
        complex(qp) :: qa(size(a,1), size(a,2)), qb(size(b,1), size(b,2))
        complex(qp) :: qr(size(a,1), size(b,2))
        integer :: i1, i2
        do i2 = 1, size(a,2)
            do i1 = 1, size(a,1)
                qa(i1,i2) = to_cqp(a(i1,i2))
            end do
        end do
        do i2 = 1, size(b,2)
            do i1 = 1, size(b,1)
                qb(i1,i2) = to_cqp(b(i1,i2))
            end do
        end do
        qr = matmul(qa, qb)
        do i2 = 1, size(res,2)
            do i1 = 1, size(res,1)
                res(i1,i2) = from_cqp(qr(i1,i2))
            end do
        end do
    end function

    ! cx matmul: 2D × 1D → 1D
    pure function cx_matmul_mv(a, b) result(res)
        type(complex128x2), intent(in) :: a(:,:), b(:)
        type(complex128x2) :: res(size(a,1))
        complex(qp) :: qa(size(a,1), size(a,2)), qb(size(b))
        complex(qp) :: qr(size(a,1))
        integer :: i1, i2
        do i2 = 1, size(a,2)
            do i1 = 1, size(a,1)
                qa(i1,i2) = to_cqp(a(i1,i2))
            end do
        end do
        do i1 = 1, size(b)
            qb(i1) = to_cqp(b(i1))
        end do
        qr = matmul(qa, qb)
        do i1 = 1, size(res)
            res(i1) = from_cqp(qr(i1))
        end do
    end function

    ! cx matmul: 1D × 2D → 1D
    pure function cx_matmul_vm(a, b) result(res)
        type(complex128x2), intent(in) :: a(:), b(:,:)
        type(complex128x2) :: res(size(b,2))
        complex(qp) :: qa(size(a)), qb(size(b,1), size(b,2))
        complex(qp) :: qr(size(b,2))
        integer :: i1, i2
        do i1 = 1, size(a)
            qa(i1) = to_cqp(a(i1))
        end do
        do i2 = 1, size(b,2)
            do i1 = 1, size(b,1)
                qb(i1,i2) = to_cqp(b(i1,i2))
            end do
        end do
        qr = matmul(qa, qb)
        do i1 = 1, size(res)
            res(i1) = from_cqp(qr(i1))
        end do
    end function

    ! ================================================================
    ! random_number (impure subroutine)
    ! ================================================================

    subroutine mf_random_number_0d(harvest)
        type(float64x2), intent(out) :: harvest
        real(qp) :: tmp
        call random_number(tmp)
        harvest = from_qp(tmp)
    end subroutine

    subroutine mf_random_number_1d(harvest)
        type(float64x2), intent(out) :: harvest(:)
        real(qp) :: tmp(size(harvest,1))
        integer :: i1
        call random_number(tmp)
        do i1 = 1, size(tmp, 1)
            harvest(i1) = from_qp(tmp(i1))
        end do
    end subroutine

    subroutine mf_random_number_2d(harvest)
        type(float64x2), intent(out) :: harvest(:, :)
        real(qp) :: tmp(size(harvest,1), size(harvest,2))
        integer :: i1, i2
        call random_number(tmp)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2) = from_qp(tmp(i1, i2))
        end do
        end do
    end subroutine

    subroutine mf_random_number_3d(harvest)
        type(float64x2), intent(out) :: harvest(:, :, :)
        real(qp) :: tmp(size(harvest,1), size(harvest,2), size(harvest,3))
        integer :: i1, i2, i3
        call random_number(tmp)
        do i3 = 1, size(tmp, 3)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2, i3) = from_qp(tmp(i1, i2, i3))
        end do
        end do
        end do
    end subroutine

    subroutine mf_random_number_4d(harvest)
        type(float64x2), intent(out) :: harvest(:, :, :, :)
        real(qp) :: tmp(size(harvest,1), size(harvest,2), size(harvest,3), size(harvest,4))
        integer :: i1, i2, i3, i4
        call random_number(tmp)
        do i4 = 1, size(tmp, 4)
        do i3 = 1, size(tmp, 3)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2, i3, i4) = from_qp(tmp(i1, i2, i3, i4))
        end do
        end do
        end do
        end do
    end subroutine

    subroutine mf_random_number_5d(harvest)
        type(float64x2), intent(out) :: harvest(:, :, :, :, :)
        real(qp) :: tmp(size(harvest,1), &
            size(harvest,2), &
            size(harvest,3), &
            size(harvest,4), &
            size(harvest,5))
        integer :: i1, i2, i3, i4, i5
        call random_number(tmp)
        do i5 = 1, size(tmp, 5)
        do i4 = 1, size(tmp, 4)
        do i3 = 1, size(tmp, 3)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2, i3, i4, i5) = from_qp(tmp(i1, i2, i3, i4, i5))
        end do
        end do
        end do
        end do
        end do
    end subroutine

    subroutine mf_random_number_6d(harvest)
        type(float64x2), intent(out) :: harvest(:, :, :, :, :, :)
        real(qp) :: tmp(size(harvest,1), &
            size(harvest,2), &
            size(harvest,3), &
            size(harvest,4), &
            size(harvest,5), &
            size(harvest,6))
        integer :: i1, i2, i3, i4, i5, i6
        call random_number(tmp)
        do i6 = 1, size(tmp, 6)
        do i5 = 1, size(tmp, 5)
        do i4 = 1, size(tmp, 4)
        do i3 = 1, size(tmp, 3)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2, i3, i4, i5, i6) = from_qp(tmp(i1, i2, i3, i4, i5, i6))
        end do
        end do
        end do
        end do
        end do
        end do
    end subroutine

    subroutine mf_random_number_7d(harvest)
        type(float64x2), intent(out) :: harvest(:, :, :, :, :, :, :)
        real(qp) :: tmp(size(harvest,1), &
            size(harvest,2), &
            size(harvest,3), &
            size(harvest,4), &
            size(harvest,5), &
            size(harvest,6), &
            size(harvest,7))
        integer :: i1, i2, i3, i4, i5, i6, i7
        call random_number(tmp)
        do i7 = 1, size(tmp, 7)
        do i6 = 1, size(tmp, 6)
        do i5 = 1, size(tmp, 5)
        do i4 = 1, size(tmp, 4)
        do i3 = 1, size(tmp, 3)
        do i2 = 1, size(tmp, 2)
        do i1 = 1, size(tmp, 1)
            harvest(i1, i2, i3, i4, i5, i6, i7) = from_qp(tmp(i1, i2, i3, i4, i5, i6, i7))
        end do
        end do
        end do
        end do
        end do
        end do
        end do
    end subroutine

    ! ================================================================
    ! Binary arithmetic (+, -, *, /) — elemental
    ! ================================================================

    elemental function mf_add_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + to_qp(b))
    end function

    elemental function mf_add_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    elemental function mf_add_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    elemental function mf_add_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    elemental function mf_add_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + to_cqp(b))
    end function

    elemental function mf_add_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + cmplx(b, kind=qp))
    end function

    elemental function mf_add_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + cmplx(b, kind=qp))
    end function

    elemental function dp_add_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    elemental function dp_add_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    elemental function sp_add_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    elemental function sp_add_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    elemental function int_add_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    elemental function int_add_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    elemental function cx_add_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_qp(b))
    end function

    elemental function cx_add_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    elemental function cx_add_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    elemental function cx_add_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    elemental function cx_add_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_cqp(b))
    end function

    elemental function cx_add_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + cmplx(b, kind=qp))
    end function

    elemental function cx_add_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + cmplx(b, kind=qp))
    end function

    elemental function cdp_add_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_qp(b))
    end function

    elemental function cdp_add_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_cqp(b))
    end function

    elemental function csp_add_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_qp(b))
    end function

    elemental function csp_add_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_cqp(b))
    end function

    elemental function mf_sub_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - to_qp(b))
    end function

    elemental function mf_sub_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    elemental function mf_sub_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    elemental function mf_sub_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    elemental function mf_sub_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - to_cqp(b))
    end function

    elemental function mf_sub_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - cmplx(b, kind=qp))
    end function

    elemental function mf_sub_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - cmplx(b, kind=qp))
    end function

    elemental function dp_sub_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    elemental function dp_sub_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    elemental function sp_sub_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    elemental function sp_sub_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    elemental function int_sub_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    elemental function int_sub_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    elemental function cx_sub_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_qp(b))
    end function

    elemental function cx_sub_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    elemental function cx_sub_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    elemental function cx_sub_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    elemental function cx_sub_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_cqp(b))
    end function

    elemental function cx_sub_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - cmplx(b, kind=qp))
    end function

    elemental function cx_sub_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - cmplx(b, kind=qp))
    end function

    elemental function cdp_sub_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_qp(b))
    end function

    elemental function cdp_sub_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_cqp(b))
    end function

    elemental function csp_sub_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_qp(b))
    end function

    elemental function csp_sub_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_cqp(b))
    end function

    elemental function mf_mul_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * to_qp(b))
    end function

    elemental function mf_mul_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    elemental function mf_mul_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    elemental function mf_mul_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    elemental function mf_mul_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * to_cqp(b))
    end function

    elemental function mf_mul_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * cmplx(b, kind=qp))
    end function

    elemental function mf_mul_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * cmplx(b, kind=qp))
    end function

    elemental function dp_mul_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    elemental function dp_mul_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    elemental function sp_mul_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    elemental function sp_mul_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    elemental function int_mul_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    elemental function int_mul_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    elemental function cx_mul_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_qp(b))
    end function

    elemental function cx_mul_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    elemental function cx_mul_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    elemental function cx_mul_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    elemental function cx_mul_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_cqp(b))
    end function

    elemental function cx_mul_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * cmplx(b, kind=qp))
    end function

    elemental function cx_mul_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * cmplx(b, kind=qp))
    end function

    elemental function cdp_mul_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_qp(b))
    end function

    elemental function cdp_mul_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_cqp(b))
    end function

    elemental function csp_mul_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_qp(b))
    end function

    elemental function csp_mul_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_cqp(b))
    end function

    elemental function mf_div_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / to_qp(b))
    end function

    elemental function mf_div_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    elemental function mf_div_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    elemental function mf_div_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    elemental function mf_div_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / to_cqp(b))
    end function

    elemental function mf_div_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / cmplx(b, kind=qp))
    end function

    elemental function mf_div_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / cmplx(b, kind=qp))
    end function

    elemental function dp_div_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    elemental function dp_div_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    elemental function sp_div_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    elemental function sp_div_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    elemental function int_div_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    elemental function int_div_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    elemental function cx_div_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_qp(b))
    end function

    elemental function cx_div_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    elemental function cx_div_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    elemental function cx_div_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    elemental function cx_div_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_cqp(b))
    end function

    elemental function cx_div_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / cmplx(b, kind=qp))
    end function

    elemental function cx_div_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / cmplx(b, kind=qp))
    end function

    elemental function cdp_div_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_qp(b))
    end function

    elemental function cdp_div_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_cqp(b))
    end function

    elemental function csp_div_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_qp(b))
    end function

    elemental function csp_div_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_cqp(b))
    end function

    ! ================================================================
    ! Power operator (**) ��� elemental
    ! ================================================================

    elemental function mf_pow_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** to_qp(b))
    end function

    elemental function mf_pow_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** real(b, qp))
    end function

    elemental function mf_pow_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** real(b, qp))
    end function

    elemental function mf_pow_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** b)
    end function

    elemental function mf_pow_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** to_cqp(b))
    end function

    elemental function mf_pow_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** cmplx(b, kind=qp))
    end function

    elemental function mf_pow_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** cmplx(b, kind=qp))
    end function

    elemental function dp_pow_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    elemental function dp_pow_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    elemental function sp_pow_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    elemental function sp_pow_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    elemental function int_pow_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    elemental function int_pow_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    elemental function cx_pow_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** to_qp(b))
    end function

    elemental function cx_pow_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** real(b, qp))
    end function

    elemental function cx_pow_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** real(b, qp))
    end function

    elemental function cx_pow_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** b)
    end function

    elemental function cx_pow_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** to_cqp(b))
    end function

    elemental function cx_pow_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** cmplx(b, kind=qp))
    end function

    elemental function cx_pow_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** cmplx(b, kind=qp))
    end function

    elemental function cdp_pow_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_qp(b))
    end function

    elemental function cdp_pow_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_cqp(b))
    end function

    elemental function csp_pow_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_qp(b))
    end function

    elemental function csp_pow_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_cqp(b))
    end function

    ! ================================================================
    ! Equality comparison (==, /=) — elemental
    ! ================================================================

    elemental function mf_eq_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) == to_qp(b)
    end function

    elemental function mf_eq_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    elemental function mf_eq_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    elemental function mf_eq_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    elemental function mf_eq_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) == to_cqp(b)
    end function

    elemental function mf_eq_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) == cmplx(b, kind=qp)
    end function

    elemental function mf_eq_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) == cmplx(b, kind=qp)
    end function

    elemental function dp_eq_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    elemental function dp_eq_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    elemental function sp_eq_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    elemental function sp_eq_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    elemental function int_eq_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    elemental function int_eq_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    elemental function cx_eq_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) == to_qp(b)
    end function

    elemental function cx_eq_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    elemental function cx_eq_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    elemental function cx_eq_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    elemental function cx_eq_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) == to_cqp(b)
    end function

    elemental function cx_eq_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == cmplx(b, kind=qp)
    end function

    elemental function cx_eq_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == cmplx(b, kind=qp)
    end function

    elemental function cdp_eq_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_qp(b)
    end function

    elemental function cdp_eq_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_cqp(b)
    end function

    elemental function csp_eq_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_qp(b)
    end function

    elemental function csp_eq_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_cqp(b)
    end function

    elemental function mf_ne_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) /= to_qp(b)
    end function

    elemental function mf_ne_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    elemental function mf_ne_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    elemental function mf_ne_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    elemental function mf_ne_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) /= to_cqp(b)
    end function

    elemental function mf_ne_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= cmplx(b, kind=qp)
    end function

    elemental function mf_ne_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= cmplx(b, kind=qp)
    end function

    elemental function dp_ne_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    elemental function dp_ne_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    elemental function sp_ne_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    elemental function sp_ne_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    elemental function int_ne_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    elemental function int_ne_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    elemental function cx_ne_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= to_qp(b)
    end function

    elemental function cx_ne_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    elemental function cx_ne_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    elemental function cx_ne_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    elemental function cx_ne_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= to_cqp(b)
    end function

    elemental function cx_ne_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= cmplx(b, kind=qp)
    end function

    elemental function cx_ne_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= cmplx(b, kind=qp)
    end function

    elemental function cdp_ne_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_qp(b)
    end function

    elemental function cdp_ne_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_cqp(b)
    end function

    elemental function csp_ne_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_qp(b)
    end function

    elemental function csp_ne_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_cqp(b)
    end function

    ! ================================================================
    ! Ordered comparison (<, >, <=, >=) — elemental
    ! ================================================================

    elemental function mf_lt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) < to_qp(b)
    end function

    elemental function mf_lt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    elemental function mf_lt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    elemental function mf_lt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    elemental function dp_lt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    elemental function sp_lt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    elemental function int_lt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    elemental function mf_gt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) > to_qp(b)
    end function

    elemental function mf_gt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    elemental function mf_gt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    elemental function mf_gt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    elemental function dp_gt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    elemental function sp_gt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    elemental function int_gt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    elemental function mf_le_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) <= to_qp(b)
    end function

    elemental function mf_le_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    elemental function mf_le_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    elemental function mf_le_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    elemental function dp_le_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    elemental function sp_le_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    elemental function int_le_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    elemental function mf_ge_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) >= to_qp(b)
    end function

    elemental function mf_ge_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    elemental function mf_ge_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    elemental function mf_ge_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    elemental function dp_ge_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    elemental function sp_ge_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    elemental function int_ge_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

end module multifloats
