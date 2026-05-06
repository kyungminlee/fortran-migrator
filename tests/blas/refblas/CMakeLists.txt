# refblas_quad — quad-precision Netlib BLAS reference.
#
# Built from the vendored Netlib BLAS sources (copied to
# ${staging}/_refblas_src by `migrator stage`) compiled with gfortran's
# -freal-8-real-16, which promotes REAL(KIND=8) and DOUBLE PRECISION
# entities to REAL(KIND=16) without renaming routines. The resulting
# library exposes the standard D/Z routine names (DGEMM, DDOT, ZAXPY,
# ...) operating at quad precision — used as the universal reference
# for the differential precision tests.
#
# Selection (from parent tests/blas/CMakeLists.txt):
#   -DBLAS_REF=netlib  (default) → build from vendored sources
#   -DBLAS_REF=system            → link system -lblas (double, capped)

if(BLAS_REF STREQUAL "system")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

# Vendored-source mode requires gfortran for the -freal-8-real-16 flag.
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "tests/blas/refblas: BLAS_REF=netlib needs gfortran "
                   "for -freal-8-real-16 promotion. Falling back to "
                   "system -lblas.")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

set(_refsrc_dir ${CMAKE_SOURCE_DIR}/_refblas_src)
if(NOT IS_DIRECTORY ${_refsrc_dir})
    message(STATUS "tests/blas/refblas: ${_refsrc_dir} not found "
                   "(expected to be copied by `migrator stage`); "
                   "falling back to system -lblas.")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

file(GLOB _refsrcs_f   CONFIGURE_DEPENDS ${_refsrc_dir}/*.f)
file(GLOB _refsrcs_f90 CONFIGURE_DEPENDS ${_refsrc_dir}/*.f90)
set(_refsrcs ${_refsrcs_f} ${_refsrcs_f90})
# Skip mixed-precision routines (single+double): not promotable cleanly.
list(FILTER _refsrcs EXCLUDE REGEX "/(dsdot|sdsdot)\\.f9?0?$")

add_library(refblas_quad STATIC ${_refsrcs})
# -freal-8-real-16 promotes REAL(KIND=8) and DOUBLE PRECISION (= KIND=8
# by default) entities to REAL(KIND=16), as well as their COMPLEX
# counterparts. D/Z literals (1.0D0, (1.0D0,2.0D0)) get promoted too.
# S/C single-precision routines are unaffected (KIND=4 stays KIND=4).
# -std=legacy is gfortran-only; -w silences the avalanche of "obsolete"
# warnings from the upstream f77 sources.
target_compile_options(refblas_quad PRIVATE
    -freal-8-real-16
    -w
    -std=legacy
)
# .f90 files are free-form by extension; the rest are fixed-form.
# No per-file form override needed.
