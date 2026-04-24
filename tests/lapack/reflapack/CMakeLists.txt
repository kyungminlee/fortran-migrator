# reflapack_quad — quad-precision Netlib LAPACK reference.
#
# Built from the vendored Netlib LAPACK sources (copied to
# ${staging}/_reflapack_src by `pyengine stage`) compiled with gfortran's
# -freal-8-real-16, which promotes REAL(KIND=8) and DOUBLE PRECISION
# entities to REAL(KIND=16) without renaming routines. The resulting
# library exposes the standard D/Z routine names (DGESV, DGETRF, ZHEEV,
# ...) operating at quad precision — used as the universal reference
# for the differential precision tests.
#
# Selection (from parent tests/lapack/CMakeLists.txt):
#   -DLAPACK_REF=netlib  (default) → build from vendored sources
#   -DLAPACK_REF=system            → link system -llapack (double, capped)

if(LAPACK_REF STREQUAL "system")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

# Vendored-source mode requires gfortran for the -freal-8-real-16 flag.
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "tests/lapack/reflapack: LAPACK_REF=netlib needs gfortran "
                   "for -freal-8-real-16 promotion. Falling back to "
                   "system -llapack.")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

set(_reflap_dir ${CMAKE_SOURCE_DIR}/_reflapack_src)
if(NOT IS_DIRECTORY ${_reflap_dir})
    message(STATUS "tests/lapack/reflapack: ${_reflap_dir} not found "
                   "(expected to be copied by `pyengine stage`); "
                   "falling back to system -llapack.")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

file(GLOB _reflap_f   CONFIGURE_DEPENDS ${_reflap_dir}/*.f)
file(GLOB _reflap_f90 CONFIGURE_DEPENDS ${_reflap_dir}/*.f90)
file(GLOB _reflap_FF  CONFIGURE_DEPENDS ${_reflap_dir}/*.F90)
set(_reflap_srcs ${_reflap_f} ${_reflap_f90} ${_reflap_FF})

# Skip mixed-precision / cross-precision routines (single+double or
# complex+real32 bridges) that don't promote cleanly under
# -freal-8-real-16.
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dsgesv|zcgesv|dsposv|zcposv|dsgels|zcgels)\\.f9?0?$")
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dlag2s|slag2d|zlag2c|clag2z|dlat2s|zlat2c)\\.f9?0?$")

# Skip multifloats-specific auxiliary modules (la_*_mf / la_*_ep are
# for targets that don't use the standard -freal-8-real-16 promotion).
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_ep")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_ep")

add_library(reflapack_quad STATIC ${_reflap_srcs})
target_compile_options(reflapack_quad PRIVATE
    -freal-8-real-16
    -w
    -std=legacy
)
# LAPACK internally calls BLAS. Link the quad-promoted refblas_quad so
# every routine in the reference chain operates at KIND=16.
target_link_libraries(reflapack_quad PUBLIC refblas_quad)
