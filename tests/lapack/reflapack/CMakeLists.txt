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
file(GLOB _reflap_F   CONFIGURE_DEPENDS ${_reflap_dir}/*.F)
file(GLOB _reflap_f90 CONFIGURE_DEPENDS ${_reflap_dir}/*.f90)
file(GLOB _reflap_F90 CONFIGURE_DEPENDS ${_reflap_dir}/*.F90)
set(_reflap_srcs ${_reflap_f} ${_reflap_F} ${_reflap_f90} ${_reflap_F90})

# Skip mixed-precision / cross-precision routines (single+double or
# complex+real32 bridges) that don't promote cleanly under
# -freal-8-real-16. Match either case of the .f / .F suffix in case
# upstream Netlib ever renames a bridge file to capital F.
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dsgesv|zcgesv|dsposv|zcposv|dsgels|zcgels)\\.[fF]9?0?$")
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dlag2s|slag2d|zlag2c|clag2z|dlat2s|zlat2c)\\.[fF]9?0?$")

# Skip multifloats-specific auxiliary modules (la_*_mf / la_*_ep are
# for targets that don't use the standard -freal-8-real-16 promotion).
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_ep")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_ep")

# Skip iparam2stage.F. It is a precision-independent INTEGER function
# (no FP argument promotion), and the staged lapack/src/iparam2stage.F
# already ships a patched copy via recipes/lapack/mf_overrides/ that
# recognizes the migrator precision prefixes (Q / X / E / Y / T / V).
# If we let reflapack_quad keep its own iparam2stage_ symbol the linker
# resolves it from this archive first, which masks the patched copy in
# liblapack_common — and qsyev_2stage / vheev_2stage / etc. then get
# IPARAM2STAGE = -1 on the migrated path, undersized LWORK, and heap
# corruption inside libc free. Excluding the upstream copy here lets
# the linker pull the patched version from liblapack_common for both
# the reference (D-prefix) and migrated (Q/E/T/X/Y/V-prefix) call
# sites — they're ABI-compatible because iparam2stage takes only
# INTEGER and CHARACTER arguments.
list(FILTER _reflap_srcs EXCLUDE REGEX "/iparam2stage\\.F$")

add_library(reflapack_quad STATIC ${_reflap_srcs})
target_compile_options(reflapack_quad PRIVATE
    -freal-8-real-16
    -w
    -std=legacy
)
# LAPACK internally calls BLAS. Link the quad-promoted refblas_quad so
# every routine in the reference chain operates at KIND=16.
target_link_libraries(reflapack_quad PUBLIC refblas_quad)
