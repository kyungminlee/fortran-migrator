# Differential precision tests for migrated ScaLAPACK.
#
# Mirrors tests/pblas/CMakeLists.txt: each test program generates a
# global problem at REAL(KIND=16), distributes block-cyclicly over a
# BLACS grid, calls the migrated q/x-prefix ScaLAPACK routine (via the
# per-target wrapper), gathers the distributed result to rank 0, and
# compares against serial reflapack_quad (or refblas_quad for pure-BLAS
# auxiliaries) applied to the full global matrix.
#
# ScaLAPACK depends on BLAS + LAPACK + BLACS + PBLAS, so all four
# migrated libraries must be present alongside the reference quad
# libraries for this suite to compile.

if(NOT HAVE_REAL16)
    message(STATUS "tests/scalapack: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/scalapack: skipping — MPI Fortran not found")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/scalapack: skipping — refblas_quad not available "
                   "(tests/blas must be staged)")
    return()
endif()

if(NOT TARGET reflapack_quad)
    message(STATUS "tests/scalapack: skipping — reflapack_quad not available "
                   "(tests/lapack must be staged)")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}scalapack)
    message(STATUS "tests/scalapack: skipping — ${LIB_PREFIX}scalapack not built")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}pblas)
    message(STATUS "tests/scalapack: skipping — ${LIB_PREFIX}pblas not built")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}blacs)
    message(STATUS "tests/scalapack: skipping — ${LIB_PREFIX}blacs not built")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}lapack)
    message(STATUS "tests/scalapack: skipping — ${LIB_PREFIX}lapack not built")
    return()
endif()

# ── BLACS C-API shim ───────────────────────────────────────────────
# Same rationale as tests/pblas: PBLAS's C code calls Cblacs_* /
# C{q,e,dd}gsum2d counterparts that upstream BLACS produces by
# compiling every source twice (once with -DCallFromC). The migrated
# ${LIB_PREFIX}blacs only carries the Fortran-callable entry points,
# so we replicate the -DCallFromC pass locally. Target name differs
# from tests/pblas's to avoid collision when both suites stage.
set(_blacs_src_dir ${CMAKE_SOURCE_DIR}/blacs/src)
if(IS_DIRECTORY ${_blacs_src_dir})
    file(GLOB _blacs_c_api_supp  CONFIGURE_DEPENDS
        ${_blacs_src_dir}/blacs_*.c
        ${_blacs_src_dir}/sys2blacs_.c
        ${_blacs_src_dir}/blacs2sys_.c
        ${_blacs_src_dir}/free_handle_.c)
    file(GLOB _blacs_c_api_comm  CONFIGURE_DEPENDS
        ${_blacs_src_dir}/*gesd2d_.c ${_blacs_src_dir}/*gerv2d_.c
        ${_blacs_src_dir}/*gebs2d_.c ${_blacs_src_dir}/*gebr2d_.c
        ${_blacs_src_dir}/*gsum2d_.c ${_blacs_src_dir}/*gamx2d_.c
        ${_blacs_src_dir}/*gamn2d_.c
        ${_blacs_src_dir}/*trsd2d_.c ${_blacs_src_dir}/*trrv2d_.c
        ${_blacs_src_dir}/*trbs2d_.c ${_blacs_src_dir}/*trbr2d_.c)
    add_library(scalapack_test_blacs_c_api STATIC
        ${_blacs_c_api_supp} ${_blacs_c_api_comm})
    target_compile_definitions(scalapack_test_blacs_c_api PRIVATE CallFromC)
    target_include_directories(scalapack_test_blacs_c_api PRIVATE ${_blacs_src_dir})
    if(MPI_C_FOUND)
        target_link_libraries(scalapack_test_blacs_c_api PRIVATE MPI::MPI_C)
    endif()
    if(C_AS_CXX)
        set_source_files_properties(${_blacs_c_api_supp} ${_blacs_c_api_comm}
            PROPERTIES LANGUAGE CXX)
        target_compile_features(scalapack_test_blacs_c_api PRIVATE cxx_std_17)
        target_compile_options(scalapack_test_blacs_c_api PRIVATE
            -fpermissive -Wno-write-strings)
    endif()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# tests/scalapack/common reuses the same module names as
# tests/pblas/common (prec_kinds, compare, pblas_prec_report,
# test_data, ref_quad_blas, pblas_grid, pblas_distrib). Placing the
# ScaLAPACK variants' .mod files in their own subdirectory under fmod
# prevents collisions when both test trees build together.
set(_sp_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/scalapack_tests")
file(MAKE_DIRECTORY ${_sp_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(scalapack_test_common STATIC
    common/prec_kinds.f90
    common/compare.f90
    common/prec_report.f90
    common/test_data.f90
    common/ref_quad_blas.f90
    common/ref_quad_lapack.f90
    common/pblas_grid.f90
    common/pblas_distrib.f90
)
set_target_properties(scalapack_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_sp_moddir}"
)
target_include_directories(scalapack_test_common PUBLIC
    $<BUILD_INTERFACE:${_sp_moddir}>
)
target_link_libraries(scalapack_test_common PUBLIC
    refblas_quad
    reflapack_quad
    ${LIB_PREFIX}blacs
    MPI::MPI_Fortran
)
if(TARGET scalapack_test_blacs_c_api)
    target_link_libraries(scalapack_test_common PUBLIC scalapack_test_blacs_c_api)
endif()

# pilaenv.f vendored for PBLAS (see tests/pblas/CMakeLists.txt).
add_library(scalapack_test_pilaenv STATIC common/pilaenv.f)
set_target_properties(scalapack_test_pilaenv PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_options(scalapack_test_pilaenv PRIVATE -w -std=legacy)
target_link_libraries(scalapack_test_common PUBLIC scalapack_test_pilaenv)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/scalapack: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

set(_wrapper_sources ${_wrapper_dir}/target_scalapack.f90)
if(EXISTS ${_wrapper_dir}/ptzblas_stubs.f90)
    list(APPEND _wrapper_sources ${_wrapper_dir}/ptzblas_stubs.f90)
endif()
add_library(scalapack_test_target STATIC ${_wrapper_sources})
set_target_properties(scalapack_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_sp_moddir}"
)
target_include_directories(scalapack_test_target PUBLIC
    $<BUILD_INTERFACE:${_sp_moddir}>
)
target_link_libraries(scalapack_test_target PUBLIC
    scalapack_test_common
    ${LIB_PREFIX}scalapack
    ${LIB_PREFIX}pblas
    ${LIB_PREFIX}lapack
    ${LIB_PREFIX}blas
)
foreach(_dep ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas)
    if(TARGET ${_dep})
        target_link_libraries(scalapack_test_target PUBLIC ${_dep})
    endif()
endforeach()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(scalapack_test_target PUBLIC multifloatsf)
endif()

# ── LAMCH shim (linked AFTER ${LIB_PREFIX}lapack) ───────────────────
# Same workaround as tests/lapack: the migrator does not stage
# external/lapack-3.12.1/INSTALL, so {q,e,dd}lamch / roundup_lwork
# are unresolved. ScaLAPACK routines call through to LAPACK's lamch,
# so without the shim the test executables fail to link.
if(EXISTS ${_wrapper_dir}/lamch_shim.f90)
    add_library(scalapack_test_shim STATIC ${_wrapper_dir}/lamch_shim.f90)
    set_target_properties(scalapack_test_shim PROPERTIES
        Fortran_MODULE_DIRECTORY "${_sp_moddir}"
    )
    if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
        target_link_libraries(scalapack_test_shim PUBLIC multifloatsf)
    endif()
    target_link_libraries(scalapack_test_target PUBLIC scalapack_test_shim)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

set(SCALAPACK_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/scalapack (default 4 → 2x2 grid)")

function(add_scalapack_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # -fno-lto: same LTO-confusion issue documented in tests/pblas.
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # ScaLAPACK ↔ LAPACK ↔ PBLAS ↔ PTZBLAS ↔ BLACS ↔ BLAS form a
    # cyclic archive graph. Wrap the whole cycle in LINK_GROUP(RESCAN)
    # so ld re-scans member archives until all cross-references resolve.
    set(_sp_cycle
        scalapack_test_target scalapack_test_common
        ${LIB_PREFIX}scalapack ${LIB_PREFIX}pblas
        ${LIB_PREFIX}ptzblas  ${LIB_PREFIX}pbblas
        ${LIB_PREFIX}lapack   ${LIB_PREFIX}blas
        ${LIB_PREFIX}blacs
        reflapack_quad refblas_quad)
    foreach(_common blas_common blacs_common ptzblas_common pblas_common lapack_common scalapack_common)
        if(TARGET ${_common})
            list(APPEND _sp_cycle ${_common})
        endif()
    endforeach()
    if(TARGET scalapack_test_shim)
        list(APPEND _sp_cycle scalapack_test_shim)
    endif()
    if(TARGET scalapack_test_blacs_c_api)
        list(APPEND _sp_cycle scalapack_test_blacs_c_api)
    endif()
    if(TARGET scalapack_test_pilaenv)
        list(APPEND _sp_cycle scalapack_test_pilaenv)
    endif()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_sp_cycle}>")
    add_test(NAME scalapack_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SCALAPACK_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_group linear_solve factorization eigenvalue svd auxiliary)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_group}/test_*.f90)
    foreach(_t ${_tests})
        add_scalapack_precision_test(${_t})
    endforeach()
endforeach()
