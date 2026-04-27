# Differential precision tests for migrated BLACS.
#
# BLACS is a thin Fortran-callable layer over MPI providing process
# grid management plus per-type point-to-point send/recv, broadcast,
# and combine (sum/amx/amn) primitives. Communication is exact; we
# verify byte-equivalent payloads at REAL(KIND=16) using a 0.0
# tolerance (no roundoff).
#
# Each test program runs on a 2x2 BLACS grid (default mpiexec -n 4),
# prepares known per-rank data, calls the migrated entry point through
# a per-target wrapper, and asserts exact equality.

if(NOT HAVE_REAL16)
    message(STATUS "tests/blacs: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/blacs: skipping — MPI Fortran not found")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}blacs)
    message(STATUS "tests/blacs: skipping — ${LIB_PREFIX}blacs not built")
    return()
endif()

# ── Common helper modules ───────────────────────────────────────────
# prec_kinds / compare / test_data / pblas_grid are byte-for-byte
# identical to the copies in tests/pbblas; shared module names produce
# the same .mod in the project-wide fmod directory, no collision.
add_library(blacs_test_common STATIC
    common/prec_kinds.f90
    common/compare.f90
    common/blacs_prec_report.f90
    common/test_data.f90
    common/pblas_grid.f90
)
fortran_module_layout(blacs_test_common)
target_link_libraries(blacs_test_common PUBLIC
    ${LIB_PREFIX}blacs
    MPI::MPI_Fortran
)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/blacs: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_blacs_wrapper_src
    ${_wrapper_dir}/target_blacs.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(blacs_test_target STATIC ${_blacs_wrapper_src})
fortran_module_layout(blacs_test_target)
target_link_libraries(blacs_test_target PUBLIC
    blacs_test_common
    target_conv
    ${LIB_PREFIX}blacs
)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(blacs_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# Default 4 ranks → 2x2 grid. Override at configure time:
# -DBLACS_TEST_NPROC=<count>.
set(BLACS_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/blacs (default 4 → 2x2 grid)")

function(add_blacs_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    set(_blacs_cycle
        blacs_test_target blacs_test_common
        ${LIB_PREFIX}blacs)
    foreach(_common blas_common blacs_common)
        if(TARGET ${_common})
            list(APPEND _blacs_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_blacs_cycle}>")
    add_test(NAME blacs_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${BLACS_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

file(GLOB _tests CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.f90)
foreach(_t ${_tests})
    add_blacs_precision_test(${_t})
endforeach()
