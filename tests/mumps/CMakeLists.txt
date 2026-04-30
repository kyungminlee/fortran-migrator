# Differential precision tests for migrated MUMPS.
#
# Mirrors tests/lapack/CMakeLists.txt with two extras:
#   1. C-side bridge — upstream mumps_c.c is compiled twice (Q + X arith)
#      against tests/mumps/c/include/'s quad-precision type shadow, plus
#      type-agnostic upstream helpers (mumps_common.c, mumps_addr.c, IO
#      / save / thread / metis / pord / scotch). See tests/mumps/TODO.md
#      B2 for the override design.
#   2. Every test runs via ``mpiexec -n 1`` — the migrated qmumps Fortran
#      calls MPI primitives (mpi_init_, mpi_bcast_, ...) regardless of
#      NPROCS. Real MPI is required at runtime; libmpiseq exists as an
#      alternative (B3 in TODO.md) but its full sequential-link path is
#      not yet wired in.
#
# Outputs: per-test JSON reports under
# ``${CMAKE_BINARY_DIR}/precision_reports/<test_name>.<target>.json``.

# ── Gate checks ──────────────────────────────────────────────────────
if(NOT TARGET ${LIB_PREFIX}mumps)
    message(STATUS "tests/mumps: skipping — ${LIB_PREFIX}mumps not built")
    return()
endif()
if(NOT HAVE_REAL16)
    message(STATUS "tests/mumps: skipping — HAVE_REAL16 required")
    return()
endif()
if(NOT TARGET reflapack_quad)
    message(STATUS "tests/mumps: skipping — reflapack_quad not available "
                   "(tests/lapack must be enabled)")
    return()
endif()
if(NOT MPI_Fortran_FOUND OR NOT MPI_C_FOUND)
    message(STATUS "tests/mumps: skipping — MPI Fortran+C not found")
    return()
endif()

set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/mumps: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

set(_MUMPS_C_SRC ${PROJECT_SOURCE_DIR}/_mumps_upstream_src)
set(_MUMPS_C_INC ${PROJECT_SOURCE_DIR}/_mumps_upstream_include)
if(NOT IS_DIRECTORY ${_MUMPS_C_SRC} OR NOT IS_DIRECTORY ${_MUMPS_C_INC})
    message(STATUS "tests/mumps: skipping — upstream MUMPS C sources not "
                   "staged at ${_MUMPS_C_SRC} / ${_MUMPS_C_INC}")
    return()
endif()

set(_BRIDGE_INC ${CMAKE_CURRENT_SOURCE_DIR}/c/include)

# ── Output directory for JSON reports ────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# ── Dedicated module directory ───────────────────────────────────────
# Same prec_kinds / compare / prec_report / test_data module names as
# tests/blas and tests/lapack, so isolate the .mod outputs to avoid
# collision when multiple test trees build together.
set(_mp_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/mumps_tests")
file(MAKE_DIRECTORY ${_mp_moddir})

# ── C-side bridge ─────────────────────────────────────────────────────
# Per-arithmetic dispatch — compile mumps_c.c twice (Q and X) as
# OBJECT libraries so they can be folded into the single combined
# archive below. Override Fortran symbols so qmumps_c() lands in the
# bridge but the dmumps_assign_*/dmumps_set_tmp_ptr_c_ helpers stay
# under their D-prefix (because the migrated qmumps_f77_ Fortran
# expects them that way — same convention for Z/X).
add_library(qmumps_c_bridge OBJECT ${_MUMPS_C_SRC}/mumps_c.c)
target_include_directories(qmumps_c_bridge PRIVATE
    ${_BRIDGE_INC} ${_MUMPS_C_INC} ${_MUMPS_C_SRC}
)
target_compile_definitions(qmumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_d Add_
    dmumps_f77_=qmumps_f77_
    dmumps_set_tmp_ptr_=qmumps_set_tmp_ptr_
    dmumps_c=qmumps_c
    DMUMPS_STRUC_C=QMUMPS_STRUC_C
)
# Force-include our quad-precision mumps_c_types.h shadow BEFORE
# mumps_c.c gets a chance to pull in upstream's via mumps_common.h's
# `#include "mumps_c_types.h"` (which would otherwise resolve to
# upstream's via same-directory search).
target_compile_options(qmumps_c_bridge PRIVATE
    -include mumps_c_types.h)

add_library(xmumps_c_bridge OBJECT ${_MUMPS_C_SRC}/mumps_c.c)
target_include_directories(xmumps_c_bridge PRIVATE
    ${_BRIDGE_INC} ${_MUMPS_C_INC} ${_MUMPS_C_SRC}
)
target_compile_definitions(xmumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_z Add_
    zmumps_f77_=xmumps_f77_
    zmumps_set_tmp_ptr_=xmumps_set_tmp_ptr_
    zmumps_c=xmumps_c
    ZMUMPS_STRUC_C=XMUMPS_STRUC_C
)
target_compile_options(xmumps_c_bridge PRIVATE
    -include mumps_c_types.h)

# mumps_addr.c stays at upstream double precision because the migrated
# Fortran callers (qtools.F's QMUMPS_UPDATE_PROGRESS, xtools.F's
# XMUMPS_UPDATE_PROGRESS) pass DOUBLE PRECISION OPELI — the keep-kind
# manifest preserved DP for those callers. Compiling mumps_addr.c with
# DMUMPS_REAL=__float128 overrides would make MUMPS_SETDVAL_ADDR_C
# write 16 bytes into the 8-byte caller buffer (silent corruption →
# SIGSEGV inside qmumps_update_progress_).
add_library(mumps_c_addr_obj OBJECT ${_MUMPS_C_SRC}/mumps_addr.c)
# Include path order matters: list _MUMPS_C_INC FIRST so upstream's
# `double`-typed mumps_c_types.h gets resolved by mumps_addr.c's
# `#include "mumps_c_types.h"` chain. _BRIDGE_INC is still on the
# search path so the upstream's `#include "mumps_int_def.h"` falls
# through to our generated stand-in (upstream ships only a .in
# template).
target_include_directories(mumps_c_addr_obj PRIVATE
    ${_MUMPS_C_INC} ${_MUMPS_C_SRC} ${_BRIDGE_INC})
target_compile_definitions(mumps_c_addr_obj PRIVATE Add_)

# Combined C-side static archive. Bundles the per-arithmetic bridge
# objects with the type-agnostic upstream runtime (mumps_common.c,
# IO, save/restore, thread, numa, flytes, config, metis, pord, scotch)
# into ONE archive. Doing this rather than two archives sidesteps a
# static-link ordering problem: the migrated qmumps Fortran (built
# upstream of tests/) calls into both the bridge's
# dmumps_set_tmp_ptr_c_ and the runtime's mumps_size_c_ /
# mumps_low_level_*, and a separate runtime archive listed before
# libqmumps would have its objects skipped at scan time and not
# re-pulled when libqmumps's later objects request them.
add_library(mumps_c_bridge STATIC
    $<TARGET_OBJECTS:qmumps_c_bridge>
    $<TARGET_OBJECTS:xmumps_c_bridge>
    $<TARGET_OBJECTS:mumps_c_addr_obj>
    ${CMAKE_CURRENT_SOURCE_DIR}/c/c_parity_helpers.c
    ${_MUMPS_C_SRC}/mumps_common.c
    ${_MUMPS_C_SRC}/mumps_io.c
    ${_MUMPS_C_SRC}/mumps_io_basic.c
    ${_MUMPS_C_SRC}/mumps_io_err.c
    ${_MUMPS_C_SRC}/mumps_io_thread.c
    ${_MUMPS_C_SRC}/mumps_save_restore_C.c
    ${_MUMPS_C_SRC}/mumps_register_thread.c
    ${_MUMPS_C_SRC}/mumps_thread.c
    ${_MUMPS_C_SRC}/mumps_thread_affinity.c
    ${_MUMPS_C_SRC}/mumps_numa.c
    ${_MUMPS_C_SRC}/mumps_flytes.c
    ${_MUMPS_C_SRC}/mumps_config_file_C.c
    ${_MUMPS_C_SRC}/mumps_metis.c
    ${_MUMPS_C_SRC}/mumps_metis64.c
    ${_MUMPS_C_SRC}/mumps_metis_int.c
    ${_MUMPS_C_SRC}/mumps_pord.c
    ${_MUMPS_C_SRC}/mumps_scotch.c
    ${_MUMPS_C_SRC}/mumps_scotch64.c
    ${_MUMPS_C_SRC}/mumps_scotch_int.c
)
target_include_directories(mumps_c_bridge PUBLIC
    $<BUILD_INTERFACE:${_BRIDGE_INC}>      # FIRST — quad-precision type shadow
    $<BUILD_INTERFACE:${_MUMPS_C_INC}>     # then upstream include/
    $<BUILD_INTERFACE:${_MUMPS_C_SRC}>     # then upstream src/ for mumps_common.h etc.
)
target_compile_definitions(mumps_c_bridge PRIVATE Add_)
# Override applies to everything in mumps_c_bridge EXCEPT the
# already-built mumps_c_addr_obj (which we added via $<TARGET_OBJECTS>
# so its compile flags are frozen). The non-_addr files are either
# integer-only (mumps_common.c, IO/save/thread/etc. helpers — the
# override is harmless) or they're the per-arithmetic bridge sources
# that NEED __float128.
target_compile_options(mumps_c_bridge PRIVATE
    -include mumps_c_types.h)
target_link_libraries(mumps_c_bridge PUBLIC MPI::MPI_C)

# ── Common helper modules ────────────────────────────────────────────
add_library(mumps_test_common STATIC
    common/prec_kinds.f90
    common/compare.f90
    common/prec_report.f90
    common/test_data_mumps.f90
    common/ref_quad_lapack_solve.f90
)
set_target_properties(mumps_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_mp_moddir}"
)
target_include_directories(mumps_test_common PUBLIC
    $<BUILD_INTERFACE:${_mp_moddir}>
)
target_link_libraries(mumps_test_common PUBLIC reflapack_quad)

# ── Per-target wrapper module (kind16-only today) ────────────────────
add_fypp_source(_mumps_wrapper_src
    ${_wrapper_dir}/target_mumps.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(mumps_test_target STATIC ${_mumps_wrapper_src})
set_target_properties(mumps_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_mp_moddir}"
)
target_include_directories(mumps_test_target PUBLIC
    $<BUILD_INTERFACE:${_mp_moddir}>
)
target_link_libraries(mumps_test_target PUBLIC mumps_test_common)
# qmumps + scalapack + lapack + blas + mumps_c_bridge are NOT attached
# as PUBLIC dependencies because cmake's PUBLIC propagation would place
# them OUTSIDE the per-test LINK_GROUP wrapper. They're listed
# explicitly inside add_mumps_*_test()'s LINK_GROUP so every archive
# carrying cross-references (qmumps Fortran ↔ bridge's
# dmumps_set_tmp_ptr_c_ / mumps_size_c_ / mumps_low_level_*; bridge ↔
# qmumps_f77_) gets re-scanned until everything resolves. The wrapper
# Fortran source still needs qmumps's .mod files at compile time, so
# pull just the include dir without taking a link dependency.
get_target_property(_qmumps_moddir ${LIB_PREFIX}mumps Fortran_MODULE_DIRECTORY)
if(_qmumps_moddir)
    target_include_directories(mumps_test_target PUBLIC
        $<BUILD_INTERFACE:${_qmumps_moddir}>)
endif()
# Build-order dependency: qmumps_struc_def.mod / xmumps_struc_def.mod
# need to exist before target_mumps.f90 compiles. We didn't set
# qmumps as a link dependency (PUBLIC propagation would break the
# LINK_GROUP wrapping), so make it explicit here.
add_dependencies(mumps_test_target ${LIB_PREFIX}mumps)

# ── Test executables ──────────────────────────────────────────────────
function(_mumps_register_test name)
    add_test(NAME mumps_${name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}>
                     ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
    set_tests_properties(mumps_${name} PROPERTIES RUN_SERIAL TRUE)
endfunction()

# All static archives that carry cross-references between the migrated
# qmumps Fortran and the C bridge — collected here so each test wraps
# them in a single LINK_GROUP RESCAN.
set(_mumps_link_group_archives
    mumps_test_target
    mumps_c_bridge
    ${LIB_PREFIX}mumps
    ${LIB_PREFIX}scalapack
    ${LIB_PREFIX}scalapack_c
    ${LIB_PREFIX}pblas
    ${LIB_PREFIX}pbblas
    ${LIB_PREFIX}ptzblas
    ${LIB_PREFIX}blacs
    ${LIB_PREFIX}lapack
    ${LIB_PREFIX}blas
)

function(add_mumps_fortran_test test_src)
    get_filename_component(name ${test_src} NAME_WE)
    add_executable(${name} ${test_src})
    target_link_libraries(${name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_mumps_link_group_archives}>"
        MPI::MPI_Fortran
    )
    # Use the native objects in the LTO-fat archives (built with
    # -flto=auto -ffat-lto-objects upstream). Cross-archive LTO
    # partitioning produces spurious "undefined reference to
    # xmumps_mtransi_" / "pxtrsm_" etc. when linking via the LTO
    # plugin, even though those symbols ARE in the archives. Native
    # link sidesteps this.
    target_link_options(${name} PRIVATE -fno-lto)
    _mumps_register_test(${name})
endfunction()

function(add_mumps_c_test test_src)
    get_filename_component(name ${test_src} NAME_WE)
    add_executable(${name} ${test_src})
    # Same Fortran-linker / no-LTO setup as the Fortran tests — the
    # bridge archive contains gfortran-LTO bitcode pulled in via
    # libqmumps's LTO partitioning, and ``cc`` doesn't have the
    # gfortran LTO plugin loaded.
    set_target_properties(${name} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_mumps_link_group_archives}>"
        MPI::MPI_Fortran
        MPI::MPI_C
    )
    target_link_options(${name} PRIVATE -fno-lto)
    target_include_directories(${name} PRIVATE
        ${_BRIDGE_INC} ${_MUMPS_C_INC}
    )
    _mumps_register_test(${name})
endfunction()

file(GLOB _ftests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_*.f90)
foreach(_t ${_ftests})
    add_mumps_fortran_test(${_t})
endforeach()

file(GLOB _ctests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/c/test_*.c)
foreach(_t ${_ctests})
    add_mumps_c_test(${_t})
endforeach()
