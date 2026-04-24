# Differential precision tests for migrated PBLAS.
#
# Each test program generates a global problem at REAL(KIND=16),
# distributes block-cyclicly over a BLACS grid, calls the migrated
# q/e/dd-prefix PBLAS routine (via the per-target wrapper), gathers
# the distributed result to rank 0, and compares against serial
# refblas_quad applied to the full global matrix.
#
# The serial BLAS reference is supplied by the refblas_quad target
# built under tests/blas/refblas — i.e. tests/blas must also be
# staged for the PBLAS tests to compile.

if(NOT HAVE_REAL16)
    message(STATUS "tests/pblas: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/pblas: skipping — MPI Fortran not found "
                   "(required for BLACS grid init and rank-0 gather)")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/pblas: skipping — refblas_quad target not built "
                   "(tests/blas must be staged alongside tests/pblas)")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}pblas)
    message(STATUS "tests/pblas: skipping — ${LIB_PREFIX}pblas not built")
    return()
endif()

if(NOT TARGET ${LIB_PREFIX}blacs)
    message(STATUS "tests/pblas: skipping — ${LIB_PREFIX}blacs not built")
    return()
endif()

# ── BLACS C-API shim (Cblacs_* / C{q,e,dd}gsum2d / etc.) ────────────
# The migrated ${LIB_PREFIX}blacs compiles BLACS sources once and gets
# only the Fortran-callable `blacs_*_` / `{q,e,dd}gsum2d_` entry points.
# PBLAS's C code calls the `Cblacs_*` / `C{q,e,dd}gsum2d` counterparts,
# which upstream BLACS produces by compiling every source twice (once
# with -DCallFromC). We replicate that second pass as an OBJECT library
# local to the tests — keeps the change confined to tests/pblas/ and
# leaves the canonical ${LIB_PREFIX}blacs unchanged.
set(_blacs_src_dir ${CMAKE_SOURCE_DIR}/blacs/src)
if(IS_DIRECTORY ${_blacs_src_dir})
    # Upstream BLACS CMakeLists.txt (external/scalapack-2.2.3/BLACS/SRC)
    # compiles this set of "support" files twice — once with CallFromC,
    # once without — so both Cblacs_* and blacs_*_ symbols exist. We
    # match that set here (excluding dcputime00 / dwalltime00 which have
    # no dual-interface gate, and k{send,recv,bs,br}id which are BLACS
    # internal keys, not dual-interface).
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
    add_library(pblas_test_blacs_c_api STATIC
        ${_blacs_c_api_supp} ${_blacs_c_api_comm})
    target_compile_definitions(pblas_test_blacs_c_api PRIVATE CallFromC)
    target_include_directories(pblas_test_blacs_c_api PRIVATE ${_blacs_src_dir})
    if(MPI_C_FOUND)
        target_link_libraries(pblas_test_blacs_c_api PRIVATE MPI::MPI_C)
    endif()
    if(C_AS_CXX)
        # Multifloats target compiles C as C++ (for libmfc operator
        # overloading); match that so the BLACS C sources see the same
        # language settings as ${LIB_PREFIX}blacs does.
        set_source_files_properties(${_blacs_c_api_supp} ${_blacs_c_api_comm}
            PROPERTIES LANGUAGE CXX)
        target_compile_features(pblas_test_blacs_c_api PRIVATE cxx_std_17)
        target_compile_options(pblas_test_blacs_c_api PRIVATE
            -fpermissive -Wno-write-strings)
    endif()
endif()

# ── Common helper modules ───────────────────────────────────────────
add_library(pblas_test_common STATIC
    common/prec_kinds.f90
    common/compare.f90
    common/prec_report.f90
    common/test_data.f90
    common/ref_quad_blas.f90
    common/pblas_grid.f90
    common/pblas_distrib.f90
)
fortran_module_layout(pblas_test_common)
target_link_libraries(pblas_test_common PUBLIC
    refblas_quad
    ${LIB_PREFIX}blacs
    MPI::MPI_Fortran
)
if(TARGET pblas_test_blacs_c_api)
    target_link_libraries(pblas_test_common PUBLIC pblas_test_blacs_c_api)
endif()

# pilaenv.f is the sole Fortran file in PBLAS/SRC (precision-independent
# integer block-size query). The pblas recipe restricts extensions to
# .c/.h, so nothing in the main build ships pilaenv_. PBLAS's C code
# calls it, so we vendor a copy under common/ and compile it here.
add_library(pblas_test_pilaenv STATIC common/pilaenv.f)
set_target_properties(pblas_test_pilaenv PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_options(pblas_test_pilaenv PRIVATE -w -std=legacy)
target_link_libraries(pblas_test_common PUBLIC pblas_test_pilaenv)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/pblas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

set(_wrapper_sources ${_wrapper_dir}/target_pblas.f90)
if(EXISTS ${_wrapper_dir}/ptzblas_stubs.f90)
    # Link-time stubs for routines referenced by the migrated PTZBLAS
    # that the migrator leaves dangling (CALL ZZDOTC etc.). See the
    # stub file for the full rationale.
    list(APPEND _wrapper_sources ${_wrapper_dir}/ptzblas_stubs.f90)
endif()
add_library(pblas_test_target STATIC ${_wrapper_sources})
fortran_module_layout(pblas_test_target)
target_link_libraries(pblas_test_target PUBLIC
    pblas_test_common
    ${LIB_PREFIX}pblas
)
# PBLAS C code calls into PTZBLAS (qvasum_, pairs like qvdot_) and
# PBBLAS (pbqmatadd etc.). The main build doesn't wire these as
# public deps of ${LIB_PREFIX}pblas, so we add them here when present.
foreach(_pblas_dep ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas)
    if(TARGET ${_pblas_dep})
        target_link_libraries(pblas_test_target PUBLIC ${_pblas_dep})
    endif()
endforeach()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(pblas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# Default 4 ranks → 2×2 grid (even for NP=1 this degenerates to 1×1
# and every test still passes, just as serial BLAS). Override at
# configure time: -DPBLAS_TEST_NPROC=<count>.
set(PBLAS_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/pblas (default 4 → 2x2 grid)")

function(add_pblas_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # Turn off LTO for the final link. Ubuntu's default -flto=auto
    # confuses the linker about which archive members satisfy which
    # cross-archive references (we've seen Cblacs_pnum / dwalltime00_
    # reported undefined even when --start-group is in effect). Fat
    # LTO objects still contain native code, so disabling LTO only
    # changes which representation the linker chooses.
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # PBLAS / PTZBLAS / BLACS / C-API shim form a cyclic static-archive
    # graph — e.g. libqpblas references Cblacs_gridinfo (in the shim),
    # the shim references MPI, PTZBLAS references BLAS/BLACS, PBLAS
    # references PTZBLAS. GNU ld in its default single-pass mode
    # leaves some symbols unresolved; LINK_GROUP(RESCAN) wraps the
    # set in --start-group/--end-group so the linker re-scans every
    # member archive until all cross-references settle.
    # LINK_GROUP(RESCAN, ...) only wraps the explicitly-named targets in
    # --start-group/--end-group — transitive deps end up outside the
    # group. We therefore enumerate every archive that participates in
    # the BLACS↔PBLAS↔PTZBLAS dependency cycle.
    set(_pblas_cycle
        pblas_test_target pblas_test_common
        ${LIB_PREFIX}pblas ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas
        ${LIB_PREFIX}blacs ${LIB_PREFIX}blas refblas_quad)
    # pblas_common / blacs_common / ptzblas_common hold the common
    # (precision-independent) halves split out by _install_library_pair;
    # they are referenced inside the cycle too.
    foreach(_common blas_common blacs_common ptzblas_common pblas_common)
        if(TARGET ${_common})
            list(APPEND _pblas_cycle ${_common})
        endif()
    endforeach()
    if(TARGET pblas_test_blacs_c_api)
        list(APPEND _pblas_cycle pblas_test_blacs_c_api)
    endif()
    if(TARGET pblas_test_pilaenv)
        list(APPEND _pblas_cycle pblas_test_pilaenv)
    endif()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_pblas_cycle}>")
    add_test(NAME pblas_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PBLAS_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_lvl level1 level2 level3)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_lvl}/test_*.f90)
    foreach(_t ${_tests})
        add_pblas_precision_test(${_t})
    endforeach()
endforeach()
