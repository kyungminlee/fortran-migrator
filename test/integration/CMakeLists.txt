# Top-level test directory — added by cmake/CMakeLists.txt when
# BUILD_TESTING=ON. Each library has its own subdir; we only descend
# into ones whose precision library was actually built (i.e., the
# library was migrated for this target) AND whose CMakeLists.txt
# exists in the staging copy.

# fypp is used to expand the per-framework target_*_body.fypp template
# into per-target target_*.f90 sources. Available via `uv run fypp`
# from the project venv.
find_program(FYPP fypp REQUIRED)

# Helper for per-target wrapper shims:
#   add_fypp_source(<out_var> <shim.fypp> <body_dir>)
# generates <build>/.../<basename>.f90 from <shim.fypp>, with -I
# pointing at <body_dir> so #:include resolves the shared body file.
function(add_fypp_source out_var shim_path body_dir)
    get_filename_component(_stem ${shim_path} NAME_WE)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/${_stem}.f90)
    # Depend on every .fypp in body_dir so edits to the shared body
    # (or any future #:include sibling) trigger regeneration.
    file(GLOB _body_fypps CONFIGURE_DEPENDS ${body_dir}/*.fypp)
    add_custom_command(
        OUTPUT ${_out}
        COMMAND ${FYPP} -I ${body_dir} ${shim_path} ${_out}
        DEPENDS ${shim_path} ${_body_fypps}
        COMMENT "fypp ${shim_path} -> ${_out}"
        VERBATIM
    )
    set(${out_var} ${_out} PARENT_SCOPE)
endfunction()

# Export-rename step for the quad reference archives:
#   refquad_rename_exports(<target> [SIBLING <archive-target>])
#
# Renames <target>'s Fortran-mangled exports from `<name>_` to
# `<name>_quad_` (scripts/refquad_rename_archive.sh as a POST_BUILD
# step) so the quad reference doesn't collide with the standard
# ``blas`` / ``lapack`` archive's native-precision symbols when the
# kind4 / kind8 baseline tests link both. The consumer side lives in
# tests/common/ref_quad_*.f90: every routine in the interface block is
# suffixed ``_quad``, with a thin module procedure under the original
# name forwarding to the renamed external, so test files keep
# ``use ref_quad_blas, only: dasum`` unchanged.
#
# SIBLING names a second archive whose defs join the rename set — and
# which <target> is made to depend on. reflapack_quad needs it: its U
# refs to dnrm2_ etc. must be renamed in lockstep with refblas_quad's
# defs of the same names.
#
# The `-x` excludes are the migrated ``*_common`` archives. Their T/W
# defs (la_constants_ep, iparam2stage_ep_, mumps_abort_, …) are
# referenced by the quad references and must stay un-renamed, or the
# migrated link cycle can no longer resolve them.
#
# The script is looked for both in-tree (scripts/) and one level up —
# the staging tree, where ``migrator stage`` plants the cmake/ files
# alongside ``scripts/``. Skipped with a warning if neither resolves.
function(refquad_rename_exports target)
    cmake_parse_arguments(PARSE_ARGV 1 _rqr "" "SIBLING" "")

    find_file(REFQUAD_RENAME_SCRIPT refquad_rename_archive.sh
        PATHS
            ${CMAKE_SOURCE_DIR}/scripts
            ${CMAKE_SOURCE_DIR}/../scripts
            ${PROJECT_SOURCE_DIR}/scripts
            ${PROJECT_SOURCE_DIR}/../scripts
        NO_DEFAULT_PATH)
    if(NOT REFQUAD_RENAME_SCRIPT)
        message(WARNING "${target}: rename script not found; "
                        "kind4/kind8 baseline tests will see symbol collisions")
        return()
    endif()

    set(_sibling_archive)
    if(NOT "${_rqr_SIBLING}" STREQUAL "")
        set(_sibling_archive $<TARGET_FILE:${_rqr_SIBLING}>)
    endif()

    set(_refquad_excludes)
    foreach(_c blas_common lapack_common scalapack_common
               pblas_common pbblas_common ptzblas_common blacs_common
               mumps_common xblas_common)
        if(TARGET ${_c})
            list(APPEND _refquad_excludes -x $<TARGET_FILE:${_c}>)
        endif()
    endforeach()

    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND sh ${REFQUAD_RENAME_SCRIPT}
                $<TARGET_FILE:${target}>
                ${_sibling_archive}
                ${_refquad_excludes}
        COMMENT "${target}: rename Fortran exports -> *_quad_"
        VERBATIM)

    if(NOT "${_rqr_SIBLING}" STREQUAL "")
        add_dependencies(${target} ${_rqr_SIBLING})
    endif()
endfunction()

# Disable the free-form line-length limit for the whole test tree (this
# directory and every add_subdirectory below it). Test sources are free
# form (.f90) and some legitimately run past column 132 — ASCII argument
# tables (e.g. blas/level3/test_zgemmtr.f90) that hard-error when the cut
# lands inside a character constant, and UTF-8 box-drawing comments whose
# bytes overflow 132 and get silently truncated. The migrated libraries
# get the analogous fixed-form treatment in cmake/CMakeLists.txt; this is
# the free-form counterpart, scoped to tests so library binaries are
# untouched. Intel free form is unbounded by default, so no flag needed.
if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran"
   OR FORTRAN_COMPILER_FAMILY STREQUAL "flang")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()

# Task 44: extended stagings link the ep_-privatized BLACS/ScaLAPACK-
# common engine, whose family-independent entry points are ep_-prefixed
# so they never collide with an MKL in the same consumer. Preprocessed
# test sources (.F90) include tests/common/privatized_blacs_names.h,
# which maps the pristine names onto the ep_ twins when this macro is
# defined. Baseline stagings write an empty PRIVATIZED_LIBRARIES and
# keep pristine names.
if("blacs" IN_LIST PRIVATIZED_LIBRARIES)
    add_compile_definitions(EPLINALG_TEST_EP_BLACS)
endif()

# Per-target conversion library shared by every framework's test wrappers.
add_subdirectory(common)

if(BASELINE_BUILD)
    # Baseline (kind4/kind8): no migrated archives exist; STAGED_LIBRARIES
    # is empty by design. Iterate the on-disk test framework dirs and
    # descend into each whose corresponding std archive built. blas/lapack
    # always; pblas/scalapack only when MPI was found.
    foreach(_lib blas lapack pblas scalapack)
        set(_libdir ${CMAKE_CURRENT_SOURCE_DIR}/${_lib})
        if(TARGET ${_lib} AND EXISTS ${_libdir}/CMakeLists.txt)
            add_subdirectory(${_lib})
        endif()
    endforeach()
else()
    # Non-baseline builds must carry a non-empty LIB_PAIR_PREFIX; otherwise
    # the gate below would collapse to ``TARGET <lib>`` and accidentally
    # admit the std archives, running each test tree's non-baseline
    # branch against the wrong (standard-precision) target.
    if(LIB_PAIR_PREFIX STREQUAL "")
        message(FATAL_ERROR
            "tests/CMakeLists.txt: LIB_PAIR_PREFIX is empty in a "
            "non-baseline build — refusing to descend (would run "
            "migrated test branches against std archives). Check "
            "target_config.cmake and the active TARGET_NAME.")
    endif()
    foreach(_lib ${STAGED_LIBRARIES})
        set(_libdir ${CMAKE_CURRENT_SOURCE_DIR}/${_lib})
        if(TARGET ${LIB_PAIR_PREFIX}${_lib}
                AND EXISTS ${_libdir}/CMakeLists.txt)
            add_subdirectory(${_lib})
        endif()
    endforeach()

endif()
