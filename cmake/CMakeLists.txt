# Unified CMake build for migrated extended-precision libraries.
#
# This file is copied into the staging directory by:
#   python -m pyengine stage <staging_dir> --target <target>
#
# Usage:
#   cmake -S <staging_dir> -B <staging_dir>/build -DCMAKE_BUILD_TYPE=Release
#   cmake --build <staging_dir>/build -j
#   cmake --install <staging_dir>/build --prefix /usr/local
#
cmake_minimum_required(VERSION 3.20)

# Load target configuration (written by the Python migration tool)
include(${CMAKE_CURRENT_SOURCE_DIR}/target_config.cmake)

project(${TARGET_NAME}_libraries Fortran C CXX)

# ── Compiler-aware install infrastructure ────────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include(FortranCompiler)

# ── Global compiler flags ────────────────────────────────────────────
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
endif()
set(CMAKE_Fortran_PREPROCESS ON)

# Detect 80-bit extended precision (KIND=10) support.
#
# We probe the literal ``REAL(KIND=10)`` on purpose, not
# ``selected_real_kind(18, 4931)``. On compilers that lack x87 80-bit
# reals (ifx on Linux/x86-64), the latter returns 16 and silently
# aliases the KIND=10 path onto KIND=16, which makes the EISNAN and
# QISNAN overloads in LA_XISNAN_EP collapse to the same argument type
# and triggers an ambiguous-generic-interface error. Asking for
# KIND=10 directly fails to compile on such compilers, which is the
# behavior we want.
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("
    program test_real10
    real(kind=10) :: x
    x = 1.0_10
    end program
" HAVE_REAL10 SRC_EXT F90)
if(HAVE_REAL10)
    add_compile_definitions(HAVE_REAL10)
endif()

# Detect 128-bit quadruple precision (IEEE binary128) support.
#
# We probe ``selected_real_kind(33, 4931)`` (the pattern used inside
# la_constants_ep to pick ``qp``), not the literal ``REAL(KIND=16)``.
# LLVM Flang 20 accepts the KIND=16 literal but then returns -1 for
# ``selected_real_kind(33)``, so the module ends up with ``qp = -1``
# and every subsequent ``real(qp)`` declaration fails. Probing the
# real selector catches this case; the module's qp block is then
# gated behind HAVE_REAL16 for flang's benefit.
check_fortran_source_compiles("
    program test_real16
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp) :: x
    x = 1.0_qp
    end program
" HAVE_REAL16 SRC_EXT F90)
if(HAVE_REAL16)
    add_compile_definitions(HAVE_REAL16)
endif()

# ── MPI (needed by BLACS, PBLAS, ScaLAPACK C) ───────────────────────
find_package(MPI COMPONENTS C CXX Fortran QUIET)

# ── Multifloats (conditional) ────────────────────────────────────────
if(NEEDS_MULTIFLOATS)
    find_program(FYPP_EXECUTABLE fypp REQUIRED)

    # Use the multifloats top-level directory (not src/ and fsrc/
    # separately): upstream's src/CMakeLists.txt uses CMAKE_SOURCE_DIR
    # which is only correct when multifloats is the top-level project.
    # EXCLUDE_FROM_ALL keeps us from building multifloats's install
    # rules / tests by default.
    #
    # BUILD_TESTING/MULTIFLOATS_BUILD_BENCH are set as plain (non-cache)
    # variables so they only suppress tests *inside* the multifloats
    # sub-build. A `CACHE … FORCE` here would also disable our own
    # `tests/` directory at the end of this file.
    set(BUILD_TESTING OFF)
    set(MULTIFLOATS_BUILD_BENCH OFF)
    if(DEFINED MULTIFLOATS_DIR)
        message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
        add_subdirectory(${MULTIFLOATS_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
    else()
        include(FetchContent)
        set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
            CACHE STRING "Git URL for the multifloats library")
        set(MULTIFLOATS_GIT_TAG "main"
            CACHE STRING "Git tag/branch/commit for multifloats")
        message(STATUS "Fetching multifloats from ${MULTIFLOATS_GIT_REPO} (${MULTIFLOATS_GIT_TAG})")
        FetchContent_Declare(multifloats_fetch
            GIT_REPOSITORY ${MULTIFLOATS_GIT_REPO}
            GIT_TAG        ${MULTIFLOATS_GIT_TAG}
        )
        FetchContent_Populate(multifloats_fetch)
        # Upstream src/CMakeLists.txt uses CMAKE_SOURCE_DIR/include which
        # only resolves correctly when multifloats is the top-level build.
        # Patch it to PROJECT_SOURCE_DIR (the multifloats sub-project root)
        # so the add_subdirectory inclusion below can find float64x2.h.
        file(READ "${multifloats_fetch_SOURCE_DIR}/src/CMakeLists.txt" _mf_src_cmake)
        string(REPLACE
            "\${CMAKE_SOURCE_DIR}/include"
            "\${PROJECT_SOURCE_DIR}/../include"
            _mf_src_cmake "${_mf_src_cmake}")
        file(WRITE "${multifloats_fetch_SOURCE_DIR}/src/CMakeLists.txt" "${_mf_src_cmake}")
        # Same symptom in fsrc/CMakeLists.txt: ``${CMAKE_SOURCE_DIR}/cmake/
        # multifloatsfConfig.cmake.in`` resolves against our staging dir,
        # not the fetched multifloats source root. Redirect to the
        # multifloats sub-project root so configure_package_config_file
        # finds the template where FetchContent put it.
        if(EXISTS "${multifloats_fetch_SOURCE_DIR}/fsrc/CMakeLists.txt")
            file(READ "${multifloats_fetch_SOURCE_DIR}/fsrc/CMakeLists.txt"
                 _mf_fsrc_cmake)
            string(REPLACE
                "\${CMAKE_SOURCE_DIR}/cmake/"
                "\${PROJECT_SOURCE_DIR}/../cmake/"
                _mf_fsrc_cmake "${_mf_fsrc_cmake}")
            file(WRITE "${multifloats_fetch_SOURCE_DIR}/fsrc/CMakeLists.txt"
                 "${_mf_fsrc_cmake}")
        endif()
        add_subdirectory(
            ${multifloats_fetch_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
    endif()

    # la_constants_mf helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_constants_mf.f90")
        add_library(la_constants_mf STATIC "${MF_HELPERS_DIR}/la_constants_mf.f90")
        fortran_module_layout(la_constants_mf)
        target_link_libraries(la_constants_mf PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()

    # la_xisnan_mf helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        add_library(la_xisnan_mf STATIC "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        fortran_module_layout(la_xisnan_mf)
        target_link_libraries(la_xisnan_mf PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()

    # multifloats_mpi — C++ MPI handle registration + bridge header
    if(EXISTS "${MF_HELPERS_DIR}/multifloats_mpi.cpp" AND MPI_CXX_FOUND)
        add_library(multifloats_mpi STATIC "${MF_HELPERS_DIR}/multifloats_mpi.cpp")
        target_include_directories(multifloats_mpi PUBLIC
            $<BUILD_INTERFACE:${MF_HELPERS_DIR}>)
        target_link_libraries(multifloats_mpi PUBLIC MPI::MPI_CXX multifloats)
        target_compile_features(multifloats_mpi PRIVATE cxx_std_17)
    elseif(EXISTS "${MF_HELPERS_DIR}/multifloats_bridge.h")
        # No MPI — still provide the bridge header via an interface target
        add_library(multifloats_mpi INTERFACE)
        target_include_directories(multifloats_mpi INTERFACE
            $<BUILD_INTERFACE:${MF_HELPERS_DIR}>)
    endif()
endif()

# ── Helper functions ─────────────────────────────────────────────────

function(add_migrated_fortran_library lib_name)
    set(_dir ${CMAKE_CURRENT_SOURCE_DIR}/${lib_name})
    if(NOT EXISTS "${_dir}/manifest.cmake")
        return()
    endif()
    include(${_dir}/manifest.cmake)

    # Prefix source paths with library directory
    set(_common "")
    foreach(f ${${lib_name}_COMMON_SOURCES})
        list(APPEND _common "${_dir}/${f}")
    endforeach()
    set(_precision "")
    foreach(f ${${lib_name}_PRECISION_SOURCES})
        list(APPEND _precision "${_dir}/${f}")
    endforeach()

    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        fortran_module_layout(${_common_target})
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    fortran_module_layout(${_precision_target})
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()

    # Link multifloats helpers for Fortran libraries. multifloatsf
    # comes from a fetched sub-project that owns its own install/export
    # set; wrap the link in ``$<BUILD_INTERFACE:>`` so the dependency
    # is visible during the build (tests, locally-built downstreams)
    # but isn't propagated into our install(EXPORT) — otherwise CMake
    # rejects the generate step because multifloatsf is "not in any
    # export set" known to us.
    if(NEEDS_MULTIFLOATS)
        target_link_libraries(${_precision_target} PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
        if(TARGET la_constants_mf)
            target_link_libraries(${_precision_target} PUBLIC la_constants_mf)
        endif()
        if(TARGET la_xisnan_mf)
            target_link_libraries(${_precision_target} PUBLIC la_xisnan_mf)
        endif()
    endif()
endfunction()


function(add_migrated_c_library lib_name)
    set(_dir ${CMAKE_CURRENT_SOURCE_DIR}/${lib_name})
    if(NOT EXISTS "${_dir}/manifest.cmake")
        return()
    endif()
    include(${_dir}/manifest.cmake)

    # Prefix source paths with library directory
    set(_common "")
    foreach(f ${${lib_name}_COMMON_SOURCES})
        list(APPEND _common "${_dir}/${f}")
    endforeach()
    set(_precision "")
    foreach(f ${${lib_name}_PRECISION_SOURCES})
        list(APPEND _precision "${_dir}/${f}")
    endforeach()
    set(_dual "")
    foreach(f ${${lib_name}_DUAL_INTERFACE_SOURCES})
        list(APPEND _dual "${_dir}/${f}")
    endforeach()

    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        target_include_directories(${_common_target} PUBLIC $<BUILD_INTERFACE:${_dir}/src>)
        if(MPI_C_FOUND)
            # Link MPI privately (needed for compiling the C sources) and
            # propagate it via INTERFACE link-only, so downstream Fortran
            # targets inherit the link dependency without MPI's C-only
            # compile flags (e.g. mpich's ``-ffat-lto-objects`` is rejected
            # by flang-new-20).
            target_link_libraries(${_common_target} PRIVATE MPI::MPI_C)
            target_link_libraries(${_common_target} INTERFACE $<LINK_ONLY:MPI::MPI_C>)
        endif()
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    target_include_directories(${_precision_target} PUBLIC $<BUILD_INTERFACE:${_dir}/src>)
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()
    if(MPI_C_FOUND)
        target_link_libraries(${_precision_target} PRIVATE MPI::MPI_C)
        target_link_libraries(${_precision_target} INTERFACE $<LINK_ONLY:MPI::MPI_C>)
    endif()

    # Dual-interface compile: files listed in
    # ``${lib_name}_DUAL_INTERFACE_SOURCES`` are already compiled once
    # above (they live in _precision / _common, since every such file
    # is also one of BLACS's regular sources) and expose the
    # Fortran-callable entry points by default. Compile them a second
    # time as an OBJECT library with ``-DCallFromC`` so the C-callable
    # entry points (``Cblacs_*``, ``Cdgesd2d`` etc.) are also present
    # in the final static archive. Empty ``_dual`` list is a no-op.
    if(_dual)
        set(_dual_target "${_precision_target}_c_api")
        add_library(${_dual_target} OBJECT ${_dual})
        target_compile_definitions(${_dual_target} PRIVATE CallFromC)
        target_include_directories(${_dual_target} PRIVATE $<BUILD_INTERFACE:${_dir}/src>)
        if(MPI_C_FOUND)
            target_link_libraries(${_dual_target} PRIVATE MPI::MPI_C)
        endif()
        if(C_AS_CXX)
            set_source_files_properties(${_dual} PROPERTIES LANGUAGE CXX)
            target_compile_features(${_dual_target} PRIVATE cxx_std_17)
            target_compile_options(${_dual_target} PRIVATE -fpermissive -Wno-write-strings)
            # Multifloats: ``Bdef.h`` pulls in ``multifloats_bridge.h`` under
            # multifloats mode. That header is exposed through the
            # ``multifloats_mpi`` target (INTERFACE or STATIC depending on
            # MPI availability); propagate its include directories so the
            # dual-interface compile finds it just like the main target.
            if(TARGET multifloats_mpi)
                target_link_libraries(${_dual_target} PRIVATE multifloats_mpi)
            endif()
        endif()
        target_sources(${_precision_target}
                       PRIVATE $<TARGET_OBJECTS:${_dual_target}>)
    endif()

    # Multifloats: compile C as C++ for operator overloading support.
    # ``copy_files``-staged Fortran helpers (e.g. pblas/pilaenv.f) ride in
    # COMMON_SOURCES but must retain their native language; filter them
    # out of the language override.
    if(C_AS_CXX)
        set(_all_sources "")
        foreach(_s ${_common} ${_precision})
            get_filename_component(_ext "${_s}" EXT)
            string(TOLOWER "${_ext}" _ext)
            if(NOT _ext STREQUAL ".f" AND NOT _ext STREQUAL ".f90"
                AND NOT _ext STREQUAL ".f90" AND NOT _ext STREQUAL ".f03"
                AND NOT _ext STREQUAL ".for")
                list(APPEND _all_sources "${_s}")
            endif()
        endforeach()
        set_source_files_properties(${_all_sources} PROPERTIES LANGUAGE CXX)
        if(TARGET ${_common_target})
            target_compile_features(${_common_target} PRIVATE cxx_std_17)
            target_compile_options(${_common_target} PRIVATE -fpermissive -Wno-write-strings)
        endif()
        target_compile_features(${_precision_target} PRIVATE cxx_std_17)
        target_compile_options(${_precision_target} PRIVATE -fpermissive -Wno-write-strings)
        # Link multifloats + bridge for C libraries in multifloats mode
        if(TARGET multifloats)
            target_link_libraries(${_precision_target} PUBLIC multifloats)
        endif()
        if(TARGET multifloats_mpi)
            target_link_libraries(${_precision_target} PUBLIC multifloats_mpi)
            if(TARGET ${_common_target})
                target_link_libraries(${_common_target} PUBLIC multifloats_mpi)
            endif()
        endif()
    endif()
endfunction()

# ── Libraries in dependency order ────────────────────────────────────

# 1. BLAS (no deps)
add_migrated_fortran_library(blas)

# 2. BLACS (no library deps, needs MPI)
add_migrated_c_library(blacs)
if(TARGET ${LIB_PREFIX}blacs AND NEEDS_MULTIFLOATS AND TARGET multifloats_mpi)
    target_link_libraries(${LIB_PREFIX}blacs PUBLIC multifloats_mpi)
endif()

# 3. LAPACK (depends on BLAS)
add_migrated_fortran_library(lapack)
if(TARGET ${LIB_PREFIX}lapack AND TARGET ${LIB_PREFIX}blas)
    target_link_libraries(${LIB_PREFIX}lapack PUBLIC ${LIB_PREFIX}blas)
endif()

# 4. PTZBLAS (depends on BLAS)
add_migrated_fortran_library(ptzblas)
if(TARGET ${LIB_PREFIX}ptzblas AND TARGET ${LIB_PREFIX}blas)
    target_link_libraries(${LIB_PREFIX}ptzblas PUBLIC ${LIB_PREFIX}blas)
endif()

# 5. PBBLAS (depends on BLAS, BLACS)
add_migrated_fortran_library(pbblas)
if(TARGET ${LIB_PREFIX}pbblas)
    foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs)
        if(TARGET ${_dep})
            target_link_libraries(${LIB_PREFIX}pbblas PUBLIC ${_dep})
        endif()
    endforeach()
endif()

# 6. PBLAS (depends on BLAS, BLACS, PTZBLAS, PBBLAS; C library)
# PBLAS C entry points call into PTZBLAS kernels (e.g. qvasum_, qvvdotc_)
# and PBBLAS kernels (pbqmatadd_, pbqtran_) at runtime, so both must be
# in the link graph alongside BLAS and BLACS.
add_migrated_c_library(pblas)
if(TARGET ${LIB_PREFIX}pblas)
    foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs
                 ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas)
        if(TARGET ${_dep})
            target_link_libraries(${LIB_PREFIX}pblas PUBLIC ${_dep})
        endif()
    endforeach()
endif()

# 7. ScaLAPACK (depends on LAPACK, BLACS, PBLAS)
add_migrated_fortran_library(scalapack)
if(TARGET ${LIB_PREFIX}scalapack)
    foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas)
        if(TARGET ${_dep})
            target_link_libraries(${LIB_PREFIX}scalapack PUBLIC ${_dep})
        endif()
    endforeach()
endif()

# 8. ScaLAPACK C (depends on LAPACK, BLACS, PBLAS, ScaLAPACK).
# Skip for multifloats: scalapack_c includes REDIST/SRC sources whose
# K&R-style function definitions don't survive C-as-C++ compilation.
# Patching every K&R def in those files is the proper fix (project
# policy: assume ANSI C); until that lands, the multifloats target
# simply does not build scalapack_c. tests/scalapack ships no
# target_multifloats wrapper today, so nothing depends on the
# resulting symbols on this target.
if(NOT C_AS_CXX)
    add_migrated_c_library(scalapack_c)
endif()
if(TARGET ${LIB_PREFIX}scalapack_c)
    # REDIST/SRC sources select their Fortran-call entry-point naming
    # convention (trailing underscore vs UpCase vs no-change) at compile
    # time via ``Add_`` / ``UpCase`` / ``NoChange``. Upstream's
    # SLmake.inc.example sets ``-DAdd_`` to get the lowercase + trailing
    # underscore convention that gfortran emits, so we match.
    target_compile_definitions(${LIB_PREFIX}scalapack_c PRIVATE Add_)
    if(TARGET scalapack_c_common)
        target_compile_definitions(scalapack_c_common PRIVATE Add_)
    endif()
    foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas ${LIB_PREFIX}scalapack)
        if(TARGET ${_dep})
            target_link_libraries(${LIB_PREFIX}scalapack_c PUBLIC ${_dep})
        endif()
    endforeach()
    # Reverse dep: ScaLAPACK's Fortran routines call qlamov_/xlamov_
    # (and the redistribution entry points pqgemr2d_/pxgemr2d_) which
    # live in ${LIB_PREFIX}scalapack_c. Downstream test link cycles
    # then need only to name ${LIB_PREFIX}scalapack to pick both up.
    if(TARGET ${LIB_PREFIX}scalapack)
        target_link_libraries(${LIB_PREFIX}scalapack PUBLIC
            ${LIB_PREFIX}scalapack_c)
    endif()
endif()

# ── Install ──────────────────────────────────────────────────────────
# Each precision library gets its own export set and Config.cmake,
# so consumers can do: find_package(qblas), find_package(qlapack), etc.
# The common library is bundled into the same export as its precision library.

function(_install_library_pair lib_name)
    set(_precision_target "${LIB_PREFIX}${lib_name}")
    set(_common_target "${lib_name}_common")
    set(_export "${_precision_target}Targets")

    if(NOT TARGET ${_precision_target})
        return()
    endif()

    # Install precision-specific library
    fortran_install_library(${_precision_target} EXPORT ${_export})
    get_target_property(_fc_moddir ${_precision_target} Fortran_MODULE_DIRECTORY)
    if(_fc_moddir)
        fortran_install_modules(${_precision_target})
    endif()

    # Bundle common library into the same export
    if(TARGET ${_common_target})
        fortran_install_library(${_common_target} EXPORT ${_export})
        get_target_property(_fc_moddir ${_common_target} Fortran_MODULE_DIRECTORY)
        if(_fc_moddir)
            fortran_install_modules(${_common_target})
        endif()
    endif()
endfunction()

foreach(_lib ${STAGED_LIBRARIES})
    _install_library_pair(${_lib})
endforeach()

if(NEEDS_MULTIFLOATS)
    # multifloats defines its own install rules; only install local helpers
    foreach(_t la_constants_mf la_xisnan_mf multifloats_mpi)
        if(TARGET ${_t})
            fortran_install_library(${_t} EXPORT ${_t}Targets)
            get_target_property(_fc_moddir ${_t} Fortran_MODULE_DIRECTORY)
            if(_fc_moddir)
                fortran_install_modules(${_t})
            endif()
        endif()
    endforeach()
endif()

# ── Tests ────────────────────────────────────────────────────────────
# Differential precision test suite. Each library has its own subdir
# under tests/ which is added when its tests/<lib>/CMakeLists.txt is
# present in the staging copy.
#
# We deliberately do NOT gate on BUILD_TESTING: the multifloats
# sub-build sets it locally to OFF to suppress its own tests, and
# that local value shadows the cache value from include(CTest) at
# the parent scope. Always-on with a directory check gives us the
# single signal "is the staging dir self-contained?" — drop tests/
# from staging to skip them.
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    enable_testing()
    add_subdirectory(tests)
endif()
