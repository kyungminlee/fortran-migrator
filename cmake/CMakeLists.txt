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

# Detect 80-bit extended precision (KIND=10) support
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("
    program test_real10
    integer, parameter :: ep = selected_real_kind(18, 4931)
    real(ep) :: x
    x = 1.0_ep
    end program
" HAVE_REAL10 SRC_EXT F90)
if(HAVE_REAL10)
    add_compile_definitions(HAVE_REAL10)
endif()

# ── MPI (needed by BLACS, PBLAS, ScaLAPACK C) ───────────────────────
find_package(MPI COMPONENTS C CXX Fortran QUIET)

# ── Multifloats (conditional) ────────────────────────────────────────
if(NEEDS_MULTIFLOATS)
    if(DEFINED MULTIFLOATS_DIR)
        # Use a local directory (for offline builds or development)
        message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
        add_subdirectory(${MULTIFLOATS_DIR} ${CMAKE_CURRENT_BINARY_DIR}/_multifloats)
    else()
        # Fetch from GitHub
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
        FetchContent_MakeAvailable(multifloats_fetch)
    endif()

    # la_constants_mf helper module
    if(EXISTS "${MF_HELPERS_DIR}/la_constants_mf.f90")
        add_library(la_constants_mf STATIC "${MF_HELPERS_DIR}/la_constants_mf.f90")
        fortran_module_layout(la_constants_mf)
        target_link_libraries(la_constants_mf PUBLIC multifloats)
    endif()

    # la_xisnan_mf helper module
    if(EXISTS "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        add_library(la_xisnan_mf STATIC "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        fortran_module_layout(la_xisnan_mf)
        target_link_libraries(la_xisnan_mf PUBLIC multifloats)
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

    # Link multifloats helpers for Fortran libraries
    if(NEEDS_MULTIFLOATS)
        target_link_libraries(${_precision_target} PUBLIC multifloats)
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

    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        target_include_directories(${_common_target} PUBLIC ${_dir}/src)
        if(MPI_C_FOUND)
            target_link_libraries(${_common_target} PUBLIC MPI::MPI_C)
        endif()
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    target_include_directories(${_precision_target} PUBLIC ${_dir}/src)
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()
    if(MPI_C_FOUND)
        target_link_libraries(${_precision_target} PUBLIC MPI::MPI_C)
    endif()

    # Multifloats: compile C as C++ for operator overloading support
    if(C_AS_CXX)
        set(_all_sources ${_common} ${_precision})
        set_source_files_properties(${_all_sources} PROPERTIES LANGUAGE CXX)
        if(TARGET ${_common_target})
            target_compile_features(${_common_target} PRIVATE cxx_std_17)
            target_compile_options(${_common_target} PRIVATE -fpermissive -Wno-write-strings)
        endif()
        target_compile_features(${_precision_target} PRIVATE cxx_std_17)
        target_compile_options(${_precision_target} PRIVATE -fpermissive -Wno-write-strings)
        # Link multifloats for C libraries in multifloats mode
        if(TARGET multifloats)
            target_link_libraries(${_precision_target} PUBLIC multifloats)
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

# 6. PBLAS (depends on BLAS, BLACS; C library)
add_migrated_c_library(pblas)
if(TARGET ${LIB_PREFIX}pblas)
    foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs)
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

# 8. ScaLAPACK C (depends on LAPACK, BLACS, PBLAS, ScaLAPACK)
add_migrated_c_library(scalapack_c)
if(TARGET ${LIB_PREFIX}scalapack_c)
    foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas ${LIB_PREFIX}scalapack)
        if(TARGET ${_dep})
            target_link_libraries(${LIB_PREFIX}scalapack_c PUBLIC ${_dep})
        endif()
    endforeach()
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
    foreach(_t multifloats la_constants_mf la_xisnan_mf)
        if(TARGET ${_t})
            fortran_install_library(${_t} EXPORT ${_t}Targets)
            get_target_property(_fc_moddir ${_t} Fortran_MODULE_DIRECTORY)
            if(_fc_moddir)
                fortran_install_modules(${_t})
            endif()
        endif()
    endforeach()
endif()
