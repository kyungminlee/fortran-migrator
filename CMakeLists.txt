cmake_minimum_required(VERSION 3.27)
project(fortran-migrator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(FORTRAN_MIGRATOR_USE_FLANG
       "Link against LLVM/Flang parser libraries for AST-based analysis" OFF)
option(FORTRAN_MIGRATOR_BUILD_TESTS
       "Build tests" ON)

# ---------------------------------------------------------------------------
# Core library (everything except main)
# ---------------------------------------------------------------------------
add_library(fortran-migrator-core STATIC
    src/symbol_database.cpp
    src/source_scanner.cpp
    src/transformation.cpp
    src/rewriter.cpp
    src/fortran_format.cpp
)
target_include_directories(fortran-migrator-core PUBLIC src)

# ---------------------------------------------------------------------------
# Flang parser integration (optional)
# ---------------------------------------------------------------------------
if(FORTRAN_MIGRATOR_USE_FLANG)
    find_package(LLVM REQUIRED CONFIG)
    message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
    message(STATUS "LLVM install prefix: ${LLVM_INSTALL_PREFIX}")

    # Flang parser/support/common libraries ship alongside LLVM when Flang
    # is built as part of the monorepo.  The library names changed between
    # LLVM versions (FortranParser → flangParser, etc.).
    foreach(_component Parser Common Support)
        string(TOUPPER ${_component} _UPPER)
        find_library(FLANG_${_UPPER}_LIB
            NAMES Fortran${_component} flang${_component}
            HINTS ${LLVM_LIBRARY_DIRS}
        )
        if(NOT FLANG_${_UPPER}_LIB)
            message(FATAL_ERROR
                "Could not find Flang ${_component} library in ${LLVM_LIBRARY_DIRS}")
        endif()
        message(STATUS "Flang ${_component}: ${FLANG_${_UPPER}_LIB}")
    endforeach()

    target_compile_definitions(fortran-migrator-core PUBLIC USE_FLANG_PARSER)
    target_include_directories(fortran-migrator-core PUBLIC ${LLVM_INCLUDE_DIRS})
    target_link_libraries(fortran-migrator-core PUBLIC
        ${FLANG_PARSER_LIB}
        ${FLANG_COMMON_LIB}
        ${FLANG_SUPPORT_LIB}
    )
endif()

# ---------------------------------------------------------------------------
# Main executable
# ---------------------------------------------------------------------------
add_executable(fortran-migrator src/main.cpp)
target_link_libraries(fortran-migrator PRIVATE fortran-migrator-core)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(FORTRAN_MIGRATOR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
