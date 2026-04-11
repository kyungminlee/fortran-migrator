cmake_minimum_required(VERSION 3.10)
project(multifloats LANGUAGES Fortran C)

add_library(multifloats multifloats.f90)
target_include_directories(multifloats PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# --------------------------------------------------------------------- #
# libmfc -- C companion to the Fortran multifloats module.
#
# Provides plain-C struct types matching the Fortran SEQUENCE layout and
# runtime MPI datatype/op registration used by the migrated BLACS and
# PBLAS C sources. Required only by the multifloats target of the
# migrator; the Fortran-only targets do not use it.
# --------------------------------------------------------------------- #
find_package(MPI COMPONENTS C)

if(MPI_C_FOUND)
    add_library(mfc STATIC src/multifloats_mpi.c)
    target_include_directories(mfc PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(mfc PUBLIC MPI::MPI_C)
    set_target_properties(mfc PROPERTIES
        C_STANDARD 99
        C_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE ON
    )
else()
    message(STATUS "MPI not found; libmfc (C companion to multifloats) "
                   "will not be built. BLACS/PBLAS multifloats migration "
                   "requires MPI.")
endif()
