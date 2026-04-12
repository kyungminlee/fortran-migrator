cmake_minimum_required(VERSION 3.10)
project(multifloats LANGUAGES Fortran CXX)

add_library(multifloats multifloats.f90)
target_include_directories(multifloats PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# --------------------------------------------------------------------- #
# multifloats_mpi -- C++ companion for MPI handle registration.
#
# Provides runtime MPI_FLOAT64X2 / MPI_COMPLEX128X2 datatype and
# MPI_DD_SUM / MPI_ZZ_SUM / MPI_DD_AMX / MPI_DD_AMN / MPI_ZZ_AMX /
# MPI_ZZ_AMN op registration. Compiled as C++ (uses multifloats.hh
# operator overloads in the user-op callbacks). Linked into both the
# migrated BLACS and PBLAS archives.
# --------------------------------------------------------------------- #
find_package(MPI COMPONENTS CXX)

if(MPI_CXX_FOUND)
    add_library(multifloats_mpi STATIC src/multifloats_mpi.cpp)
    target_include_directories(multifloats_mpi PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/..    # for multifloats.hh
    )
    target_link_libraries(multifloats_mpi PUBLIC MPI::MPI_CXX)
    set_target_properties(multifloats_mpi PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE ON
    )
else()
    message(STATUS "MPI not found; multifloats_mpi (C++ companion) "
                   "will not be built. BLACS/PBLAS multifloats migration "
                   "requires MPI.")
endif()
