# Fortran-only consumer of an installed eplinalg prefix.
#
# LANGUAGES is deliberately just Fortran. This is the consumer shape
# doc/user/ advertises as the normal one, and it is the shape that was
# broken: the generated Config used to call
#
#     find_package(MPI COMPONENTS C)
#
# unconditionally, and FindMPI reports a component NOT-FOUND — "MPI
# component 'C' was requested, but language C is not enabled" — when the
# consuming project() never enabled that language. Every MPI-bound package
# was therefore unfindable from here, while the in-tree build (one
# project() enabling Fortran, C and C++) never noticed.
#
# Adding C to the line below would make this file pass without testing
# anything. Don't.
cmake_minimum_required(VERSION 3.20)
project(consume_fortran LANGUAGES Fortran)

# Compiler-tagged only: exercises the Fortran ABI tag round-trip.
find_package(eplinalgStdBlas REQUIRED)

# MPI-bound: exercises consumer-side MPI flavor detection with no C
# compiler to probe with (the tag has to come from reading mpi.h), and
# the MPI::MPI_C forwarding shim, since this archive's exported interface
# carries $<LINK_ONLY:MPI::MPI_C> that nothing else can resolve here.
find_package(eplinalgStdScalapack REQUIRED)

if(NOT TARGET MPI::MPI_C)
  message(FATAL_ERROR
    "eplinalgStdScalapack's interface references MPI::MPI_C, but no such "
    "target is in scope. The forwarding shim in the package's -deps- file "
    "did not fire.")
endif()

add_executable(consume_fortran main.f90)
target_link_libraries(consume_fortran PRIVATE
  eplinalg::scalapack
  eplinalg::blas)
