# C-only consumer of a Fortran package, via the -DEPLINALG_FORTRAN_TAG
# escape hatch.
#
# The other two consumers cover the packages a single-language project
# can find on its own. This one covers the case that has no automatic
# answer: the archives carry Fortran objects, so their ABI tag depends on
# a Fortran compiler this project does not have. The consumer names the
# tag itself, and the package then has to supply the Fortran runtime —
# there is no Fortran driver here to do it, and eplinalg::blas drags in
# xerbla, whose WRITE and STOP need _gfortran_st_write and
# _gfortran_stop_string.
#
# Two things must NOT be added below, both of which would make this file
# pass while testing nothing:
#   - Fortran in LANGUAGES. That restores the Fortran driver, and the
#     runtime arrives whether or not the package declares it.
#   - gfortran (or -lgfortran) on target_link_libraries. That is exactly
#     the manual step this path exists to remove.
#
# The tag is passed by the caller (-DEPLINALG_FORTRAN_TAG=gfortran-15),
# not hardcoded, so the same project can be pointed at another build.
cmake_minimum_required(VERSION 3.20)
project(consume_c_fortran_tag LANGUAGES C)

if(NOT DEFINED EPLINALG_FORTRAN_TAG)
  message(FATAL_ERROR
    "This consumer must be configured with -DEPLINALG_FORTRAN_TAG=<tag>; "
    "it has no Fortran compiler to derive one from. That is the point of "
    "the test.")
endif()

find_package(eplinalgStdBlas REQUIRED)

add_executable(consume_c_fortran_tag main.c)
target_link_libraries(consume_c_fortran_tag PRIVATE eplinalg::blas)
