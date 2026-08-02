# C-only consumer of an installed eplinalg prefix.
#
# LANGUAGES is deliberately just C, and the package picked is the one
# whose archives are pure C — so nothing drags in libgfortran. That is
# exactly the condition under which the defect this guards is visible:
# METIS calls log/pow/sqrt, nothing declared libm on the exported
# interface, and every other consumer in the tree links through the
# gfortran driver, where libgfortran's own `NEEDED libm.so.6` resolves
# those symbols by accident. Here there is no such accident.
#
# It is also the guard for eplinalgOrdering being a normally generated
# package: it used to ship a hand-written Config assembled with a bare
# include(), bypassing the generator every other package goes through.
#
# Adding Fortran to the line below would make this file pass without
# testing anything. Don't.
cmake_minimum_required(VERSION 3.20)
project(consume_c LANGUAGES C)

find_package(eplinalgOrdering REQUIRED)

add_executable(consume_c main.c)
target_link_libraries(consume_c PRIVATE eplinalg::metis)
