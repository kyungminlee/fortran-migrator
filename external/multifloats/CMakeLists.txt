cmake_minimum_required(VERSION 3.10)
project(multifloats LANGUAGES Fortran)

add_library(multifloats multifloats.f90)
target_include_directories(multifloats PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
