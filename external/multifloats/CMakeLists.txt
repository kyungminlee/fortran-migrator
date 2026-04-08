cmake_minimum_required(VERSION 3.10)
project(multifloats LANGUAGES Fortran)

add_library(multifloats multifloats.f90)
target_include_directories(multifloats PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_executable(multifloats_unittest multifloats_unittest.f90)
target_link_libraries(multifloats_unittest PRIVATE multifloats)

enable_testing()
add_test(NAME multifloats_unittest COMMAND multifloats_unittest)
