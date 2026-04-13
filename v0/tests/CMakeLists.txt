# Test executables link against the core library.

add_executable(test_symbol_database test_symbol_database.cpp)
target_link_libraries(test_symbol_database PRIVATE fortran-migrator-core)
add_test(NAME symbol_database COMMAND test_symbol_database)

add_executable(test_transformation test_transformation.cpp)
target_link_libraries(test_transformation PRIVATE fortran-migrator-core)
add_test(NAME transformation COMMAND test_transformation)

add_executable(test_fortran_format test_fortran_format.cpp)
target_link_libraries(test_fortran_format PRIVATE fortran-migrator-core)
add_test(NAME fortran_format COMMAND test_fortran_format)

add_executable(test_rewriter test_rewriter.cpp)
target_link_libraries(test_rewriter PRIVATE fortran-migrator-core)
add_test(NAME rewriter COMMAND test_rewriter)
