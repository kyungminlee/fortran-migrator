# Top-level test directory — added by cmake/CMakeLists.txt when
# BUILD_TESTING=ON. Each library has its own subdir; we only descend
# into ones whose precision library was actually built (i.e., the
# library was migrated for this target) AND whose CMakeLists.txt
# exists in the staging copy.

foreach(_lib ${STAGED_LIBRARIES})
    set(_libdir ${CMAKE_CURRENT_SOURCE_DIR}/${_lib})
    if(TARGET ${LIB_PREFIX}${_lib}
            AND EXISTS ${_libdir}/CMakeLists.txt)
        add_subdirectory(${_lib})
    endif()
endforeach()
