
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds are not allowed. Please make a new directory (called a build directory) and run CMake from there. In case of error, try to delete generated cmake files from current directory.")
endif()


