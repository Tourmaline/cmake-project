# HELP: 
# CMake Useful Variables: http://www.cmake.org/Wiki/CMake_Useful_Variables 
# Example: https://github.com/dealii/dealii

# stop if cmake version below 2.8
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(PROJECT ergo)

# project name
project(${PROJECT})

# project version
set(PROJECT_VERSION 3.51)

# set variables for config.h
set(ERGO_SPREFIX     	\"${PROJECT_SOURCE_DIR}\")
set(PACKAGE_VERSION  	\"${PROJECT_VERSION}\")
set(VERSION 		${PACKAGE_VERSION})
set(PACKAGE 		\"${PROJECT_NAME}\")
set(NAME 		ERGO)
set(PACKAGE_NAME 	\"${NAME}\")
set(PACKAGE_BUGREPORT 	\"info@ergoscf.org\")
set(PACKAGE_STRING 	"\"${NAME} ${PROJECT_VERSION}\"")
set(PACKAGE_TARNAME 	${PACKAGE})
set(PACKAGE_URL 	\"\")


# enable c and c++ languages
enable_language(C CXX)


# location of additional cmake modules
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

# guard against in-source builds and bad build-type strings
include(ConfigSafeGuards)
# testing
include(ConfigTesting)

# detect operating system 
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
# detect host processor
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

#-----------------------------
# CHECK OPTIONS AND LIBRARIES
#-----------------------------

option(ENABLE_OMP 		"Enable OpenMP" OFF)
option(USE_LINALG_TEMPLATES 	"Use C++ Linear Algebra Template Library (BLAS+LAPACK)" OFF)
option(USE_SSE_INTRINSICS 	"Use SSE intrinsics to generate efficient computational kernels" OFF)
option(PRECISION_SINGLE  	"Use single precision" OFF)
option(PRECISION_DOUBLE		"Use single precision" ON)
option(PRECISION_LONG_DOUBLE 	"Use single precision" OFF)
option(STDC_HEADERS 		"ANSI C header files are used" ON)
option(USE_BLAS_IN_XC 		"Use BLAS for small matrices in XC integration" OFF)
option(YYTEXT_POINTER  		"`lex' declares `yytext' as a `char *' by default, not a `char[]'" ON)




if(PRECISION_SINGLE OR PRECISION_LONG_DOUBLE)
  set(PRECISION_DOUBLE OFF)
endif()


if(PRECISION_SINGLE)
  set(FUNC_PRECISION 1)
endif()
if(PRECISION_DOUBLE)
  set(FUNC_PRECISION 2)
endif()
if(PRECISION_LONG_DOUBLE)
  set(FUNC_PRECISION 3)
endif()



# CHECK FOR PTHREADS
find_package(Threads REQUIRED)
set(THREAD_LIBRARIES pthread)

# CHECK FOR MATH FUNCTIONS

#set(CMAKE_REQUIRED_INCLUDES math.h)
#set(CMAKE_REQUIRED_LIBRARIES m)

include(CheckFunctionExists)
check_function_exists(atexit  HAVE_ATEXIT)

check_function_exists ( vprintf HAVE_VPRINTF )
if ( NOT HAVE_VPRINTF )
check_function_exists ( _doprnt HAVE_DOPRNT )
endif ( )

include(CheckLibraryExists)
check_library_exists(m sqrtl  "" HAVE_SQRTL)
check_library_exists(m erf  "" HAVE_ERF)
check_library_exists(m erfc  "" HAVE_ERFC)
check_library_exists(m erfcf  "" HAVE_ERFCF)
check_library_exists(m erfcl "" HAVE_ERFCL)
check_library_exists(m erff "" HAVE_ERFF)
check_library_exists(m erfl "" HAVE_ERFL)

include(CheckIncludeFiles)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stddef.h HAVE_STDDEF_H)


# CHECK FOR OPENMP
if(ENABLE_OMP)
  add_definitions(-DHAVE_OPENMP) # used in the code
endif()


# CHECK FOR Linear Algebra Template Library 

 # if user specify to use internal library
 if(USE_LINALG_TEMPLATES)
   if(BLAS_LIBRARIES OR BLA_VENDOR) # if user also ask for the external BLAS library, write error
     message(FATAL_ERROR "option USE_LINALG_TEMPLATES conflicts with option BLAS_LIBRARIES
	Try to remove cache file: rm CMakeCache.txt")    
   else()
     # use internal BLAS library
     # add_definitions(-DUSE_LINALG_TEMPLATES=ON)
     message(STATUS "Use internal C++ Linear Algebra Template Library (BLAS+LAPACK)") 
     set(LAPACK_LIBRARIES templateblas templatelapack)
  endif()
 endif()
 
 # if user do not specify to use internal library
 if(NOT USE_LINALG_TEMPLATES)
   # needed for FindLAPACK module ?
   #enable_language(Fortran)
 
   # module check for existence of BLAS library and compatibility of LAPACK and BLAS libraries
   # (if it can link them together)
   find_package(LAPACK)
 
   if(LAPACK_FOUND) 
     message (STATUS "Found working BLAS/LAPACK libraries: ${LAPACK_LIBRARIES}.")
     # set linker options
     set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}" 
       CACHE STRING "Linker flags for executables" FORCE )
     set ( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}" 
       CACHE STRING "Linker flags for modules" FORCE )
     set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}" 
       CACHE STRING "Linker flags for shared libs" FORCE )
    
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} templateblas templatelapack)
   else() # BLAS or LAPACK is not found
     message (FATAL_ERROR "
	Cannot find working BLAS/LAPACK libraries.
     	Specify BLAS vendor (cmake -DBLA_VENDOR=$vendor$ /path/to/source)
     	where $vendor$ can be for example Goto, ATLAS, Generic, Intel,
     	or specify path to the library (cmake -DBLAS_LIBRARIES=/path/to/library /path/to/source).
     	Consider checking USE_LINALG_TEMPLATES (cmake -DUSE_LINALG_TEMPLATES=ON /path/to/source)." )
   endif()
 endif()

# set include directories for internal library 
set(INTERNAL_BLAS_INCLUDE_DIRS   "${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/blas")
set(INTERNAL_LAPACK_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/lapack")




# CHECK FOR SSE
if(USE_SSE_INTRINSICS)
FIND_PACKAGE(SSE)
IF(C_SSE2_FOUND)
  SET(CMAKE_C_FLAGS "${C_SSE2_FLAGS} -DUSE_SSE2 ${CMAKE_C_FLAGS}")
ENDIF(C_SSE2_FOUND)
IF(C_SSE3_FOUND)
  SET(CMAKE_C_FLAGS "${C_SSE3_FLAGS} -DUSE_SSE3 ${CMAKE_C_FLAGS}")
ENDIF(C_SSE3_FOUND)
IF(C_SSE4_1_FOUND)
  SET(CMAKE_C_FLAGS "${C_SSE4_1_FLAGS} -DUSE_SSE4_1 ${CMAKE_C_FLAGS}")
ENDIF(C_SSE4_1_FOUND)
IF(C_SSE4_2_FOUND)
  SET(CMAKE_C_FLAGS "${C_SSE4_2_FLAGS} -DUSE_SSE4_2 ${CMAKE_C_FLAGS}")
ENDIF(C_SSE4_2_FOUND)
endif()






# CREATE CONFIGURE FILE 
configure_file (
    ${PROJECT_SOURCE_DIR}/cmake_config.h.in
    ${PROJECT_SOURCE_DIR}/source/config.h
    @ONLY # this means process only @VARIABLES@ and not ${variables}
    )



#-------------------------------
# COMPILATION AND INSTALLATION  
#-------------------------------

# we default to Release build type
# can be None, Release, Debug, RelWithDebInfo, MinSizeRel  
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS           "${CMAKE_C_FLAGS} -Wall ")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall ")
    set(CMAKE_C_FLAGS_RELEASE   "-O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG     "-O0 -g3 -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -fprofile-arcs -ftest-coverage")
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG    "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE  "-O3")
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif()
endif()





# setup CPack
# cpack -D CPACK_GENERATOR="ZIP;TGZ" /path/to/build/tree
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/ergo_license_long.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION}")
set(CPACK_PACKAGE_CONTACT "ergoscf.org")
include(CPack)


# add subdirectory
add_subdirectory(source)
file(COPY test DESTINATION ${PROJECT_BINARY_DIR})

