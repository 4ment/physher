cmake_minimum_required(VERSION 3.5)

SET(DEFAULT_SSE_LEVEL "sse3")

PROJECT(phyc C CXX)

SET( LIBRARY_NAME "PHYC" )

SET( ${LIBRARY_NAME}_MAJOR_VERSION 2 )
SET( ${LIBRARY_NAME}_MINOR_VERSION 0 )
SET( ${LIBRARY_NAME}_PATCH_LEVEL 0 )


#####################################################################
# Settings
#####################################################################
# library *nix style versioning
SET( ${LIBRARY_NAME}_SOVERSION "${${LIBRARY_NAME}_MAJOR_VERSION}.${${LIBRARY_NAME}_MINOR_VERSION}" )
SET( ${LIBRARY_NAME}_VERSION "${${LIBRARY_NAME}_SOVERSION}.${${LIBRARY_NAME}_PATCH_LEVEL}" )

IF(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type ${CMAKE_BUILD_TYPE}")

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static library (default shared)" ON )

OPTION( BUILD_CPP_WRAPPER "Build C++ wrapper (default off)" OFF )

OPTION(USE_SSE_SUPPORT "Use SSE if available." ON)
OPTION(USE_OPENMP_SUPPORT "Use OpenMP if available." OFF)
OPTION(USE_PTHREAD_SUPPORT "Use PThread if available." ON)

OPTION(LIBRARY_ONLY "Compile library only." OFF)
OPTION(DISABLE_GSL "Compile without GSL." OFF)


#####################################################################
# 32/64 bit
#####################################################################
IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	MESSAGE( STATUS "32 bit architecture detected" )
ELSEIF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
	MESSAGE( STATUS "64 bit architecture detected" )

	IF( BUILD_32BIT_COMPATIBLE )
    	# only works with cmake >= 2.6
		SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
		SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32" )
		MESSAGE( STATUS "Build in 32 bit compatibility mode" )
	ELSE()
		MESSAGE( STATUS "Build in native 64 bit mode" )
	ENDIF()
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)

#####################################################################
# GCC
#####################################################################
IF(CMAKE_COMPILER_IS_GNUCC)
    EXECUTE_PROCESS(
        COMMAND ${CMAKE_C_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
    )
    
    # There are extra optimization capabilities with these compilers
    IF(${GCC_VERSION} VERSION_GREATER 4.7)
    	add_definitions("-DGCC_GT_47")
    	MESSAGE("-- Using GGC >= 4.7")
        SET(OPENMP_TASK_ENABLED "1")
    ENDIF(${GCC_VERSION} VERSION_GREATER 4.7)
    
    IF(${GCC_VERSION} VERSION_LESS 4.4)
        MESSAGE(STATUS "Disabling OpenMP because GCC >= 4.4 is needed.")
        SET(USE_OPENMP_SUPPORT OFF)
    ENDIF(${GCC_VERSION} VERSION_LESS 4.4)
ENDIF(CMAKE_COMPILER_IS_GNUCC)


#####################################################################
# PTHREAD
#####################################################################

SET(${LIBRARY_NAME}_PTHREAD_ENABLED "0")

IF(USE_PTHREAD_SUPPORT)

    find_package(Threads)

    IF(THREADS_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
        add_definitions("-DPTHREAD_ENABLED")
        SET(USE_OPENMP_SUPPORT OFF)
        SET(${LIBRARY_NAME}_PTHREAD_ENABLED "1")
    ENDIF(THREADS_FOUND)
ENDIF(USE_PTHREAD_SUPPORT)

#####################################################################
# OPENMP
#####################################################################

SET(${LIBRARY_NAME}_OPENMP_ENABLED "0")

IF(USE_OPENMP_SUPPORT)

    #ENABLE_LANGUAGE(CXX)
    find_package(OpenMP)
            
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(${LIBRARY_NAME}_OPENMP_ENABLED "1")
        IF(OPENMP_TASK_ENABLED)
            add_definitions("-DOPENMP_TASK_ENABLED")
            MESSAGE("-- Using OpenMP TASKS")
        ENDIF()
    ELSE()
    message(STATUS "OPENMP is not available on this machine. Try to use pthreads instead. cmake -DUSE_PTHREAD_SUPPORT ..")
    ENDIF()
ENDIF()


#####################################################################
# AVX
#####################################################################

SET(${LIBRARY_NAME}_AVX_ENABLED "0")

IF(USE_AVX_SUPPORT)
    message(STATUS "Using AVX")

    add_definitions("-DAVX_ENABLED")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
	SET(${LIBRARY_NAME}_AVX_ENABLED "1")
        
    SET(USE_SSE_SUPPORT OFF)
ENDIF(USE_AVX_SUPPORT)

#####################################################################
# SSE
# MacPro: sse4_2, sse4_1, sse3, sse2, sse
# Mac mini server: none
# MacBook Pro: sse4_2, sse4_1, sse3, sse2, sse
#####################################################################

IF(USE_SSE_SUPPORT)
	IF(NOT SSE_LEVEL)
		SET(SSE_LEVEL ${DEFAULT_SSE_LEVEL})
	ENDIF(NOT SSE_LEVEL)

	add_definitions("-DSSE3_ENABLED")

	IF(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
		MESSAGE(STATUS "Using sse2neon for arm64 Apple M1")
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv8-a+fp+simd+crypto+crc")
	ELSE()
		MESSAGE(STATUS "Using SSE ${SSE_LEVEL}")
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${SSE_LEVEL}")
	ENDIF()

	SET(${LIBRARY_NAME}_SSE_ENABLED "1")
	SET(${LIBRARY_NAME}_SSE_LEVEL ${SSE_LEVEL})
ELSE(USE_SSE_SUPPORT)
	SET(${LIBRARY_NAME}_SSE_ENABLED "0")
	SET(${LIBRARY_NAME}_SSE_LEVEL "0")
ENDIF(USE_SSE_SUPPORT)


#####################################################################
# GSL
#####################################################################

IF(NOT DISABLE_GSL)
	find_package(GSL REQUIRED)
	include_directories(BEFORE SYSTEM ${GSL_INCLUDE_DIRS})
	link_directories(${GSL_LIBRARY_DIRS})
ELSE()
	add_definitions("-DGSL_DISABLED")
	SET(LIBRARY_ONLY ON)
ENDIF()

#####################################################################
# Config header file
#####################################################################

execute_process(
  COMMAND git tag --points-at HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_HEAD_TAG
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

configure_file (
	"${PROJECT_SOURCE_DIR}/src/phyc/PhyCConfig.h.in"
	"${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h"
)
  
#####################################################################

# add_compile_definitions("-std=gnu99")

link_libraries(m)

IF(DISABLE_GSL)
	SET(SRC_FILES
		src/phyc/branchmodel.c
		src/phyc/ctmcscale.c
		src/phyc/datatype.c
		src/phyc/dayhoff.c
		src/phyc/demographicmodels.c
		src/phyc/discreteparameter.c
		src/phyc/distancematrix.c
		src/phyc/distmodel.c
		src/phyc/eigen.c
		src/phyc/errorfuntion.c
		src/phyc/exponential.c
		src/phyc/f81.c
		src/phyc/filereader.c
		src/phyc/gamma.c
		src/phyc/gausslaguerre.c
		src/phyc/gaussian.c
		src/phyc/geneticcode.c
		src/phyc/gensubst.c
		src/phyc/gtr.c
		src/phyc/hashfunctions.c
		src/phyc/hashtable.c
		src/phyc/hessenberg.c
		src/phyc/hky.c
		src/phyc/jc69.c
		src/phyc/K80.c
		src/phyc/lg.c
		src/phyc/lognormal.c
		src/phyc/matrix.c
		src/phyc/mjson.c
		src/phyc/mstring.c
		src/phyc/mt19937ar.c # no .h
		src/phyc/nj.c
		src/phyc/node.c
		src/phyc/nonstat.c
		src/phyc/nucsubst.c
		src/phyc/parameters.c
		src/phyc/random.c
		src/phyc/sequence.c
		src/phyc/sequenceio.c
		src/phyc/simplex.c
		src/phyc/sitemodel.c
		src/phyc/sitepattern.c
		src/phyc/solve.c
		src/phyc/statistics.c
		src/phyc/sitemodel.c
		src/phyc/substmodel.c
		src/phyc/transforms.c
		src/phyc/tree.c
		src/phyc/treeio.c
		src/phyc/treelikelihood.c
		src/phyc/treelikelihood4.c
		src/phyc/treelikelihood20.c
		src/phyc/treelikelihoodCodon.c
		src/phyc/treelikelihoodX.c
		src/phyc/treetransform.c
		src/phyc/unrest.c
		src/phyc/upgma.c
		src/phyc/utils.c
		src/phyc/wag.c
	)
	foreach (SRC_FILE ${SRC_FILES})
		string(REGEX REPLACE "[.]c$" ".h" OBJ ${SRC_FILE})
		set(HEADER_FILES ${HEADER_FILES} ${OBJ})
	endforeach ()
	list(REMOVE_ITEM HEADER_FILES src/phyc/mt19937ar.h)
	list(APPEND HEADER_FILES src/phyc/gradient.h src/phyc/mathconstant.h src/phyc/substmodels.h)

ELSE()
	file(GLOB SRC_FILES src/phyc/*.c)
	file(GLOB HEADER_FILES src/phyc/*.h)

	list(REMOVE_ITEM HEADER_FILES ${PROJECT_SOURCE_DIR}/src/phyc/vb.h ${PROJECT_SOURCE_DIR}/src/phyc/opvb.h)
	list(REMOVE_ITEM SRC_FILES ${PROJECT_SOURCE_DIR}/src/phyc/vb.c ${PROJECT_SOURCE_DIR}/src/phyc/opvb.c)
ENDIF(DISABLE_GSL)

include_directories("${PROJECT_SOURCE_DIR}/contrib")
SET(HEADER_FILES "${HEADER_FILES};${PROJECT_BINARY_DIR}/phyc/neon2sse.h")
file(COPY ${PROJECT_SOURCE_DIR}/contrib/neon2sse.h
  DESTINATION ${PROJECT_BINARY_DIR}/phyc)

SET(HEADER_FILES "${HEADER_FILES};${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h")

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(
    phyc
    ${SRC_FILES}
)
set_property(TARGET phyc PROPERTY C_STANDARD 99)

set_target_properties(
	phyc
	PROPERTIES
	# create *nix style library versions + symbolic links
    VERSION ${${LIBRARY_NAME}_VERSION}
    SOVERSION ${${LIBRARY_NAME}_SOVERSION}
    
	CLEAN_DIRECT_OUTPUT 1
)
  
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_property(
        TARGET phyc
        APPEND PROPERTY COMPILE_DEFINITIONS __MACH__
    )
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set_property(
        TARGET phyc
        APPEND PROPERTY COMPILE_DEFINITIONS __LINUX__ _XOPEN_SOURCE=600
    )
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_property(
        TARGET phyc
        APPEND PROPERTY COMPILE_DEFINITIONS __WIN32__
    )
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


ADD_CUSTOM_COMMAND(
		TARGET phyc
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/phyc"
)

IF(NOT DISABLE_GSL)
	target_link_libraries(phyc GSL::gsl GSL::gslcblas)
ENDIF()

IF(NOT NO_LIBRARY)
	add_executable(physher src/physher.c ${SRC_FILES})
	target_link_libraries(physher GSL::gsl GSL::gslcblas)
	install(TARGETS physher RUNTIME DESTINATION bin)
ELSEIF(NOT LIBRARY_ONLY AND NOT DISABLE_GSL)
	add_executable(physher src/physher.c)
	target_link_libraries(physher phyc)

	
	install(TARGETS physher RUNTIME DESTINATION bin)
ENDIF()

# Install libs
install(TARGETS phyc DESTINATION lib EXPORT phyc-targets)

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include/phyc)

install(EXPORT phyc-targets DESTINATION lib/phyc)
install(FILES phyc-config.cmake DESTINATION lib/phyc)

if(BUILD_CPP_WRAPPER)
    add_library(
        phycpp
        "src/phycpp/physher.cpp"
    )
    set_target_properties(phycpp PROPERTIES LINKER_LANGUAGE CXX)
    set_property(TARGET phycpp PROPERTY CXX_STANDARD 17)
    target_link_libraries(phycpp phyc)
    target_include_directories(phycpp PRIVATE "src/")
    set_target_properties(
        phycpp
        PROPERTIES
        # create *nix style library versions + symbolic links
        VERSION ${${LIBRARY_NAME}_VERSION}
        SOVERSION ${${LIBRARY_NAME}_SOVERSION}
        
        CLEAN_DIRECT_OUTPUT 1
    )
	install(TARGETS phycpp DESTINATION lib EXPORT phycpp-targets)
	install(FILES "src/phycpp/physher.hpp" DESTINATION include/phycpp)
	install(EXPORT phycpp-targets DESTINATION lib/phycpp)
ENDIF(BUILD_CPP_WRAPPER)

add_custom_target(uninstall COMMAND xargs rm -vf <install_manifest.txt )

#
# Tests
#

IF (BUILD_TESTING)
	enable_testing()
	add_subdirectory(tests)
ENDIF (BUILD_TESTING)

IF (BUILD_BENCHMARKING)
	add_subdirectory(examples)
ENDIF (BUILD_BENCHMARKING)