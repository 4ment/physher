cmake_minimum_required(VERSION 2.6)

SET(DEFAULT_SSE_LEVEL "sse3")

PROJECT(phyc C)

SET( LIBRARY_NAME "PHYC" )

SET( ${LIBRARY_NAME}_MAJOR_VERSION 1 )
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


OPTION(USE_SSE_SUPPORT "Use SSE if available." ON)
OPTION(USE_OPENMP_SUPPORT "Use OpenMP if available." OFF)
OPTION(USE_PTHREAD_SUPPORT "Use PThread if available." ON)

OPTION(LIBRARY_ONLY "Compile library only." OFF)


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
	
	IF(AUTO_SSE_SUPPORT)
		SET(SSE_C_TESTS "sse4_2" "sse4_1" "sse4" "sse3" "sse2" "sse")
		
		message(STATUS "Checking SSE instructions support by current CPU")
		foreach(sse_test ${SSE_C_TESTS})
			#IF(NOT SSE_FOUND)
				IF(WIN32)
					SET(SSE_CHECK_COMMAND "FAILURE")
				elseif(APPLE)
					SET(SSE_CHECK_COMMAND "sysctl -a | grep ${sse_test}")
				else()#other os
					SET(SSE_CHECK_COMMAND "grep ${sse_test} /proc/cpuinfo")
				ENDIF(WIN32)
				
				#execute_process(COMMAND ${SSE_CHECK_COMMAND} RESULT_VARIABLE ret_var OUTPUT_VARIABLE out_var)
				execute_process(COMMAND sh -c ${SSE_CHECK_COMMAND} OUTPUT_VARIABLE out_var)
	
				IF(NOT out_var STREQUAL "")
					SET(SSE_FOUND TRUE)
					MESSAGE(STATUS "${sse_test} detected and working")
					set(USE_SSE_SET ${sse_test})
					
					STRING(REGEX REPLACE "_" "." USE_SSE_SET ${USE_SSE_SET}) 
					SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${USE_SSE_SET}")
					
					#string(TOUPPER ${USE_SSE_SET} UPPER)
					#SET(UPPER "__${UPPER}__")
					#add_definitions("-D${UPPER}")
					#MESSAGE(STATUS "${UPPER}")
			
				ENDIF(NOT out_var STREQUAL "")
	
			#ENDIF(NOT SSE_FOUND)
		endforeach(sse_test)
		
		IF(NOT SSE_FOUND)
			message(STATUS "SSE not supported")
			SET(${LIBRARY_NAME}_SSE_ENABLED, "0")
		ELSE(NOT SSE_FOUND)
			add_definitions("-DSSE3_ENABLED")
			
			STRING(REGEX REPLACE "_" "." USE_SSE_SET ${USE_SSE_SET}) 
			SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${USE_SSE_SET}")
			SET(${LIBRARY_NAME}_SSE_ENABLED "1")
			SET(${LIBRARY_NAME}_SSE_LEVEL ${USE_SSE_SET})
		ENDIF(NOT SSE_FOUND)
	ELSE(AUTO_SSE_SUPPORT)		
		IF(NOT SSE_LEVEL)
			SET(SSE_LEVEL ${DEFAULT_SSE_LEVEL})
		ENDIF(NOT SSE_LEVEL)
		
		message(STATUS "Using SSE ${SSE_LEVEL}")

		add_definitions("-DSSE3_ENABLED")
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${SSE_LEVEL}")
		SET(${LIBRARY_NAME}_SSE_ENABLED "1")
		SET(${LIBRARY_NAME}_SSE_LEVEL ${SSE_LEVEL})
	ENDIF(AUTO_SSE_SUPPORT)
ELSE(USE_SSE_SUPPORT)
	SET(${LIBRARY_NAME}_SSE_ENABLED "0")
	SET(${LIBRARY_NAME}_SSE_LEVEL "0")
ENDIF(USE_SSE_SUPPORT)


#####################################################################
# Config header file
#####################################################################

configure_file (
	"${PROJECT_SOURCE_DIR}/src/phyc/PhyCConfig.h.in"
	"${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h"
)
  
#####################################################################

add_definitions("-std=c99")

add_definitions(-DNO_TRANSPOSE)

link_libraries(m)

file(GLOB SRC_FILES src/phyc/*.c)
file(GLOB HEADER_FILES src/phyc/*.h)

SET(HEADER_FILES "${HEADER_FILES};${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h")

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(
    phyc
    ${SRC_FILES}
)

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

IF(NOT LIBRARY_ONLY)	
	add_executable(physher src/physher.c)
	target_link_libraries(physher phyc)

	add_executable(simultron src/simultron.c)
	target_link_libraries(simultron phyc)

	add_executable(bootstrap src/bootstrap.c)
	target_link_libraries(bootstrap phyc)

	add_executable(modelAveraging src/modelAveraging.c)
	target_link_libraries(modelAveraging phyc)
	
	install(TARGETS physher RUNTIME DESTINATION bin)
	install(TARGETS simultron RUNTIME DESTINATION bin)
	install(TARGETS modelAveraging RUNTIME DESTINATION bin)
	install(TARGETS bootstrap RUNTIME DESTINATION bin)
ENDIF()

# Install libs
install(TARGETS phyc DESTINATION lib EXPORT phyc-targets)

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include/phyc)

install(EXPORT phyc-targets DESTINATION lib/phyc)
install(FILES phyc-config.cmake DESTINATION lib/phyc)
