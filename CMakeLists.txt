cmake_minimum_required(VERSION 3.13...3.28)

project(phyc
    VERSION 2.0.0
    DESCRIPTION "Phylogenetics library for estimating evolutionary rates and divergence times"
    LANGUAGES C
)

#####################################################################
# Version (mirrored into PhyCConfig.h)
#####################################################################
set(PHYC_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(PHYC_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(PHYC_PATCH_LEVEL   ${PROJECT_VERSION_PATCH})
set(PHYC_SOVERSION "${PHYC_MAJOR_VERSION}.${PHYC_MINOR_VERSION}")
set(PHYC_VERSION   "${PHYC_SOVERSION}.${PHYC_PATCH_LEVEL}")

#####################################################################
# Build type / options
#####################################################################
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(BUILD_SHARED_LIBS    "Build shared library (OFF for static)"       ON)
option(BUILD_CPP_WRAPPER    "Build the C++ wrapper (libphycpp)"           OFF)
option(BUILD_BENCHMARKING   "Build the benchmarking executables"          OFF)

option(USE_SSE_SUPPORT      "Use SSE if available"                        ON)
option(USE_AVX_SUPPORT      "Use AVX if available (overrides SSE)"        OFF)
option(USE_OPENMP_SUPPORT   "Use OpenMP if available"                     OFF)
option(USE_PTHREAD_SUPPORT  "Use pthreads if available (overrides OpenMP)" ON)

option(LIBRARY_ONLY         "Build the library only (no physher binary)"  OFF)
option(BUILD_LIBRARY        "Build/install libphyc (OFF: physher only)"   ON)
option(DISABLE_GSL          "Build without GSL"                           OFF)
option(ENABLE_ASAN          "Enable AddressSanitizer"                     OFF)

if(NOT BUILD_LIBRARY)
    if(LIBRARY_ONLY)
        message(FATAL_ERROR "LIBRARY_ONLY and -DBUILD_LIBRARY=off are mutually exclusive")
    endif()
    if(BUILD_CPP_WRAPPER)
        message(FATAL_ERROR "BUILD_CPP_WRAPPER needs libphyc; drop -DBUILD_LIBRARY=off")
    endif()
    message(STATUS "Building the physher executable only (no libphyc)")
endif()

set(DEFAULT_SSE_LEVEL "sse3")

# Accumulated, then applied to the phyc target once it exists.
set(PHYC_COMPILE_DEFINITIONS "")
set(PHYC_COMPILE_OPTIONS "")
set(PHYC_LINK_LIBRARIES "")

# PhyCConfig.h status flags (defaults; flipped on as features are enabled).
set(PHYC_OPENMP_ENABLED  0)
set(PHYC_PTHREAD_ENABLED 0)
set(PHYC_SSE_ENABLED     0)
set(PHYC_SSE_LEVEL       0)
set(PHYC_AVX_ENABLED     0)

#####################################################################
# Compiler-specific tuning
#####################################################################
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_C_COMPILER_VERSION VERSION_GREATER 4.7)
        list(APPEND PHYC_COMPILE_DEFINITIONS GCC_GT_47)
        set(OPENMP_TASK_ENABLED ON)
        message(STATUS "Using GCC >= 4.7")
    endif()
    if(CMAKE_C_COMPILER_VERSION VERSION_LESS 4.4)
        message(STATUS "Disabling OpenMP (needs GCC >= 4.4)")
        set(USE_OPENMP_SUPPORT OFF)
    endif()
endif()

#####################################################################
# Threading: pthreads (preferred) then OpenMP
#####################################################################
if(USE_PTHREAD_SUPPORT)
    find_package(Threads)
    if(Threads_FOUND)
        list(APPEND PHYC_COMPILE_DEFINITIONS PTHREAD_ENABLED)
        list(APPEND PHYC_LINK_LIBRARIES Threads::Threads)
        set(USE_OPENMP_SUPPORT OFF)
        set(PHYC_PTHREAD_ENABLED 1)
        message(STATUS "Using pthreads")
    endif()
endif()

if(USE_OPENMP_SUPPORT)
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        list(APPEND PHYC_LINK_LIBRARIES OpenMP::OpenMP_C)
        set(PHYC_OPENMP_ENABLED 1)
        if(OPENMP_TASK_ENABLED)
            list(APPEND PHYC_COMPILE_DEFINITIONS OPENMP_TASK_ENABLED)
            message(STATUS "Using OpenMP tasks")
        endif()
    else()
        message(STATUS "OpenMP not available; try -DUSE_PTHREAD_SUPPORT=ON")
    endif()
endif()

#####################################################################
# SIMD: AVX (overrides SSE) or SSE
#####################################################################
if(USE_AVX_SUPPORT)
    message(STATUS "Using AVX")
    list(APPEND PHYC_COMPILE_DEFINITIONS AVX_ENABLED)
    list(APPEND PHYC_COMPILE_OPTIONS -mavx)
    set(PHYC_AVX_ENABLED 1)
    set(USE_SSE_SUPPORT OFF)
endif()

if(USE_SSE_SUPPORT)
    if(NOT SSE_LEVEL)
        set(SSE_LEVEL ${DEFAULT_SSE_LEVEL})
    endif()
    list(APPEND PHYC_COMPILE_DEFINITIONS SSE3_ENABLED)
    if(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        message(STATUS "Using sse2neon for arm64 Apple Silicon")
        list(APPEND PHYC_COMPILE_OPTIONS -march=armv8-a+fp+simd+crypto+crc)
    else()
        message(STATUS "Using SSE ${SSE_LEVEL}")
        list(APPEND PHYC_COMPILE_OPTIONS -m${SSE_LEVEL})
    endif()
    set(PHYC_SSE_ENABLED 1)
    set(PHYC_SSE_LEVEL ${SSE_LEVEL})
endif()

#####################################################################
# GSL
#####################################################################
if(NOT DISABLE_GSL)
    find_package(GSL REQUIRED)
    list(APPEND PHYC_LINK_LIBRARIES GSL::gsl GSL::gslcblas)
else()
    list(APPEND PHYC_COMPILE_DEFINITIONS GSL_DISABLED)
    if(NOT BUILD_LIBRARY)
        message(FATAL_ERROR "DISABLE_GSL builds the library only; drop -DBUILD_LIBRARY=off")
    endif()
    set(LIBRARY_ONLY ON)
endif()

#####################################################################
# Config header (git provenance + feature flags)
#####################################################################
execute_process(
    COMMAND git tag --points-at HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HEAD_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

configure_file(
    "${PROJECT_SOURCE_DIR}/src/phyc/PhyCConfig.h.in"
    "${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h"
)

#####################################################################
# AddressSanitizer (applies to every target, including tests)
#####################################################################
if(ENABLE_ASAN AND CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
	message(STATUS "Using AddressSanitizer")
endif()

#####################################################################
# Sources
#####################################################################
if(DISABLE_GSL)
    set(SRC_FILES
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
        src/phyc/model.c
        src/phyc/mt19937ar.c # no .h
        src/phyc/nj.c
        src/phyc/node.c
        src/phyc/nonstat.c
        src/phyc/nucsubst.c
        src/phyc/parameters.c
        src/phyc/random.c
        src/phyc/sequence.c
        src/phyc/sequenceio.c
        src/phyc/sitemodel.c
        src/phyc/sitepattern.c
        src/phyc/solve.c
        src/phyc/statistics.c
        src/phyc/substmodel.c
        src/phyc/transforms.c
        src/phyc/tree.c
        src/phyc/treeio.c
        src/phyc/treelikelihood.c
        src/phyc/treelikelihood4.c
        src/phyc/treelikelihood4CAT.c
        src/phyc/treelikelihood20.c
        src/phyc/treelikelihoodCodon.c
        src/phyc/treelikelihoodX.c
        src/phyc/treetransform.c
        src/phyc/unrest.c
        src/phyc/upgma.c
        src/phyc/utils.c
        src/phyc/wag.c
    )
    set(HEADER_FILES "")
    foreach(SRC_FILE ${SRC_FILES})
        string(REGEX REPLACE "[.]c$" ".h" HDR ${SRC_FILE})
        list(APPEND HEADER_FILES ${HDR})
    endforeach()
    list(REMOVE_ITEM HEADER_FILES src/phyc/mt19937ar.h)
    list(APPEND HEADER_FILES
        src/phyc/gradient.h
        src/phyc/mathconstant.h
        src/phyc/substmodels.h)
else()
    file(GLOB SRC_FILES    CONFIGURE_DEPENDS src/phyc/*.c)
    file(GLOB HEADER_FILES CONFIGURE_DEPENDS src/phyc/*.h)

    list(REMOVE_ITEM SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/phyc/vb.c
        ${PROJECT_SOURCE_DIR}/src/phyc/opvb.c
        ${PROJECT_SOURCE_DIR}/src/phyc/klqp.c)
    list(REMOVE_ITEM HEADER_FILES
        ${PROJECT_SOURCE_DIR}/src/phyc/vb.h
        ${PROJECT_SOURCE_DIR}/src/phyc/opvb.h
        ${PROJECT_SOURCE_DIR}/src/phyc/klqp.h)
endif()

# Generated / vendored headers shipped alongside the library.
file(COPY ${PROJECT_SOURCE_DIR}/contrib/neon2sse.h
    DESTINATION ${PROJECT_BINARY_DIR}/phyc)
list(APPEND HEADER_FILES
    ${PROJECT_BINARY_DIR}/phyc/neon2sse.h
    ${PROJECT_BINARY_DIR}/phyc/PhyCConfig.h)

#####################################################################
# Library target
#####################################################################
# Without BUILD_LIBRARY the sources still go through a target (so the flags,
# includes and test wiring below are shared), but as an OBJECT library: no
# .a/.so is emitted, the objects are linked straight into each executable.
if(BUILD_LIBRARY)
    add_library(phyc ${SRC_FILES})
else()
    add_library(phyc OBJECT ${SRC_FILES})
endif()
add_library(phyc::phyc ALIAS phyc)

set_target_properties(phyc PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    CLEAN_DIRECT_OUTPUT 1
)
if(BUILD_LIBRARY)
    set_target_properties(phyc PROPERTIES
        VERSION ${PHYC_VERSION}
        SOVERSION ${PHYC_SOVERSION}
    )
endif()

target_include_directories(phyc
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/contrib
)

# Feature flags are referenced by public headers (treelikelihood.h, ...) so
# consumers need them; platform macros are internal (utils.c only).
target_compile_definitions(phyc
    PUBLIC
        ${PHYC_COMPILE_DEFINITIONS}
    PRIVATE
        $<$<PLATFORM_ID:Darwin>:__MACH__>
        $<$<PLATFORM_ID:Linux>:__LINUX__;_XOPEN_SOURCE=600>
        $<$<PLATFORM_ID:Windows>:__WIN32__>
)
target_compile_options(phyc PUBLIC ${PHYC_COMPILE_OPTIONS})
target_link_libraries(phyc PUBLIC ${PHYC_LINK_LIBRARIES} m)

#####################################################################
# physher executable
#####################################################################
if(NOT LIBRARY_ONLY)
    add_executable(physher src/physher.c)
    set_target_properties(physher PROPERTIES C_STANDARD 99)
    target_link_libraries(physher PRIVATE phyc)
    install(TARGETS physher RUNTIME DESTINATION bin)
endif()

#####################################################################
# C++ wrapper
#####################################################################
if(BUILD_CPP_WRAPPER)
    enable_language(CXX)
    add_library(phycpp src/phycpp/physher.cpp)
    target_link_libraries(phycpp PUBLIC phyc)
    target_include_directories(phycpp PRIVATE "${PROJECT_SOURCE_DIR}/src")
    set_target_properties(phycpp PROPERTIES
        LINKER_LANGUAGE CXX
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        VERSION ${PHYC_VERSION}
        SOVERSION ${PHYC_SOVERSION}
        CLEAN_DIRECT_OUTPUT 1
    )
    install(TARGETS phycpp DESTINATION lib EXPORT phycpp-targets)
    install(FILES "src/phycpp/physher.hpp" DESTINATION include/phycpp)
    install(EXPORT phycpp-targets DESTINATION lib/phycpp)
endif()

#####################################################################
# Install
#####################################################################
if(BUILD_LIBRARY)
    install(TARGETS phyc DESTINATION lib EXPORT phyc-targets)
    install(FILES ${HEADER_FILES} DESTINATION include/phyc)
    install(EXPORT phyc-targets DESTINATION lib/phyc)
    install(FILES phyc-config.cmake DESTINATION lib/phyc)
endif()

add_custom_target(uninstall COMMAND xargs rm -vf < install_manifest.txt)

#####################################################################
# Tests / benchmarks
#####################################################################
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_BENCHMARKING)
    add_subdirectory(examples)
endif()
