cmake_minimum_required(VERSION 2.8)
project(pq-fast-scan)

include(FeatureSummary)

message("PQ Fast Scan can only be built on hosts "
    "supporting SSE2, SSE3, SSSE3 and POPCNT")

option(AVX "Compile scan functions requiring AVX intrinsics" ON)
option(AVX2 "Compile scan functions requiring AVX2 intrinsics" OFF)

# Options
add_feature_info(AVX AVX
"Compile scan functions requiring AVX intrinsics")
add_feature_info(AVX2 AVX2
"Compile scan functions requiring AVX2 intrinsics (automatically enables AVX)")
if(AVX2)
    set(AVX ON)
endif()

feature_summary(WHAT ALL)

# Set GCC Flags
set(GCC_FLAGS
"-g -ggdb -O3 -m64 -march=native -ffast-math -Wall -Werror")

execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                OUTPUT_VARIABLE GCC_VERSION)

if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
    set(GCC_FLAGS "-fdiagnostics-color=auto ${GCC_FLAGS}")
endif()

message("GCC Flags: ${GCC_FLAGS}")

# C++ Flags
set(CMAKE_CXX_FLAGS "-std=c++11 ${GCC_FLAGS}")

# Configure options
configure_file("${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h")
include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}")

# Libraries
set(LINK_LIBS pfm)

# Object files
add_library(common OBJECT common.cpp)
add_library(binheap OBJECT binheap.cpp)
add_library(benchmark OBJECT benchmark.cpp perfevents.cpp)
add_library(scan_naive OBJECT pqscan/scan_naive.cpp)

# Target: pqscan - PQ Scan implementations - Synthetic data
set(PQSCAN_SRCS pqscan/scan_sse.cpp)
if(AVX)
    list(APPEND PQSCAN_SRCS pqscan/scan_avx.cpp)
endif()

if(AVX2)
    list(APPEND PQSCAN_SRCS pqscan/scan_gather.cpp)
endif()

add_executable(pqscan
    pqscan/pqscan.cpp
    pqscan/populate.cpp
    ${PQSCAN_SRCS}
    $<TARGET_OBJECTS:scan_naive>
    $<TARGET_OBJECTS:common>
    $<TARGET_OBJECTS:binheap>
    $<TARGET_OBJECTS:benchmark>)
target_link_libraries(pqscan ${LINK_LIBS})

# Target: pqfastcan - PQ Fast Scan - Real data
add_executable(pqfastscan
    pqfastscan/pqfastscan.cpp
    pqfastscan/fastscan.cpp
    pqfastscan/layout.cpp
    $<TARGET_OBJECTS:scan_naive>
    $<TARGET_OBJECTS:common>
    $<TARGET_OBJECTS:binheap>
    $<TARGET_OBJECTS:benchmark>)
target_link_libraries(pqfastscan ${LINK_LIBS})
