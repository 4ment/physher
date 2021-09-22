include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

add_executable(benchmarking benchmarking.c)
target_link_libraries(benchmarking phyc)
