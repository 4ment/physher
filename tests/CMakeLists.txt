include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

add_executable(test_clone test_clone.c)
target_link_libraries(test_clone phyc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  
add_test (NAME Test_clone
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_clone)


find_program(VALGRIND "valgrind" REQUIRED)
set( CTEST_MEMORYCHECK_COMMAND valgrind )
set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=callgrind -v" )

if(VALGRIND)
message(${VALGRIND})
	add_custom_target(valgrind
	COMMAND "${VALGRIND}" --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes $<TARGET_FILE:test_clone>)
endif()
