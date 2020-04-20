include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

#add_definitions("-DUNIT_TESTING")

add_executable(test_clone test_clone.c)
target_link_libraries(test_clone phyc)

add_executable(test_distributions test_distributions.c)
target_link_libraries(test_distributions phyc)

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
	COMMAND "${VALGRIND}" --tool=memcheck --leak-check=yes --error-exitcode=1 --show-possibly-lost=no --show-reachable=no --num-callers=20 --log-file=minimalraw.log --track-fds=yes $<TARGET_FILE:test_clone>)
endif()



message(${JSON_FILES})
  
