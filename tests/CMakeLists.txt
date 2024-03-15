include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

#add_definitions("-DUNIT_TESTING")

add_executable(test_clone test_clone.c)
target_link_libraries(test_clone phyc)

add_executable(test_distributions test_distributions.c)
target_link_libraries(test_distributions phyc)

add_executable(test_coalescent test_coalescent.c)
target_link_libraries(test_coalescent phyc)

add_executable(test_simplex test_simplex.c)
target_link_libraries(test_simplex phyc)

add_executable(test_tree test_tree.c)
target_link_libraries(test_tree phyc)

add_executable(test_tree_likelihood test_tree_likelihood.c)
target_link_libraries(test_tree_likelihood phyc)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  
#file(GLOB JSON_FILES ${CMAKE_CURRENT_BINARY_DIR}/data/*.json)
#EXECUTE_PROCESS(
#    	COMMAND sed -i "s@tiny.fa@data/tiny.fa@" ${JSON_FILES}
#)

add_test (NAME Test_clone
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_clone)

add_test (NAME Test_distributions
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_distributions)

add_test (NAME Test_coalescent
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_coalescent)

add_test (NAME Test_simplex
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_simplex)
add_test (NAME Test_tree
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_tree)

add_test (NAME Test_tree_likelihood
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_tree_likelihood)


find_program(VALGRIND "valgrind")
set( CTEST_MEMORYCHECK_COMMAND valgrind )
set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=callgrind -v" )

if(VALGRIND)
message(${VALGRIND})
	add_custom_target(valgrind
	COMMAND "${VALGRIND}" --tool=memcheck --leak-check=yes --error-exitcode=1 --show-possibly-lost=no --show-reachable=no --num-callers=20 --log-file=minimalraw.log --track-fds=yes $<TARGET_FILE:test_distributions>)
endif()
