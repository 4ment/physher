include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

#add_definitions("-DUNIT_TESTING")

# add_executable(test_clone test_clone.c)
# target_link_libraries(test_clone phyc)

add_executable(test_json test_json.c)
target_link_libraries(test_json phyc)

add_executable(test_distributions test_distributions.c)
target_link_libraries(test_distributions phyc)

add_executable(test_coalescent test_coalescent.c)
target_link_libraries(test_coalescent phyc)

add_executable(test_priors test_priors.c)
target_link_libraries(test_priors phyc)

add_executable(test_transforms test_transforms.c)
target_link_libraries(test_transforms phyc)

add_executable(test_tree test_tree.c)
target_link_libraries(test_tree phyc)

add_executable(test_tree_likelihood test_tree_likelihood.c)
target_link_libraries(test_tree_likelihood phyc)

add_executable(test_compound test_compound.c)
target_link_libraries(test_compound phyc)

add_executable(test_operator test_operator.c)
target_link_libraries(test_operator phyc)

add_executable(test_mcmc test_mcmc.c)
target_link_libraries(test_mcmc phyc)

add_executable(test_parsimony test_parsimony.c)
target_link_libraries(test_parsimony phyc)

add_executable(test_bootstrap test_bootstrap.c)
target_link_libraries(test_bootstrap phyc)


# Refresh test fixtures on every build (not just at configure time) so that
# added/edited files under data/ are picked up by an existing build tree.
add_custom_target(test_data ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
  COMMENT "Copying test data fixtures")
  
#file(GLOB JSON_FILES ${CMAKE_CURRENT_BINARY_DIR}/data/*.json)
#EXECUTE_PROCESS(
#    	COMMAND sed -i "s@tiny.fa@data/tiny.fa@" ${JSON_FILES}
#)

# add_test (NAME Test_clone
# 		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
# 		COMMAND test_clone)

add_test (NAME Test_json
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_json)

add_test (NAME Test_distributions
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_distributions)

add_test (NAME Test_coalescent
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_coalescent)

add_test (NAME Test_priors
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_priors)

add_test (NAME test_transforms
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_transforms)

add_test (NAME Test_tree
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_tree)

add_test (NAME Test_tree_likelihood
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_tree_likelihood)

add_test (NAME Test_compound
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_compound)

add_test (NAME Test_operator
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_operator)

add_test (NAME Test_mcmc
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_mcmc)

add_test (NAME Test_parsimony
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_parsimony)

add_test (NAME Test_bootstrap
		WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data"
		COMMAND test_bootstrap)


find_program(VALGRIND "valgrind")
set( CTEST_MEMORYCHECK_COMMAND valgrind )
set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=callgrind -v" )

if(VALGRIND)
message(${VALGRIND})
	add_custom_target(valgrind
	COMMAND "${VALGRIND}" --tool=memcheck --leak-check=yes --error-exitcode=1 --show-possibly-lost=no --show-reachable=no --num-callers=20 --log-file=minimalraw.log --track-fds=yes $<TARGET_FILE:test_distributions>)
endif()
