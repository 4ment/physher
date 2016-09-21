get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/phyc-targets.cmake)
get_filename_component(PHYC_INCLUDE_DIRS "${SELF_DIR}/../../include/phyc" ABSOLUTE)