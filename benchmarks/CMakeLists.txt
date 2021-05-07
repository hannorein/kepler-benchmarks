file(GLOB_RECURSE all_sources RELATIVE ${CMAKE_CURRENT_LIST_DIR} *.c)

add_custom_target(results_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/results)

foreach(source ${all_sources})
  get_filename_component(source_name ${source} NAME_WE)
  add_executable(${source_name} ${source})
  if(WIN32)
    target_compile_options(${source_name} PUBLIC /Wall /bigobj)
  else()
    target_compile_options(${source_name} PUBLIC -O3 -Wall -pedantic -Wextra -Werror -ffast-math -march=native)
  endif()
  target_link_libraries(${source_name} PUBLIC ${LIBM_LIBRARIES})

  # This is a hack to force the results directory to be created; it doesn't work when added as a dep to the custom_command
  add_dependencies(${source_name} results_directory)

  # Run the benchmark after it is built
  add_custom_command(
    TARGET ${source_name} POST_BUILD
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${source_name} > ${CMAKE_CURRENT_BINARY_DIR}/results/${source_name}.csv
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/results/${source_name}.csv
  )

  message(STATUS "Adding benchmark: ${source_name}")
endforeach()