file(GLOB_RECURSE all_sources RELATIVE ${CMAKE_CURRENT_LIST_DIR} *.c)

include(CheckCSourceCompiles)
set(LIBM_TEST_SOURCE "#include<math.h>\nfloat f; int main(){sqrt(f);return 0;}")
check_c_source_compiles("${LIBM_TEST_SOURCE}" HAVE_MATH)
if(HAVE_MATH)
  set(LIBM_LIBRARIES)
else()
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_c_source_compiles("${LIBM_TEST_SOURCE}" HAVE_LIBM_MATH)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT HAVE_LIBM_MATH)
    message(FATAL_ERROR "Unable to use C math library functions")
  endif()
  set(LIBM_LIBRARIES m)
endif()

add_custom_target(results_directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/results)

foreach(source ${all_sources})
  get_filename_component(source_name ${source} NAME_WE)
  add_executable(${source_name} ${source})
  if(WIN32)
    target_compile_options(${source_name} PUBLIC /Wall /bigobj)
  else()
    target_compile_options(${source_name} PUBLIC -Wall -pedantic -Wextra -Werror -ffast-math -march=native)
  endif()
  target_link_libraries(${source_name} PUBLIC ${LIBM_LIBRARIES})

  add_custom_target(
    ${source_name}_output ALL
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${source_name}
    DEPENDS ${source_name} results_directory
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/results/${source_name}.csv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results
  )

  message(STATUS "Adding benchmark: ${source_name}")
endforeach()