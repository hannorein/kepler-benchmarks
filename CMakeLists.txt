cmake_minimum_required(VERSION 3.5)

project(KeplerBenchmarks VERSION 0.1.0 LANGUAGES C)
include(GNUInstallDirs)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message(STATUS "Build-type: ${CMAKE_BUILD_TYPE}")

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Check to see if we need to link libm
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

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
# add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/benchmarks)

# Tests
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tests)
endif()