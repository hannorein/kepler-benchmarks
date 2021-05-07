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

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/benchmarks)
