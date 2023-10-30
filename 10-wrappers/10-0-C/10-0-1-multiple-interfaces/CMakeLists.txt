cmake_minimum_required(VERSION 3.9)
project(multi-interfaces-C VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-0-1" LANGUAGES C)

set(MUI_BASE_C_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_C_LIB_DIR "" CACHE PATH "Path to MUI C library directory")
set(MUI_C_INCLUDE_DIR "" CACHE PATH "Path to MUI C include directory")

set(MUI_C_CPP_SEPERATE_PATH 0)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MUI QUIET)

if(MUI_BASE_DIR)
  set(MUI_C_LIB_DIR "${MUI_BASE_DIR}/lib")
  set(MUI_C_INCLUDE_DIR "${MUI_BASE_DIR}/include/wrappers/C")
endif()

find_library(MUI_C_WRAPPER_LIB NAMES MUI_C_wrapper PATHS ${MUI_C_LIB_DIR})

if(NOT MUI_C_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_C_wrapper.a in ${MUI_C_LIB_DIR}")
endif()

include_directories(${MUI_C_INCLUDE_DIR})

add_executable (multi-interfaces multi-interfaces.c)
target_link_libraries(multi-interfaces ${MUI_C_WRAPPER_LIB})
