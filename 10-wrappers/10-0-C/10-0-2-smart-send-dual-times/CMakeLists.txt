cmake_minimum_required(VERSION 3.9)
project(smart-send-dual-times-C VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-0-2" LANGUAGES C)

set(MUI_BASE_C_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_C_LIB_DIR "" CACHE PATH "Path to MUI C library directory")
set(MUI_C_INCLUDE_DIR "" CACHE PATH "Path to MUI C include directory")

set(MUI_C_CPP_SEPERATE_PATH 0)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

find_package(MUI QUIET)

if(MUI_BASE_DIR)
  set(MUI_C_LIB_DIR "${MUI_LIB_DIR}")
  set(MUI_C_INCLUDE_DIR "${MUI_INCLUDE_DIR}/wrappers/C")
elseif(MUI_BASE_C_DIR)
  set(MUI_C_LIB_DIR "${MUI_BASE_C_DIR}/lib")
  set(MUI_C_INCLUDE_DIR "${MUI_BASE_C_DIR}/include/wrappers/C")
endif()

find_library(MUI_C_WRAPPER_LIB NAMES MUI_C_wrapper PATHS ${MUI_C_LIB_DIR})

if(NOT MUI_C_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_C_wrapper.a in ${MUI_C_LIB_DIR}")
endif()

include_directories(${MUI_C_INCLUDE_DIR})

add_executable (code1 code1.c)
target_link_libraries(code1 ${MUI_C_WRAPPER_LIB} ${MPI_LIBRARIES})

add_executable (code2 code2.c)
target_link_libraries(code2 ${MUI_C_WRAPPER_LIB} ${MPI_LIBRARIES})