cmake_minimum_required(VERSION 3.9)
project(fixedRelaxation_basic-C VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-0-4" LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(MUI_BASE_C_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_C_LIB_DIR "" CACHE PATH "Path to MUI C library directory")
set(MUI_C_INCLUDE_DIR "" CACHE PATH "Path to MUI C include directory")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MUI REQUIRED)

if(MUI_BASE_DIR)
  set(MUI_C_LIB_DIR "${MUI_BASE_DIR}/lib")
  set(MUI_C_INCLUDE_DIR "${MUI_BASE_DIR}/include/wrappers/C")
  include_directories(${MUI_INCLUDE_DIR})
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

find_library(MUI_C_WRAPPER_LIB NAMES MUI_C_wrapper PATHS ${MUI_C_LIB_DIR})

if(NOT MUI_C_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_C_wrapper.a in ${MUI_C_LIB_DIR}")
endif()

include_directories(${MUI_C_INCLUDE_DIR})

add_executable (heat-left heat-left.cpp)
target_link_libraries (heat-left ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

add_executable (heat-right heat-right.c)
target_link_libraries(heat-right ${MUI_C_WRAPPER_LIB} ${MPI_LIBRARIES})