cmake_minimum_required(VERSION 3.9)
project(aitken-dynamic-points-Fortran VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-1-5" LANGUAGES Fortran CXX)

# Set the Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_COMPILER mpifort)
set(MUI_BASE_FORTRAN_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_FORTRAN_LIB_DIR "" CACHE PATH "Path to MUI FORTRAN library directory")
set(MUI_FORTRAN_INCLUDE_DIR "" CACHE PATH "Path to MUI FORTRAN include directory")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

find_package(MUI REQUIRED)

if(MUI_BASE_DIR)
  set(MUI_FORTRAN_LIB_DIR "${MUI_BASE_DIR}/lib")
  # Add the include directories
  include_directories(${MUI_BASE_DIR}/lib/mod)
  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_BASE_DIR}/include/wrappers/Fortran")
  include_directories(${MUI_INCLUDE_DIR})
endif()

find_library(MUI_FORTRAN_WRAPPER_LIB NAMES MUI_Fortran_wrapper PATHS ${MUI_FORTRAN_LIB_DIR})

if(NOT MUI_FORTRAN_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_Fortran_wrapper.a in ${MUI_FORTRAN_LIB_DIR}")
endif()

include_directories(${MUI_FORTRAN_INCLUDE_DIR})

add_executable (heat-left heat-left.cpp)
target_link_libraries (heat-left ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

add_executable (heat-right
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_1d.f90
				heat-right.f90)
target_link_libraries(heat-right
				${MUI_FORTRAN_WRAPPER_LIB}
				stdc++
				mpi_cxx)