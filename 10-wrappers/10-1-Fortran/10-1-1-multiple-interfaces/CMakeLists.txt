cmake_minimum_required(VERSION 3.9)
project(multi-interfaces-Fortran VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-1-1" LANGUAGES Fortran)

# Set the Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_COMPILER mpifort)

set(MUI_BASE_FORTRAN_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_FORTRAN_LIB_DIR "" CACHE PATH "Path to MUI FORTRAN library directory")
set(MUI_FORTRAN_INCLUDE_DIR "" CACHE PATH "Path to MUI FORTRAN include directory")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MUI QUIET)

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

add_executable (multi-interfaces 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_1d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_2d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_3d.f90 
				multi-interfaces.f90)

# Add the libraries to link
target_link_libraries(multi-interfaces
				${MUI_FORTRAN_WRAPPER_LIB}
				stdc++
				mpi_cxx)
