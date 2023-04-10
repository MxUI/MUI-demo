cmake_minimum_required(VERSION 3.9)
project(smart-send-dual-times-Fortran VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-1-2" LANGUAGES Fortran)

# Set the Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_COMPILER mpifort)

set(MUI_BASE_FORTRAN_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_FORTRAN_LIB_DIR "" CACHE PATH "Path to MUI FORTRAN library directory")
set(MUI_FORTRAN_INCLUDE_DIR "" CACHE PATH "Path to MUI FORTRAN include directory")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

find_package(MUI QUIET)

if(MUI_BASE_DIR)
  set(MUI_FORTRAN_LIB_DIR "${MUI_LIB_DIR}")
  # Add the include directories
  include_directories(${MUI_LIB_DIR}/mod)
  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_INCLUDE_DIR}/wrappers/Fortran")
  include_directories(${MUI_INCLUDE_DIR})
elseif(MUI_BASE_FORTRAN_DIR)
  set(MUI_FORTRAN_LIB_DIR "${MUI_BASE_FORTRAN_DIR}/lib")
  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_BASE_FORTRAN_DIR}/include/wrappers/Fortran")
endif()

find_library(MUI_FORTRAN_WRAPPER_LIB NAMES MUI_Fortran_wrapper PATHS ${MUI_FORTRAN_LIB_DIR})

if(NOT MUI_FORTRAN_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_Fortran_wrapper.a in ${MUI_FORTRAN_LIB_DIR}")
endif()

include_directories(${MUI_FORTRAN_INCLUDE_DIR})

add_executable (code1
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_1d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_2d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_3d.f90 
				code1.f90)

# Add the libraries to link
target_link_libraries(code1
				${MUI_FORTRAN_WRAPPER_LIB}
				stdc++
				mpi_cxx)

add_executable (code2
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_1d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_2d.f90 
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_3d.f90 
				code2.f90)

# Add the libraries to link
target_link_libraries(code2
				${MUI_FORTRAN_WRAPPER_LIB}
				stdc++
				mpi_cxx)
