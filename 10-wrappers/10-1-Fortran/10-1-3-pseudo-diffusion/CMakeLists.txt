cmake_minimum_required(VERSION 3.9)
project(pseudo-diffusion-Fortran VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-1-3" LANGUAGES Fortran CXX)

# Set the Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_COMPILER mpifort)
set(MUI_BASE_FORTRAN_DIR "" CACHE PATH "Path to MUI base directory")
set(MUI_FORTRAN_LIB_DIR "" CACHE PATH "Path to MUI FORTRAN library directory")
set(MUI_FORTRAN_INCLUDE_DIR "" CACHE PATH "Path to MUI FORTRAN include directory")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MUI_FORTRAN_CPP_SEPERATE_PATH 0)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

find_package(MUI REQUIRED)

if(${MUI_FORTRAN_CPP_SEPERATE_PATH} EQUAL 1)
	if(MUI_BASE_DIR)
	  include_directories(${MUI_INCLUDE_DIR})
	elseif(MUI_BASE_FORTRAN_DIR)
	  set(MUI_FORTRAN_LIB_DIR "${MUI_BASE_FORTRAN_DIR}/lib")
	  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_BASE_FORTRAN_DIR}/include/wrappers/Fortran")
	  include_directories(${MUI_BASE_FORTRAN_DIR})
	endif()
else()
	if(MUI_BASE_DIR)
	  set(MUI_FORTRAN_LIB_DIR "${MUI_LIB_DIR}")
	  # Add the include directories
	  include_directories(${MUI_LIB_DIR}/mod)
	  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_INCLUDE_DIR}/wrappers/Fortran")
	  include_directories(${MUI_INCLUDE_DIR})
	elseif(MUI_BASE_FORTRAN_DIR)
	  set(MUI_FORTRAN_LIB_DIR "${MUI_BASE_FORTRAN_DIR}/lib")
	  set(MUI_FORTRAN_INCLUDE_DIR "${MUI_BASE_FORTRAN_DIR}/include/wrappers/Fortran")
	  include_directories(${MUI_BASE_FORTRAN_DIR})
	endif()
endif()

find_library(MUI_FORTRAN_WRAPPER_LIB NAMES MUI_Fortran_wrapper PATHS ${MUI_FORTRAN_LIB_DIR})

if(NOT MUI_FORTRAN_WRAPPER_LIB)
  message(FATAL_ERROR "Cannot find libMUI_Fortran_wrapper.a in ${MUI_FORTRAN_LIB_DIR}")
endif()

include_directories(${MUI_FORTRAN_INCLUDE_DIR})

add_executable (3D-pseudo-diffusion-fine 3D-pseudo-diffusion-fine.cpp)
target_link_libraries (3D-pseudo-diffusion-fine ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

add_executable (3D-pseudo-diffusion-coarse
				${MUI_FORTRAN_INCLUDE_DIR}/mui_f_wrapper_2d.f90
				3D-pseudo-diffusion-coarse.f90)
target_link_libraries(3D-pseudo-diffusion-coarse
				${MUI_FORTRAN_WRAPPER_LIB}
				stdc++
				mpi_cxx)