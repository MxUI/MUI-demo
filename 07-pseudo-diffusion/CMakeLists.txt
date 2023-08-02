cmake_minimum_required(VERSION 3.0)
project(pseudo-diffusion VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 07" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (3D-pseudo-diffusion-coarse 3D-pseudo-diffusion-coarse.cpp)
add_executable (3D-pseudo-diffusion-fine 3D-pseudo-diffusion-fine.cpp)
target_link_libraries(3D-pseudo-diffusion-coarse ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
target_link_libraries(3D-pseudo-diffusion-fine ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
