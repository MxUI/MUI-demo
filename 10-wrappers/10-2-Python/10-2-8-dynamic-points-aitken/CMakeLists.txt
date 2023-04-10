cmake_minimum_required(VERSION 3.9)
project(python_aitken_dyanmic_points VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-2-8" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (heat-left heat-left.cpp)
target_link_libraries (heat-left ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})