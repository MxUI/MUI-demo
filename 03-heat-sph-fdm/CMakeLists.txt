cmake_minimum_required(VERSION 3.0)
project(heat-sph-fdm VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 03" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (sph-standalone sph-standalone.cpp)
add_executable (fdm fdm.cpp)
add_executable (sph sph.cpp)
target_link_libraries(fdm ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
target_link_libraries(sph ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
