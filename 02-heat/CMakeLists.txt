cmake_minimum_required(VERSION 3.0)
project(ping_pong VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 01" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (heat heat.cpp)
add_executable (heat-coarse heat-coarse.cpp)
add_executable (heat-fine heat-fine.cpp)
target_link_libraries(heat-coarse ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
target_link_libraries(heat-fine ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
