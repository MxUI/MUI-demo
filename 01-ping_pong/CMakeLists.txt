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

add_executable (mui-ping mui-ping.cpp)
add_executable (mui-pong mui-pong.cpp)
target_link_libraries(mui-ping ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
target_link_libraries(mui-pong ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
