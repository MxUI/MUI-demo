cmake_minimum_required(VERSION 3.9)
project(Python_Wrapper-dual-time-types VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 10-2-5" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (PUSHER_FETCHER_0 PUSHER_FETCHER_0.cpp)
target_link_libraries (PUSHER_FETCHER_0 ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})