cmake_minimum_required(VERSION 3.0)
project(2d-exchange VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 06" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (domain-coarse domain-coarse.cpp)
add_executable (domain-refine domain-refine.cpp)
target_link_libraries(domain-coarse ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
target_link_libraries(domain-refine ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
