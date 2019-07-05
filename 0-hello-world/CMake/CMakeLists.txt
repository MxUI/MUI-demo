cmake_minimum_required(VERSION 3.9)
project(hello VERSION 1.0 DESCRIPTION "Multiscale Universal Interface Demo 0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(mui REQUIRED)

add_executable (hello ./hello.cpp )
target_link_libraries (hello ${MPI_LIBRARIES} mui)
install(TARGETS hello DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
