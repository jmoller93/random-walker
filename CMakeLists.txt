cmake_minimum_required(VERSION 2.8.12)
project(walkers)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include/pybind11/include")
set(SOURCES src/walker.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

add_subdirectory(include/pybind11)
pybind11_add_module(walkers ${SOURCES} src/bindings.cpp)
