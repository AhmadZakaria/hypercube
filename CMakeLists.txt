project(TNN-G)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg  -O3 -fopenmp -march=native -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -Wall  -debug")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -Wall")

cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

