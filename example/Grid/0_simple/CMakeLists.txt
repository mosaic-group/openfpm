cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(grid_simple LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(grid main.cpp)
target_link_libraries(grid PRIVATE openfpm::binary_config)
