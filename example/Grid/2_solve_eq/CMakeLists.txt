cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(grid_solve_equation LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(periodic main.cpp)
target_link_libraries(periodic PRIVATE openfpm::binary_config)
