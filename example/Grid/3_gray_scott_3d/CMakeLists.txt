cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(grid_gray_scott_3d LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(gray_scott main.cpp)
target_link_libraries(gray_scott PRIVATE openfpm::binary_config)
