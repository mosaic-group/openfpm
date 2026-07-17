cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(grid_stencil LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(stencil main.cpp)
target_link_libraries(stencil PRIVATE openfpm::binary_config)
