cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sparse_grid_filling_benchmark LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(gray_scott_sparse_opt main.cpp)
target_link_libraries(gray_scott_sparse_opt PRIVATE openfpm::binary_config)
