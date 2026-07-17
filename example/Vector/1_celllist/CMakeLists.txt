cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_celllist LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(cell main.cpp)
target_link_libraries(cell PRIVATE openfpm::binary_config)
