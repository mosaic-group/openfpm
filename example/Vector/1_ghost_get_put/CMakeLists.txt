cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_ghost_get_put LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(ghost main.cpp)
target_link_libraries(ghost PRIVATE openfpm::binary_config)
