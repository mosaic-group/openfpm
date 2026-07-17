cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_simple LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(vect main.cpp)
target_link_libraries(vect PRIVATE openfpm::binary_config)
