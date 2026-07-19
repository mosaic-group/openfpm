cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_complex_usage LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(complex_use main.cpp)
target_link_libraries(complex_use PRIVATE openfpm::binary_config)
