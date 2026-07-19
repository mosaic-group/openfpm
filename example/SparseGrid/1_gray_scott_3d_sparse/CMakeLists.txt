cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(gray_scott_sparse LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(gray_scott_sparse main.cpp)
target_link_libraries(gray_scott_sparse PRIVATE openfpm::binary_config)
