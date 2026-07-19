cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(simple_graph LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(plot main.cpp)
target_link_libraries(plot PRIVATE openfpm::binary_config)
