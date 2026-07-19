cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(multiphase_celllist_verlet LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(multip main.cpp)
target_link_libraries(multip PRIVATE openfpm::binary_config)
