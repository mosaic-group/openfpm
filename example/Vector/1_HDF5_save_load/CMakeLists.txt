cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_hdf5 LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(hdf5 main.cpp)
target_link_libraries(hdf5 PRIVATE openfpm::binary_config)
