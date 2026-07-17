cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vcluster_simple LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(vcluster main.cpp)
target_link_libraries(vcluster PRIVATE openfpm::binary_config)
