cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(molecular_dynamics_symmetric_crs LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(md_sym main.cpp)
target_link_libraries(md_sym PRIVATE openfpm::binary_config)
