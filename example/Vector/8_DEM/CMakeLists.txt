cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(discrete_element_method LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(dem main.cpp)
target_link_libraries(dem PRIVATE openfpm::numerics)
