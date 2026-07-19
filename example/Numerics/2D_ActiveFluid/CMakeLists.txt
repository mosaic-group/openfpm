cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(active_fluid_2d LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(Active2d main.cpp)
target_link_libraries(Active2d PRIVATE openfpm::numerics)
