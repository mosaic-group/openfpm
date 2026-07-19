cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(particle_level_set LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(particle_level_set main.cpp)
target_link_libraries(particle_level_set PRIVATE openfpm::numerics)
