cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sussman_redistancing_sphere LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(example_sussman_sphere main.cpp)
target_link_libraries(example_sussman_sphere PRIVATE openfpm::numerics)
