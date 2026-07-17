cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sphere_redistancing LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(sphere_redistance main.cpp)
target_link_libraries(sphere_redistance PRIVATE openfpm::numerics)
