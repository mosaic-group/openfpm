cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_expressions LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(expr main.cpp)
target_link_libraries(expr PRIVATE openfpm::numerics)
