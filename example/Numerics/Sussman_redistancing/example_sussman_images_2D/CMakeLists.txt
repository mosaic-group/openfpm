cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sussman_redistancing_images_2d LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(example_sussman_images main.cpp)
target_link_libraries(example_sussman_images PRIVATE openfpm::numerics)
