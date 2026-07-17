cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(pse_diffusion_1d LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(diff_1d main.cpp)
target_link_libraries(diff_1d PRIVATE openfpm::numerics)
