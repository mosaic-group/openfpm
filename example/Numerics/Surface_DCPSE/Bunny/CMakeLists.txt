cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(surface_dcpse_bunny LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(Bunny Bunny.cpp)
target_link_libraries(Bunny PRIVATE openfpm::numerics)
