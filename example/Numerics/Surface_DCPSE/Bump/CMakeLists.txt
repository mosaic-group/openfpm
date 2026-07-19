cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(surface_dcpse_bump LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(Bump Bump.cpp)
target_link_libraries(Bump PRIVATE openfpm::numerics)
