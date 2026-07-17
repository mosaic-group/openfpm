cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vortex_in_cell LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(vic_petsc main_vic_petsc.cpp)
target_link_libraries(vic_petsc PRIVATE openfpm::numerics)

add_executable(vic_petsc_opt main_vic_petsc_opt.cpp)
target_link_libraries(vic_petsc_opt PRIVATE openfpm::numerics)
