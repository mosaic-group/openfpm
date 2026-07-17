cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(stokes_2d_incompressible LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(stokes_2d_eigen main_eigen.cpp)
target_link_libraries(stokes_2d_eigen PRIVATE openfpm::numerics)

add_executable(stokes_2d_petsc main_petsc.cpp)
target_link_libraries(stokes_2d_petsc PRIVATE openfpm::numerics)
