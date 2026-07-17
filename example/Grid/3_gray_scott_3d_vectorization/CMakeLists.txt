cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(grid_gray_scott_3d_vectorization LANGUAGES CXX Fortran)

find_package(openfpm CONFIG REQUIRED)

add_executable(gray_scott main.cpp update_new.f90)
target_link_libraries(gray_scott PRIVATE openfpm::binary_config)
target_compile_options(gray_scott PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:
        -ffree-line-length-none -fno-range-check -fno-second-underscore
        -fimplicit-none>)
