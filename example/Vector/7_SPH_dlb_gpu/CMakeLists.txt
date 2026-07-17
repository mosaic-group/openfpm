cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sph_dlb_gpu LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(sph_dlb_gpu SOURCES main.cu)
target_link_libraries(sph_dlb_gpu PRIVATE openfpm::numerics)
