cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(diffusion_caco3 LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(diff_CaCO3 SOURCES main.cu)
target_link_libraries(diff_CaCO3 PRIVATE openfpm::numerics)
