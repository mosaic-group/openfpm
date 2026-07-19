cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(heat_conduction_ceramics LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(heat_ceram SOURCES main.cu)
target_link_libraries(heat_ceram PRIVATE openfpm::numerics)
