cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(nbody_benchmark LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(nbody SOURCES main.cu)
target_link_libraries(nbody PRIVATE openfpm::binary_config)
