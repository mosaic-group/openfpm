cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(gpu_first_step LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(gpu_fstep SOURCES main.cu)
target_link_libraries(gpu_fstep PRIVATE openfpm::binary_config)
