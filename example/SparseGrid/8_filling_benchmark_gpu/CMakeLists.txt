cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(sparse_grid_filling_benchmark LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(gray_scott_sparse_gpu SOURCES main.cu)
target_link_libraries(gray_scott_sparse_gpu PRIVATE openfpm::binary_config)
