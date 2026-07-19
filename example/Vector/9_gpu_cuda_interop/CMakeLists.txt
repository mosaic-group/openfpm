cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(gpu_cuda_interop LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(gpu_interop SOURCES main.cu)
target_link_libraries(gpu_interop PRIVATE openfpm::binary_config)
