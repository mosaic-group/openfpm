cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(syncthreads_benchmark LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(sync_dlb SOURCES main.cu)
target_link_libraries(sync_dlb PRIVATE openfpm::numerics)
