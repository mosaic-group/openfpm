cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(memory_bandwidth LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(memBW SOURCES main.cu)
target_link_libraries(memBW PRIVATE openfpm::binary_config)
