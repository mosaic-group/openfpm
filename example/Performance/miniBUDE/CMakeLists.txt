cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(miniBUDE LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(miniBUDE SOURCES main.cu)
target_link_libraries(miniBUDE PRIVATE openfpm::binary_config)
