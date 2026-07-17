cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(molecular_dynamics_gpu LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(md_dyn SOURCES main.cu)
target_link_libraries(md_dyn PRIVATE openfpm::binary_config)
