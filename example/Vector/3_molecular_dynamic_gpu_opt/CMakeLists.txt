cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(molecular_dynamics_gpu_opt LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(md_dyn_gpu SOURCES main_gpu.cu)
target_link_libraries(md_dyn_gpu PRIVATE openfpm::binary_config)

add_executable(md_dyn_cpu main_cpu.cpp)
target_link_libraries(md_dyn_cpu PRIVATE openfpm::binary_config)

add_executable(md_dyn_cpu_best main_cpu_best.cpp)
target_link_libraries(md_dyn_cpu_best PRIVATE openfpm::binary_config)
