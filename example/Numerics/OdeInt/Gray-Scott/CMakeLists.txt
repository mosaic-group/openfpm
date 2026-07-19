cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(gray_scott_odeint LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(gs_odeint_cuda SOURCES gs-odeint-cuda.cu)
target_link_libraries(gs_odeint_cuda PRIVATE openfpm::numerics)

openfpm_add_gpu_executable(rk4_odeint_cuda SOURCES rk4-odeint-cuda.cu)
target_link_libraries(rk4_odeint_cuda PRIVATE openfpm::numerics)
