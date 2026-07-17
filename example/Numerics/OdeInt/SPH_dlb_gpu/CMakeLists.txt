cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(odeint_sph_dlb_gpu LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

openfpm_add_gpu_executable(sph_dlb SOURCES main.cu)
target_link_libraries(sph_dlb PRIVATE openfpm::numerics)

openfpm_add_gpu_executable(sph_dlb_odeint SOURCES main2.cu)
target_link_libraries(sph_dlb_odeint PRIVATE openfpm::numerics)
