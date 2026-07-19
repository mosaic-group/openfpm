cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(odeint_sph_dlb LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(sph_dlb main.cpp)
target_link_libraries(sph_dlb PRIVATE openfpm::numerics)

add_executable(sph_dlb_odeint main2.cpp)
target_link_libraries(sph_dlb_odeint PRIVATE openfpm::numerics)
