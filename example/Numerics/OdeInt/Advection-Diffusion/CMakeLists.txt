cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(odeint_advection_diffusion LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(example_odeint main.cpp)
target_link_libraries(example_odeint PRIVATE openfpm::numerics)

add_executable(example_odeint2 main2.cpp)
target_link_libraries(example_odeint2 PRIVATE openfpm::numerics)
