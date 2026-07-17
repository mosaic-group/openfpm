cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(stokes_flow_ball_3d LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(Stokes3dBall main.cpp)
target_link_libraries(Stokes3dBall PRIVATE openfpm::numerics)
