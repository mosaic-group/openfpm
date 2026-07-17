cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(surface_dcpse_circle LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(Circ CircLap.cpp)
target_link_libraries(Circ PRIVATE openfpm::numerics)

add_executable(CircPoisson CircPoisson.cpp)
target_link_libraries(CircPoisson PRIVATE openfpm::numerics)
