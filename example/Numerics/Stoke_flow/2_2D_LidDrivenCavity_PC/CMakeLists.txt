cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(lid_driven_cavity LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(DcpseLid mainDCPSE.cpp)
target_link_libraries(DcpseLid PRIVATE openfpm::numerics)

add_executable(FDLid mainFD.cpp)
target_link_libraries(FDLid PRIVATE openfpm::numerics)
