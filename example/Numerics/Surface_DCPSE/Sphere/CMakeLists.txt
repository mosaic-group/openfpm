cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(surface_dcpse_sphere LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(SphereE Sph.cpp)
target_link_libraries(SphereE PRIVATE openfpm::numerics)

add_executable(SphereBench SphBench.cpp)
target_link_libraries(SphereBench PRIVATE openfpm::numerics)
