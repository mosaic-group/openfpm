cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(ps_cma_es LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(ps_cma_es main.cpp)
target_link_libraries(ps_cma_es PRIVATE openfpm::numerics)
target_include_directories(ps_cma_es PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
