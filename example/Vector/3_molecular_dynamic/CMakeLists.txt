cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(molecular_dynamics LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

add_executable(md_dyn main.cpp)
target_link_libraries(md_dyn PRIVATE openfpm::numerics)

add_executable(md_dyn_expr main_expr.cpp)
target_link_libraries(md_dyn_expr PRIVATE openfpm::numerics)

add_executable(md_dyn_vl main_vl.cpp)
target_link_libraries(md_dyn_vl PRIVATE openfpm::numerics)

add_executable(md_dyn_pap main_expr_paper.cpp)
target_link_libraries(md_dyn_pap PRIVATE openfpm::numerics)
