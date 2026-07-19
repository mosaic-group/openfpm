cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_reorder LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(md_data_ord main_data_ord.cpp)
target_link_libraries(md_data_ord PRIVATE openfpm::binary_config)

add_executable(md_comp_ord main_comp_ord.cpp)
target_link_libraries(md_comp_ord PRIVATE openfpm::binary_config)
