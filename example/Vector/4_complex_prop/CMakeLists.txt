cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vector_complex_properties LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(vect_cp main.cpp)
target_link_libraries(vect_cp PRIVATE openfpm::binary_config)

add_executable(vect_ser main_ser.cpp)
target_link_libraries(vect_ser PRIVATE openfpm::binary_config)
