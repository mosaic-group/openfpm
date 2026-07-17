cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(vcluster_serial_parallel LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

add_executable(serial_parallel main.cpp)
target_link_libraries(serial_parallel PRIVATE openfpm::binary_config)
