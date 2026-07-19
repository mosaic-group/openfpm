cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(pse_derivative_1d_multiprecision LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

option(OPENFPM_BUILD_QUADMATH_EXAMPLE
    "Build the optional libquadmath multiprecision example" OFF)
if (OPENFPM_BUILD_QUADMATH_EXAMPLE)
    find_library(QUADMATH_LIBRARY quadmath)
    if (NOT QUADMATH_LIBRARY)
        message(FATAL_ERROR "OPENFPM_BUILD_QUADMATH_EXAMPLE requires libquadmath")
    endif()
    add_executable(pse_1d_128 main_float128.cpp)
    target_link_libraries(pse_1d_128 PRIVATE
        openfpm::numerics ${QUADMATH_LIBRARY})
endif()
