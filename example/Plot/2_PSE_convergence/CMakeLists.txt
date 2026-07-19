cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(pse_convergence_plot LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)

option(OPENFPM_BUILD_QUADMATH_EXAMPLE
    "Build the optional libquadmath convergence example" OFF)
if (OPENFPM_BUILD_QUADMATH_EXAMPLE)
    find_library(QUADMATH_LIBRARY quadmath)
    if (NOT QUADMATH_LIBRARY)
        message(FATAL_ERROR "OPENFPM_BUILD_QUADMATH_EXAMPLE requires libquadmath")
    endif()
    add_executable(conv_p main.cpp)
    target_link_libraries(conv_p PRIVATE
        openfpm::numerics ${QUADMATH_LIBRARY})
endif()
